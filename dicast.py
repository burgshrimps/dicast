import argparse
import logging
import sys
from collections import OrderedDict
from joblib import Parallel, delayed
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import re
import resource
from contextlib import contextmanager
from time import perf_counter
from datetime import datetime
import vcfpy
import pysam
import pickle


@contextmanager
def stage_timer(stage_name, rows):
    """Record wall-time, self+children CPU-time and peak RSS delta for a stage.

    `rows` is a list that receives dicts; callers aggregate them at the end."""
    ru0_self = resource.getrusage(resource.RUSAGE_SELF)
    ru0_children = resource.getrusage(resource.RUSAGE_CHILDREN)
    t0 = perf_counter()
    try:
        yield
    finally:
        wall = perf_counter() - t0
        ru1_self = resource.getrusage(resource.RUSAGE_SELF)
        ru1_children = resource.getrusage(resource.RUSAGE_CHILDREN)
        cpu_self = (ru1_self.ru_utime + ru1_self.ru_stime) - (ru0_self.ru_utime + ru0_self.ru_stime)
        cpu_children = (ru1_children.ru_utime + ru1_children.ru_stime) - (ru0_children.ru_utime + ru0_children.ru_stime)
        # ru_maxrss is a high-water mark (not a delta); Linux reports kB.
        peak_rss_kb = max(ru1_self.ru_maxrss, ru1_children.ru_maxrss)
        rows.append({
            'stage': stage_name,
            'wall_seconds': round(wall, 3),
            'cpu_self_seconds': round(cpu_self, 3),
            'cpu_children_seconds': round(cpu_children, 3),
            'peak_rss_gb': round(peak_rss_kb / (1024 * 1024), 3),
        })

from dicast_lib.parsing import parse_arguments, validate_inputs
from dicast_lib.utils import replace_filename
from dicast_lib.prepare import VariantPrep
from dicast_lib.collect_reference import ReferenceAnnotator
from dicast_lib.collect_illumina import AlignmentAnnotatorIllumina
from dicast_lib.model import Dicast
from dicast_lib.multi import find_rescue_candidates

# Suppresses stacked htslib warnings when a BAM/CRAM index predates the alignment file.
pysam.set_verbosity(0)


# List of chromosomes to process
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
          'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
          'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

# List of SV types currently supported by dicast
sv_types = ['DEL', 'DUP', 'INS', 'INV']


def combine_feature_files(sample: str, ref: str, workdir: str) -> pd.DataFrame:
    """ Combines raw, reference and alignment features into a single TSV file.

    Args:
        sample (str): Sample name
        ref (str): Reference genome name
        workdir (str): Working directory

    Returns:
        pd.DataFrame: DataFrame with combined features
    """    

    df_raw = pd.read_csv(f'{workdir}/{sample}_{ref}.SVs.raw.tsv', sep='\t', low_memory=False,
                         dtype={'sample': str, 'cohort_samples': str})
    df_ref = pd.read_csv(f'{workdir}/{sample}_{ref}.SVs.ref.tsv', sep='\t', low_memory=False)
    
    filenames_aln_ill = glob(f'{workdir}/{sample}_{ref}.SVs.aln.ill.*.*.tsv')
    
    df_aln_ill = pd.concat([pd.read_csv(f, sep='\t') for f in filenames_aln_ill], ignore_index=True)
    
    df = df_raw.merge(df_ref.drop(['sample', 'sv_type', 'chrom', 'chrom_2', 'start', 'end', 'cohort', 'technology', 'caller', 'reference'], axis=1), on='id', how='inner')
    df = df.merge(df_aln_ill.drop(['sample', 'sv_type', 'chrom', 'chrom_2', 'start', 'end', 'sv_len', 'cohort', 'technology', 'caller', 'reference'], axis=1), on='id', how='inner')

    return df


def add_info_tag_to_vcf(arguments: argparse.Namespace):
    """ Adds dicast quality to input VCFs.

    Args:
        arguments (argparse.Namespace): Parsed command line arguments
    """    

    dicast_df = pd.read_csv(f'{arguments.workdir}/{arguments.sample}_{arguments.ref}.SVs.dicast.tsv', sep='\t',
                            dtype={'sample': str, 'cohort_samples': str})

    for caller, vcf_filename in arguments.vcfs:
        vcf_in = vcfpy.Reader.from_path(vcf_filename)
        vcf_in.header.add_info_line(OrderedDict([('ID', 'DQ'),
                                                 ('Number', '1'),
                                                 ('Type', 'String'),
                                                 ('Description', 'Dicast Quality Score')]))
        vcf_basename_out = os.path.basename(vcf_filename).replace('.vcf', '.dicast.vcf').replace('.gz', '')
        vcf_filename_out = os.path.join(arguments.workdir, vcf_basename_out)
        vcf_out = vcfpy.Writer.from_path(vcf_filename_out, vcf_in.header)

        dicast_df_caller = dicast_df[dicast_df['caller']== caller].copy().reset_index(drop=True)
        dicast_df_caller = dicast_df_caller[['id', 'dicast_qual']].set_index('id').T.to_dict('list')

        for rec in vcf_in:
            if rec.ID[0] in dicast_df_caller.keys():
                qual_dicast = dicast_df_caller[rec.ID[0]][0]
            else:
                qual_dicast = -1

            rec.INFO['DQ'] = str(qual_dicast)
            vcf_out.write_record(rec)
            
        logging.info(
            f'Added DQ tag to {caller} VCF file')


def extract_reference_features(arguments: argparse.Namespace, sample: str):
    """ Extracts reference features for a given sample.

    Args:
        arguments (argparse.Namespace): Parsed command line arguments
        sample (str): Sample name
    """    
    
    reference_filenames = {
            'repeats_filename' : arguments.repeats,
            'vntrs_filename' : arguments.vntrs,
            'strs_filename' : arguments.strs,
            'cpgislands_filename' : arguments.cgis,
            'centromeres_filename' : arguments.centromeres,
            'asmb_gaps_filename' : arguments.gaps,
            'alt_haps_filename' : arguments.althaps,
            'gc_filename' : arguments.gc
    }
    
    RA = ReferenceAnnotator(reference_filenames)
    RA.load_from_csv('/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.raw.tsv']))
    RA.split_bnd()
    logging.info('# Annotation Repeats')
    RA.annotate_repeats()
    logging.info('# Annotation VNTRs')
    RA.annotate_vntrs()
    logging.info('# Annotation STRs')
    RA.annotate_strs()
    logging.info('# Annotation CGIs')
    RA.annotate_cpg_islands()
    logging.info('# Annotation Centromeres')
    RA.annotate_centromeres()
    logging.info('# Annotation Assembly Gaps')
    RA.annotate_asmb_gaps()
    logging.info('# Annotation Alternative Haplotypes')
    RA.annotate_alt_haps()
    logging.info('# Annotation GC Content')
    RA.annotate_gc_content()
    logging.info('# Aggregate Results')
    RA.aggregate_results()
    logging.info('# Save Reference Features')
    RA.to_csv('/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.ref.tsv']))


def collect_aln_features(bam_filename: str, variant_filename: str, variant_annot_filename: str, chrom: str, sv_type: str, sample: str):
    """ Collects alignment features for a given chromosome.

    Args:
        bam_filename (str): BAM file
        variant_filename (str):  TSV file with variants to annotate
        variant_annot_filename (str): Output TSV file with annotated variants
        chrom (str): Chromosome name
        sv_type (str): SV type
        sample (str): Sample name
    """

    AAI = AlignmentAnnotatorIllumina(bam_filename, chrom, sv_type, sample)
    AAI.load_from_csv(variant_filename)
    AAI.calculate_coverage_baseline()
    AAI.calculate_insertsize_baseline()
    AAI.calculate_mapping_quality_baseline()
    AAI.annotate_coverage()
    AAI.annotate_read_based_features()
    AAI.to_csv(variant_annot_filename)


def score_variants(sv_types: list, arguments: argparse.Namespace, sample: str):
    """ Uses dicast model to score variants.

    Args:
        sv_types (list): List of SV types
        arguments (argparse.Namespace): Parsed command line arguments
        sample (str): Sample name
    """    

    variant_features_filename = f'{arguments.workdir}/{sample}_{arguments.ref}.SVs.annot.tsv'
    use_pop_models = getattr(arguments, 'pop', False)
    dicast_dfs = []
    for sv_type in sv_types:
        if sv_type != 'INV':
            model_filename = f'{arguments.models}/dicast_{sv_type}.json'
            if use_pop_models and sv_type in ('DEL', 'INS'):
                pop_model_filename = f'{arguments.models}/dicast_{sv_type}_pop.json'
                if os.path.isfile(pop_model_filename):
                    model_filename = pop_model_filename
                else:
                    logging.info(f'--pop set but no population model for {sv_type} (pop models ship for DEL and INS only); falling back to {model_filename}')
            dicast = Dicast(sv_type)
            dicast.load_from_csv(variant_features_filename)
            dicast.impute_missing_values()
            dicast.load(model_filename)
            dicast.predict()
            dicast_dfs.append(dicast.to_df())   
        else:
            dicast = Dicast(sv_type)
            dicast.load_from_csv(variant_features_filename)
            dicast.impute_missing_values()
            dicast.score_inversions()
            dicast_dfs.append(dicast.to_df())
    dicast_df = pd.concat([df for df in dicast_dfs if not df.empty], ignore_index=True)
    dicast_df.to_csv(f'{arguments.workdir}/{sample}_{arguments.ref}.SVs.dicast.tsv', sep='\t', index=False, na_rep='NA')


if __name__ == '__main__':

    # Parse command line arguments
    arguments = parse_arguments()
    validate_inputs(arguments)

    # Set up logging
    logFormatter = logging.Formatter('%(asctime)s %(message)s')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    print('')
    logging.info('############### Start DICAST ###############\n')
    logging.info('CMD: python3 {0}'.format(' '.join(sys.argv)))
    print('')

    if arguments.command == 'call':

        logging.info('MODE: call')
        logging.info(f'COHORT: {arguments.cohort}')
        logging.info(f'SAMPLE: {arguments.sample}')
        logging.info(f'CHROMOSOMES: {arguments.chrom}')
        logging.info(f'REF: {arguments.ref}')
        logging.info(f'TECH: {arguments.technology}')
        logging.info(f'WORKDIR: {arguments.workdir}')
        logging.info(f'FAI: {arguments.fai}')
        logging.info(f'ANNOT-DIR: {arguments.annot_dir}')
        logging.info(f'REPEATS: {arguments.repeats}')
        logging.info(f'CGIS: {arguments.cgis}')
        logging.info(f'CENTROMERES: {arguments.centromeres}')
        logging.info(f'GAPS: {arguments.gaps}')
        logging.info(f'ALTHAPS: {arguments.althaps}')
        logging.info(f'VNTRS: {arguments.vntrs}')
        logging.info(f'STRS: {arguments.strs}')
        logging.info(f'GC: {arguments.gc}')
        logging.info(f'BAM: {arguments.bam}')
        logging.info(f'THREADS: {arguments.threads}')
        logging.info(f'MODELS: {arguments.models}')
        logging.info(f'POP: {arguments.pop}')
        logging.info(f'SV CALLERS: {", ".join([caller for caller, _ in arguments.vcfs])}')
        logging.info(f'VCFs: {", ".join([vcf for _, vcf in arguments.vcfs])}')
        print('')

        # Restrict feature extraction to a single chromosome if specified
        if arguments.chrom != 'all':
            chroms = arguments.chrom

        # Restrict SV types if specified (default: use the module-level sv_types)
        if arguments.sv_types:
            sv_types = arguments.sv_types
            logging.info(f'SV TYPES (restricted): {sv_types}')

        benchmark_rows = []
        total_t0 = perf_counter()

        def _write_benchmark():
            if not arguments.benchmark:
                return
            total_wall = round(perf_counter() - total_t0, 3)
            ru_self = resource.getrusage(resource.RUSAGE_SELF)
            ru_children = resource.getrusage(resource.RUSAGE_CHILDREN)
            total_cpu_self = round(ru_self.ru_utime + ru_self.ru_stime, 3)
            total_cpu_children = round(ru_children.ru_utime + ru_children.ru_stime, 3)
            total_peak_rss_gb = round(max(ru_self.ru_maxrss, ru_children.ru_maxrss) / (1024 * 1024), 3)
            rows = list(benchmark_rows) + [{
                'stage': 'total',
                'wall_seconds': total_wall,
                'cpu_self_seconds': total_cpu_self,
                'cpu_children_seconds': total_cpu_children,
                'peak_rss_gb': total_peak_rss_gb,
            }]
            os.makedirs(os.path.dirname(os.path.abspath(arguments.benchmark)), exist_ok=True)
            pd.DataFrame(rows).to_csv(arguments.benchmark, sep='\t', index=False)
            logging.info(f'# Benchmark TSV written to {arguments.benchmark}')

        try:
            with stage_timer('feature_collection', benchmark_rows):
                # Variant Preparation
                logging.info('# Create Variant DataFrame')
                VP = VariantPrep(arguments.cohort, arguments.ref, arguments.workdir,
                                 arguments.technology, chroms, arguments.fai, sv_types)
                logging.info('# Read VCFs')
                VP.read_vcf(arguments.vcfs, arguments.sample)
                logging.info('# Read Variants')
                VP.read_variants()
                logging.info('# Filter Variants')
                VP.filter_variants()
                logging.info('# Save Variants')
                VP.save_variants()

                # Reference Feature Collection
                logging.info('# Collect Reference Features')
                extract_reference_features(arguments, arguments.sample)

                # Alignment Feature Collection
                logging.info(f'# Collect Alignment Features')
                parallel_input = []
                for sv_type in sv_types:
                    for chrom in chroms:
                        variant_filename = '/'.join([arguments.workdir, arguments.sample + '_' + arguments.ref + '.SVs.raw.tsv'])
                        variant_annot_filename = '/'.join([arguments.workdir, arguments.sample + '_' + arguments.ref + '.SVs.aln.ill.' + chrom + '.' + sv_type + '.tsv'])
                        parallel_input.append((arguments.bam, variant_filename, variant_annot_filename,
                                            chrom, sv_type, arguments.sample))
                Parallel(n_jobs=arguments.threads)(delayed(collect_aln_features)(*args) for args in parallel_input)

                # Combination of Raw, Reference and Alignment Features
                logging.info('# Combination Output Files')
                df = combine_feature_files(arguments.sample, arguments.ref, arguments.workdir)
                df.to_csv(f'{arguments.workdir}/{arguments.sample}_{arguments.ref}.SVs.annot.tsv', sep='\t', index=False, na_rep='NA')

            with stage_timer('prediction', benchmark_rows):
                # Variant Prediction
                logging.info('# Variant Prediction')
                score_variants(sv_types, arguments, arguments.sample)

                # Add Info Tag to VCF
                logging.info('# Add Info Tag to VCF')
                add_info_tag_to_vcf(arguments)
        finally:
            _write_benchmark()

        logging.info('############### End DICAST ###############\n')

    elif arguments.command == 'multi':

        logging.info('MODE: multi')
        logging.info(f'COHORT: {arguments.cohort}')
        logging.info(f'CHROMOSOMES: {arguments.chrom}')
        logging.info(f'REF: {arguments.ref}')
        logging.info(f'TECH: {arguments.technology}')
        logging.info(f'WORKDIR: {arguments.workdir}')
        logging.info(f'FAI: {arguments.fai}')
        logging.info(f'ANNOT-DIR: {arguments.annot_dir}')
        logging.info(f'MODELS: {arguments.models}')
        logging.info(f'THREADS: {arguments.threads}')
        logging.info(f'POP: {arguments.pop}')

        bam_dict = {sample: bam for sample, bam in arguments.bams}
        samples = list(bam_dict.keys())
        vcfs_by_sample = {sample: [] for sample in samples}
        for sample, caller, vcf_file in arguments.vcfs:
            vcfs_by_sample[sample].append([caller, vcf_file])

        logging.info(f'SAMPLES: {", ".join(samples)}')
        for sample in samples:
            logging.info(f'BAM ({sample}): {bam_dict[sample]}')
            logging.info(f'SV CALLERS ({sample}): {", ".join(caller for caller, _ in vcfs_by_sample[sample])}')
        print('')

        # Restrict feature extraction to a single chromosome if specified
        if arguments.chrom != 'all':
            chroms = arguments.chrom

        # Restrict SV types if specified (default: use the module-level sv_types)
        if arguments.sv_types:
            sv_types = arguments.sv_types
            logging.info(f'SV TYPES (restricted): {sv_types}')

        benchmark_rows = []
        total_t0 = perf_counter()

        def _write_benchmark():
            if not arguments.benchmark:
                return
            total_wall = round(perf_counter() - total_t0, 3)
            ru_self = resource.getrusage(resource.RUSAGE_SELF)
            ru_children = resource.getrusage(resource.RUSAGE_CHILDREN)
            total_cpu_self = round(ru_self.ru_utime + ru_self.ru_stime, 3)
            total_cpu_children = round(ru_children.ru_utime + ru_children.ru_stime, 3)
            total_peak_rss_gb = round(max(ru_self.ru_maxrss, ru_children.ru_maxrss) / (1024 * 1024), 3)
            rows = list(benchmark_rows) + [{
                'stage': 'total',
                'wall_seconds': total_wall,
                'cpu_self_seconds': total_cpu_self,
                'cpu_children_seconds': total_cpu_children,
                'peak_rss_gb': total_peak_rss_gb,
            }]
            os.makedirs(os.path.dirname(os.path.abspath(arguments.benchmark)), exist_ok=True)
            pd.DataFrame(rows).to_csv(arguments.benchmark, sep='\t', index=False)
            logging.info(f'# Benchmark TSV written to {arguments.benchmark}')

        try:
            with stage_timer('feature_collection', benchmark_rows):
                # Variant Preparation: each sample's own caller variants first
                own_variant_dfs = {}
                for sample in samples:
                    logging.info(f'# Create Variant DataFrame for {sample}')
                    VP = VariantPrep(arguments.cohort, arguments.ref, arguments.workdir,
                                     arguments.technology, chroms, arguments.fai, sv_types)
                    VP.read_vcf(vcfs_by_sample[sample], sample)
                    VP.read_variants()
                    VP.filter_variants()
                    own_variant_dfs[sample] = VP.get_variant_df()

                # Rescue: add variants found by the OTHER samples' callers but missing from
                # this sample's own callers, so they still get scored against this sample's BAM.
                # Rescued rows are tagged via 'caller' (rescue:<origin sample>:<origin caller>)
                # so they never collide with a real caller name in add_info_tag_to_vcf below.
                logging.info('# Determine cross-sample rescue candidates')
                rescue_dfs = find_rescue_candidates(own_variant_dfs)

                for sample in samples:
                    df_sample = pd.concat([own_variant_dfs[sample], rescue_dfs[sample]], ignore_index=True)
                    logging.info(f'# Save Variants for {sample} ({len(own_variant_dfs[sample])} own, {len(rescue_dfs[sample])} rescued)')
                    df_sample.to_csv(f'{arguments.workdir}/{sample}_{arguments.ref}.SVs.raw.tsv', sep='\t', index=False, na_rep='NA')

                # Reference Feature Collection
                for sample in samples:
                    logging.info(f'# Collect Reference Features for {sample}')
                    extract_reference_features(arguments, sample)

                # Alignment Feature Collection
                logging.info('# Collect Alignment Features')
                parallel_input = []
                for sample in samples:
                    for sv_type in sv_types:
                        for chrom in chroms:
                            variant_filename = '/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.raw.tsv'])
                            variant_annot_filename = '/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.aln.ill.' + chrom + '.' + sv_type + '.tsv'])
                            parallel_input.append((bam_dict[sample], variant_filename, variant_annot_filename,
                                                chrom, sv_type, sample))
                Parallel(n_jobs=arguments.threads)(delayed(collect_aln_features)(*args) for args in parallel_input)

                # Combination of Raw, Reference and Alignment Features
                logging.info('# Combination Output Files')
                for sample in samples:
                    df = combine_feature_files(sample, arguments.ref, arguments.workdir)
                    df.to_csv(f'{arguments.workdir}/{sample}_{arguments.ref}.SVs.annot.tsv', sep='\t', index=False, na_rep='NA')

            with stage_timer('prediction', benchmark_rows):
                for sample in samples:
                    # Variant Prediction
                    logging.info(f'# Variant Prediction for {sample}')
                    score_variants(sv_types, arguments, sample)

                    # Add Info Tag to VCF (only tags this sample's own caller VCFs)
                    logging.info(f'# Add Info Tag to VCF for {sample}')
                    arguments.sample = sample
                    arguments.vcfs = vcfs_by_sample[sample]
                    add_info_tag_to_vcf(arguments)
        finally:
            _write_benchmark()

        logging.info('############### End DICAST ###############\n')