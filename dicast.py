import argparse
import logging
import sys
from joblib import Parallel, delayed
import pandas as pd
from glob import glob
from tqdm import tqdm
import os
import re
from datetime import datetime
import vcfpy
import pysam
import pickle
import shutil

from dicast_lib.parsing import parse_arguments
from dicast_lib.utils import read_parameters, replace_filename
from dicast_lib.prepare import VariantPrep
from dicast_lib.collect_reference import ReferenceAnnotator
from dicast_lib.collect_illumina import AlignmentAnnotatorIllumina
from dicast_lib.model import Dicast
from dicast_lib.model_exome import DicastExome
from dicast_lib.cohort import Cohort


# List of chromosomes to process
chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
          'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
          'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

# List of SV types currently supported by dicast
sv_types = ['DEL', 'DUP', 'INS', 'INV']


def collect_aln_features(bam_filename: str, variant_filename: str, variant_annot_filename: str, chrom: str, sv_type: str, sample: str, exome_regions: pd.DataFrame=None):
    """ Collects alignment features for a given chromosome. 

    Args:
        bam_filename (str): BAM file
        variant_filename (str):  TSV file with variants to annotate
        variant_annot_filename (str): Output TSV file with annotated variants
        chrom (str): Chromosome name
        sv_type (str): SV type
        sample (str): Sample name
        exome_regions (pd.DataFrame, optional): Exome regions. Defaults to None.
    """    

    AAI = AlignmentAnnotatorIllumina(bam_filename, chrom, sv_type, sample)
    AAI.load_from_csv(variant_filename)
    AAI.calculate_coverage_baseline()
    AAI.calculate_insertsize_baseline()
    AAI.calculate_mapping_quality_baseline()
    AAI.annotate_coverage()
    AAI.annotate_read_based_features()
    AAI.to_csv(variant_annot_filename)


def combine_feature_files(sample: str, ref: str, workdir: str) -> pd.DataFrame:
    """ Combines raw, reference and alignment features into a single TSV file.

    Args:
        sample (str): Sample name
        ref (str): Reference genome name
        workdir (str): Working directory

    Returns:
        pd.DataFrame: DataFrame with combined features
    """    

    df_raw = pd.read_csv(f'{workdir}/{sample}_{ref}.SVs.raw.tsv', sep='\t', low_memory=False)
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

    dicast_df = pd.read_csv(f'{arguments.workdir}/{arguments.sample}_{arguments.ref}.SVs.dicast.tsv', sep='\t')

    for caller, vcf_filename in arguments.vcfs:
        vcf_in = vcfpy.Reader.from_path(vcf_filename)
        vcf_in.header.add_info_line(vcfpy.OrderedDict([('ID', 'DQ'), 
                                                       ('Number', '1'), 
                                                       ('Type', 'String'),  
                                                       ('Description', 'Dicast Quality Score')]))
        vcf_filename_out = vcf_filename.replace('.vcf', '.dicast.vcf').replace('.gz', '')
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
    print(RA.df_calls_annot.shape)
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
    print(RA.df_calls_annot.shape)
    RA.to_csv('/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.ref.tsv']))


def collect_aln_features(bam_filename: str, variant_filename: str, variant_annot_filename: str, chrom: str, sv_type: str, sample: str, exome_regions: pd.DataFrame=None):
    """ Collects alignment features for a given chromosome. 

    Args:
        bam_filename (str): BAM file
        variant_filename (str):  TSV file with variants to annotate
        variant_annot_filename (str): Output TSV file with annotated variants
        chrom (str): Chromosome name
        sv_type (str): SV type
        sample (str): Sample name
        exome_regions (pd.DataFrame, optional): Exome regions. Defaults to None.
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
    dicast_dfs = []
    for sv_type in sv_types:
        if sv_type != 'INV':
            model_filename = f'{arguments.models}/dicast_{sv_type}.json'
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
        logging.info(f'SV CALLERS: {", ".join([caller for caller, _ in arguments.vcfs])}')
        logging.info(f'VCFs: {", ".join([vcf for _, vcf in arguments.vcfs])}')
        if arguments.exome:
                logging.info(f'Exome Target Regions: {arguments.exome_regions}')
        print('')

        # Restrict feature extraction to a single chromosome if specified
        if arguments.chrom != 'all':
            chroms = arguments.chrom
        
        # Variant Preparation
        logging.info('# Create Variant DataFrame')
        VP = VariantPrep(arguments.cohort, arguments.sample, arguments.ref, arguments.workdir, 
                         arguments.technology, arguments.vcfs, chroms, arguments.fai, sv_types)
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
        if arguments.exome:
            logging.info('# Load Exome Target Regions')
            target_regions = pd.read_csv(arguments.exome_regions, sep='\t', names=['chr','start','end'], dtype={'chr':str,'start':int,'stop':int}, index_col=False)

        logging.info(f'# Collect Alignment Features')
        parallel_input = []
        for sv_type in sv_types:
            for chrom in chroms:   
                variant_filename = '/'.join([arguments.workdir, arguments.sample + '_' + arguments.ref + '.SVs.raw.tsv'])
                variant_annot_filename = '/'.join([arguments.workdir, arguments.sample + '_' + arguments.ref + '.SVs.aln.ill.' + chrom + '.' + sv_type + '.tsv'])
                if arguments.exome:
                    target_regions_chrom = target_regions.loc[target_regions['chr']==chrom].copy().reset_index(drop=True)
                    parallel_input.append((arguments.bam, variant_filename, variant_annot_filename, 
                                        chrom, sv_type, arguments.sample, target_regions_chrom))
                else:
                    parallel_input.append((arguments.bam, variant_filename, variant_annot_filename, 
                                        chrom, sv_type, arguments.sample))
        Parallel(n_jobs=arguments.threads)(delayed(collect_aln_features)(*args) for args in parallel_input)
        
        # Combination of Raw, Reference and Alignment Features
        logging.info('# Combination Output Files')
        df = combine_feature_files(arguments.sample, arguments.ref, arguments.workdir)
        df.to_csv(f'{arguments.workdir}/{arguments.sample}_{arguments.ref}.SVs.annot.tsv', sep='\t', index=False, na_rep='NA')
        
        # Variant Prediction
        logging.info('# Variant Prediction')
        score_variants(sv_types, arguments, arguments.sample)
        
        # Add Info Tag to VCF
        logging.info('# Add Info Tag to VCF')
        add_info_tag_to_vcf(arguments)

        logging.info('############### End DICAST ###############\n')

    elif arguments.command == 'cohort':
            
            logging.info('MODE: cohort')
            logging.info(f'COHORT: {arguments.cohort}')
            logging.info(f'CHROMOSOMES: {arguments.chrom}')
            logging.info(f'REF: {arguments.ref}')
            logging.info(f'TECH: {arguments.technology}')
            logging.info(f'WORKDIR: {arguments.workdir}')
            logging.info(f'FAI: {arguments.fai}')
            logging.info(f'REPEATS: {arguments.repeats}')
            logging.info(f'CGIS: {arguments.cgis}')
            logging.info(f'CENTROMERES: {arguments.centromeres}')
            logging.info(f'GAPS: {arguments.gaps}')
            logging.info(f'ALTHAPS: {arguments.althaps}')
            logging.info(f'VNTRS: {arguments.vntrs}')
            logging.info(f'STRS: {arguments.strs}')
            logging.info(f'GC: {arguments.gc}')
            logging.info(f'MODELS: {arguments.models}')
            logging.info(f'THREADS: {arguments.threads}')
            logging.info(f'VCFs: {", ".join([vcf for vcf in arguments.vcfs])}')
            logging.info(f'BAMs: {", ".join([bam for bam in arguments.bams])}')
            print('')
    
            # Restrict feature extration to a single chromosome if specified
            if arguments.chrom != 'all':
                chroms = arguments.chrom
            
            # Create VCF and BAM dictionaries
            vcf_dict = {vcf_path.split('/')[-2]: vcf_path for vcf_path in arguments.vcfs}
            bam_dict = {bam_path.split('/')[-2]: bam_path for bam_path in arguments.bams}
            samples = list(vcf_dict.keys())
            num_samples = len(samples)
            assert len(vcf_dict.keys()) == len(bam_dict.keys())

            # for single sample processing output the same VCF files
            if num_samples == 1:
                vcf_path = vcf_dict[samples[0]]
                vcf_out_path = vcf_path.replace('.vcf', '.regenotyped.vcf')
                # copy the VCF file to the output directory
                shutil.copy(vcf_path, vcf_out_path)
            else:
                # Variant Preparation
                logging.info('# Create Variant DataFrame')
                sample_dfs = []
                sample_dfs_unfiltered = []
                for sample, vcf_file in vcf_dict.items():
                    VP = VariantPrep(arguments.cohort, sample, arguments.ref, arguments.workdir, 
                                    arguments.technology, {sample: vcf_file}, chroms, arguments.fai, 
                                    sv_types, mode='cohort')
                    VP.read_variants()
                    sample_dfs_unfiltered.append(VP.get_variant_df())
                    VP.filter_variants()
                    VP.filter_variants_cohort(num_samples)
                    sample_dfs.append(VP.get_variant_df())
                cohort_df = pd.concat(sample_dfs, ignore_index=True)
                cohort_df_unfiltered = pd.concat(sample_dfs_unfiltered, ignore_index=True)
                cohort_df.to_csv(f'{arguments.workdir}/{arguments.cohort}_{arguments.ref}.SVs.raw.tsv', sep='\t', index=False, na_rep='NA')
                cohort_df_unfiltered.to_csv(f'{arguments.workdir}/{arguments.cohort}_{arguments.ref}.SVs.raw.unfiltered.tsv', sep='\t', index=False, na_rep='NA')
                
                # Determine which variants are missing from each sample
                logging.info('# Cohort Preparation')
                CH = Cohort(cohort_df, cohort_df_unfiltered, samples, arguments.ref, arguments.workdir, vcf_dict)
                CH.get_missing_variants()
                CH.save_missing_variants()
                
                # Collect Reference Features
                for sample in samples:
                    logging.info(f'# Collect Reference Features for {sample}')
                    extract_reference_features(arguments, sample)
                
                # Collect Alignment Features
                if arguments.exome:
                    logging.info('# Load Exome Target Regions')
                    target_regions = pd.read_csv(arguments.exome_regions, sep='\t', names=['chr','start','end'], dtype={'chr':str,'start':int,'stop':int}, index_col=False)

                logging.info(f'# Collect Alignment Features')
                parallel_input = []
                for sample in samples:
                    for sv_type in sv_types:
                        for chrom in chroms:   
                            variant_filename = '/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.raw.tsv'])
                            variant_annot_filename = '/'.join([arguments.workdir, sample + '_' + arguments.ref + '.SVs.aln.ill.' + chrom + '.' + sv_type + '.tsv'])
                            if arguments.exome:
                                target_regions_chrom = target_regions.loc[target_regions['chr']==chrom].copy().reset_index(drop=True)
                                parallel_input.append((bam_dict[sample], variant_filename, variant_annot_filename, 
                                                    chrom, sv_type, sample, target_regions_chrom))
                            else:
                                parallel_input.append((bam_dict[sample], variant_filename, variant_annot_filename, 
                                                    chrom, sv_type, sample))
                                
                Parallel(n_jobs=arguments.threads)(delayed(collect_aln_features)(*args) for args in parallel_input)

                # Combination of Raw, Reference and Alignment Features
                logging.info('# Combination Output Files')
                for sample in samples:
                    df = combine_feature_files(sample, arguments.ref, arguments.workdir)
                    df.to_csv(f'{arguments.workdir}/{sample}_{arguments.ref}.SVs.annot.tsv', sep='\t', index=False, na_rep='NA')
                
                # Variant Prediction
                logging.info('# Variant Prediction')
                for sample in samples:
                    score_variants(sv_types, arguments, sample)
                
                # Add variants to VCF files
                logging.info('# Update cohort information')
                CH.load_dicast_predictions()
                CH.update_cohort_information()

                logging.info('# Find overlapping variants')
                CH.find_overlapping_variants()

                logging.info('# Update VCF files')
                CH.update_vcf_files()

                logging.info('############### End DICAST ###############\n')

