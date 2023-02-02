import argparse
import logging
import sys
from joblib import Parallel, delayed
import pandas as pd

from lib.parsing import parse_arguments
from lib.utils import read_parameters
from lib.prepare import VariantPrep
from lib.feature_collection_reference import ReferenceAnnotator
from lib.feature_collection_illumina import AlignmentAnnotatorIllumina

def collect_aln_features(chrom):
    """ Collects alignment features for a given chromosome. """

    AAI = AlignmentAnnotatorIllumina(arguments.sample, arguments.ref, arguments.workdir, params, chrom=chrom)
    AAI.calculate_coverage_baseline()
    AAI.calculate_insertsize_baseline()
    AAI.calculate_mapping_quality_baseline()
    AAI.annotate_coverage()
    AAI.annotate_read_based_features()
    AAI.to_csv()


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
    logging.info('############### Start DICAST-PREP ###############\n')
    logging.info('CMD: python3 {0}'.format(' '.join(sys.argv)))

    if arguments.command == 'prepare':

        logging.info('MODE: prepare')
        logging.info(f'SAMPLE: {arguments.sample}')
        logging.info(f'REF: {arguments.ref}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'WORKDIR: {arguments.workdir}')
        print('')
        
        params = read_parameters(arguments.params)
        """
        VP = VariantPrep(arguments.sample, arguments.ref, params, arguments.workdir)
        logging.info('# Read Variants')
        VP.read_variants() 
        logging.info('# Filter Variants')
        VP.filter_variants()
        logging.info('# Save Variants')
        VP.save_variants()

        logging.info('# Collect Reference Features')
        RA = ReferenceAnnotator(arguments.sample, arguments.ref, arguments.workdir, params)
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
        logging.info('# Save Reference Features')
        RA.to_csv()
        """
        logging.info('# Collect Illumina Alignment Features')
        #chroms = params['chroms']
        chroms = ['chr19', 'chr20']
        Parallel(n_jobs=len(chroms))(delayed(collect_aln_features)(chrom) for chrom in chroms)

        logging.info('# Combination Output Files')
        df_raw = pd.read_csv(f'{arguments.workdir}/ensemble/{arguments.sample}_{arguments.ref}.SVs.raw.tsv', sep='\t')
        df_ref = pd.read_csv(f'{arguments.workdir}/ensemble/{arguments.sample}_{arguments.ref}.SVs.ref.tsv', sep='\t', low_memory=False)
        filenames_aln_ill = glob(f'{arguments.workdir}/ensemble//{arguments.sample}_{arguments.ref}.SVs.aln.ill.*.tsv')
        df_aln_ill = pd.concat([pd.read_csv(f, sep='\t') for f in filenames_aln_ill], ignore_index=True)
        df = df_raw.merge(df_ref.drop(['sample', 'type', 'chrom', 'chrom2', 'start', 'end'], axis=1), on='id', how='inner')
        df = df.merge(df_aln_ill.drop(['sample', 'type', 'chrom', 'chrom2', 'start', 'end'], axis=1), on='id', how='inner')
        df.to_csv(f'{arguments.workdir}/ensemble/{arguments.sample}_{arguments.ref}.SVs.annot.tsv', sep='\t', index=False, na_rep='NA')

        logging.info('############### Finished DICAST-PREP ###############\n')

    elif arguments.command == 'train':

        logging.info('MODE: train')
        logging.info(f'REF: {arguments.ref}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'WORKDIR: {arguments.workdir}')
        print('')


    