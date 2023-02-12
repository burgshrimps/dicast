import argparse
import logging
import sys
from joblib import Parallel, delayed
import pandas as pd
from glob import glob
from tqdm import tqdm

from lib.parsing import parse_arguments
from lib.utils import read_parameters
from lib.prepare import VariantPrep
from lib.collect_reference import ReferenceAnnotator
from lib.collect_illumina import AlignmentAnnotatorIllumina
from lib.model import Dicast
import os
import re

# List of chromosomes to process
CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
          'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
          'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']


def collect_aln_features(chrom):
    """ Collects alignment features for a given chromosome. 
    
    param chrom: Chromosome name """

    AAI = AlignmentAnnotatorIllumina(arguments.sample, arguments.ref, arguments.workdir, params, chrom=chrom)
    AAI.calculate_coverage_baseline()
    AAI.calculate_insertsize_baseline()
    AAI.calculate_mapping_quality_baseline()
    AAI.annotate_coverage()
    AAI.annotate_read_based_features()
    AAI.to_csv()


def combine_feature_files(sample, ref, workdir):
    """ Combines variant, reference and alignment features into one file.
    
    param sample: Sample name
    param ref: Reference name
    param workdir: Working directory 
    
    return: pandas dataframe with combined features """

    df_raw = pd.read_csv(f'{workdir}/ensemble/{sample}_{ref}.SVs.raw.tsv', sep='\t', low_memory=False)
    df_ref = pd.read_csv(f'{workdir}/ensemble/{sample}_{ref}.SVs.ref.tsv', sep='\t', low_memory=False)
    filenames_aln_ill = glob(f'{workdir}/ensemble//{sample}_{ref}.SVs.aln.ill.*.tsv')
    df_aln_ill = pd.concat([pd.read_csv(f, sep='\t') for f in filenames_aln_ill], ignore_index=True)
    df = df_raw.merge(df_ref.drop(['sample', 'type', 'chrom', 'chrom2', 'start', 'end'], axis=1), on='id', how='inner')
    df = df.merge(df_aln_ill.drop(['sample', 'type', 'chrom', 'chrom2', 'start', 'end'], axis=1), on='id', how='inner')

    return df


def create_manual_curation_set(chrom):
    """ Create a set of SVs for manual curation. Train on all chromosomes except the current one and then test on the current one.
    
    param chrom: Chromosome name """

    dicast = Dicast('curate', arguments.svtype, params, clf=arguments.clfname, clfparams=clfparams, chr_excl=[chrom], chr_incl=[chrom])
    dicast.load_data()
    dicast.train()
    dicast.predict()
    dicast.compute_qual()
    dicast.save_curation()
    print('Finished chromosome {0}'.format(chrom))


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

    if arguments.command == 'prepare':

        logging.info('MODE: prepare')
        logging.info(f'SAMPLE: {arguments.sample}')
        logging.info(f'REF: {arguments.ref}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'WORKDIR: {arguments.workdir}')
        print('')
        
        params = read_parameters(arguments.params)
        
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

        logging.info('# Collect Illumina Alignment Features')
        Parallel(n_jobs=len(CHROMS))(delayed(collect_aln_features)(chrom) for chrom in CHROMS)
        
        logging.info('# Combination Output Files')
        df = combine_feature_files(arguments.sample, arguments.ref, arguments.workdir)
        df.to_csv(f'{arguments.workdir}/ensemble/{arguments.sample}_{arguments.ref}.SVs.annot.tsv', sep='\t', index=False, na_rep='NA')

    elif arguments.command == 'train':

        logging.info('MODE: train')
        logging.info(f'TYPE: {arguments.svtype}')
        logging.info(f'CLASSIFIER: {arguments.clfname}')
        logging.info(f'CLASSIFIER PARAMS: {arguments.clfparams}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'EXCLUDED CHROMS: {arguments.chr_excl}')
        print('')

        params = read_parameters(arguments.params)
        clfparams = read_parameters(arguments.clfparams)
        output = clfparams[arguments.clfname]['directory'] + '/' + arguments.clfname + '_' + arguments.svtype.upper() + '.pkl'

        dicast = Dicast('train', arguments.svtype, params, pkl=output, clf=arguments.clfname, clfparams=clfparams, chr_excl=arguments.chr_excl)
        logging.info('# Load Training Data')
        dicast.load_data()
        logging.info('# Train Model')
        dicast.train()
        logging.info('# Save Model')
        dicast.save_model()

    elif arguments.command == 'test':

        logging.info('MODE: test')
        logging.info(f'TYPE: {arguments.svtype}')
        logging.info(f'MODEL INPUT: {arguments.clf}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'INCLUDED CHROMS: {arguments.chr_incl}')
        print('')

        params = read_parameters(arguments.params)
        dicast = Dicast('test', arguments.svtype, params, pkl=arguments.clf, chr_incl=arguments.chr_incl)
        logging.info('# Load Test Data')
        dicast.load_data()
        logging.info('# Load Model')
        dicast.load_model()
        logging.info('# Predict')
        dicast.predict()
        logging.info('# Compute Quality Scores')
        dicast.compute_qual()
        logging.info('# Save Test Results')
        dicast.save_test(arguments.workdir)


    elif arguments.command == 'predict':

        logging.info('MODE: test')
        logging.info(f'MODEL INPUT: {arguments.clfdir}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'PREDICTIONS OUTPUT: {arguments.output}')
        print('')
        params = read_parameters(arguments.params)
        print(params)
        print(arguments)
        l_pred = []
        for file in os.listdir(arguments.clfdir):
            if file.endswith(".pkl"):
                saved_model = os.path.join(arguments.clfdir , file)
                logging.info(f'READING FILE: {saved_model}')
                if len(re.split(r'[_\.]', file)) != 3:
                    logging.info('ERROR: model name must be in the format <model>_<svtype>.pkl')
                    exit()
                model, svtype, _ = re.split(r'[_\.]', file)
                dicast = Dicast('predict', svtype, params, pkl=saved_model)
                logging.info('# Load Data')
                dicast.load_data()
                logging.info('# Load Model')
                dicast.load_model()
                logging.info('# Predict')
                dicast.predict()
                l_pred.append(dicast.variants)
        logging.info('# Save Predictions')
        dicast.save_predictions_tsv(arguments.output, pd.concat(l_pred))
        logging.info('# add prediction value to vcf file')
        dicast.add_predictions_to_vcf(pd.concat(l_pred))

    elif arguments.command == 'curate':
        logging.info('MODE: curate')
        logging.info(f'TYPE: {arguments.svtype}')
        logging.info(f'CLASSIFIER: {arguments.clfname}')
        logging.info(f'CLASSIFIER PARAMS: {arguments.clfparams}')
        logging.info(f'PARAMS: {arguments.params}')
        print('')

        params = read_parameters(arguments.params)
        clfparams = read_parameters(arguments.clfparams)

        Parallel(n_jobs=len(CHROMS))(delayed(create_manual_curation_set)(chrom) for chrom in CHROMS)
    else:
        raise ValueError('Invalid command')

    print('')
    logging.info('############### Finished DICAST ###############\n')


    