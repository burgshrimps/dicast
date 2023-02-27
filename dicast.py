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
from pysam import VariantFile

from lib.parsing import parse_arguments
from lib.utils import read_parameters, replace_filename
from lib.prepare import VariantPrep
from lib.collect_reference import ReferenceAnnotator
from lib.collect_illumina import AlignmentAnnotatorIllumina
from lib.model import Dicast


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
    
    param chrom: Chromosome name 
    return curation set: pandas dataframe with SVs for manual curation """

    dicast = Dicast('curate', arguments.svtype, params, clf=arguments.clfname, clfparams=clfparams, chr_excl=[chrom], chr_incl=[chrom])
    dicast.load_data()
    dicast.train()
    dicast.predict()
    dicast.compute_qual()
    print('Finished chromosome {0}'.format(chrom))
    return dicast.get_curation_set()


def save_manual_curation_set(df, params, chunk_size=200):
    """ Splits manual curation dataframe into chunks and saves them to disk.

    param df: Manual curation dataframe
    param params: Parameters
    param chunk_size: Number of SVs per chunk """

    for cohort in params:
        ref = params[cohort]['ref']

        for sample in params[cohort]['samples']:
            df_sample = df.loc[(df['sample'] == sample)].copy().reset_index(drop=True)
            df_sample_fp = df_sample[df_sample['err_type'] == 'FP'].copy().reset_index(drop=True)
            df_sample_fn = df_sample[df_sample['err_type'] == 'FN'].copy().reset_index(drop=True)

            workdir_root_sample = replace_filename(params[cohort]['variant_curation'], sample, ref)

            # Save False Positives
            chunk_fp = 0
            for i in range(0, len(df_sample_fp), chunk_size):
                workdir_sample_fp_chunk = workdir_root_sample + '/FP/' + arguments.svtype + '/chunk' + str(chunk_fp)
                filename_sample_fp_chunk = '_'.join([datetime.today().strftime('%Y%m%d'), sample, 'FP', arguments.svtype, 'chunk' + str(chunk_fp)]) + '.tsv'
                if not os.path.exists(workdir_sample_fp_chunk):
                    os.makedirs(workdir_sample_fp_chunk)

                df_sample_fp_chunk = df_sample_fp.iloc[i:i+chunk_size].copy().reset_index(drop=True)
                df_sample_fp_chunk.to_csv('/'.join([workdir_sample_fp_chunk, filename_sample_fp_chunk]), sep='\t', index=False, na_rep='NA')

                chunk_fp += 1

            # Save False Negatives
            chunk_fn = 0
            for i in range(0, len(df_sample_fn), chunk_size):
                workdir_sample_fn_chunk = workdir_root_sample + '/FN/' + arguments.svtype + '/chunk' + str(chunk_fn)
                filename_sample_fn_chunk = '_'.join([datetime.today().strftime('%Y%m%d'), sample, 'FN', arguments.svtype, 'chunk' + str(chunk_fn)]) + '.tsv'
                if not os.path.exists(workdir_sample_fn_chunk):
                    os.makedirs(workdir_sample_fn_chunk)

                df_sample_fn_chunk = df_sample_fn.iloc[i:i+chunk_size].copy().reset_index(drop=True)
                df_sample_fn_chunk.to_csv('/'.join([workdir_sample_fn_chunk, filename_sample_fn_chunk]), sep='\t', index=False, na_rep='NA')

                chunk_fn += 1


def save_predictions(workdir, pred):
    """ Splits dicast prediction dataframe by sample and saves it in workdir. 
    
    param workdir: Output directory
    param pred: Dicast prediction dataframe """

    for sample in pred['sample'].unique():
        pred_sample = pred.loc[(pred['sample'] == sample)].copy().reset_index(drop=True)
        pred_sample = pred_sample[['id', 'sample', 'tech', 'method', 'type', 'chrom', 'chrom2', 'start', 'end', 'size', 
                                   'filter', 'qual', 'pred_dicast', 'qual_dicast']]
        pred_sample.to_csv('/'.join([workdir, sample + '.dicast.tsv']), sep='\t', index=False, na_rep='NA')


def add_predictions_vcf(params, predictions, workdir):
    """ Adds Dicast predictions as a new INFO field to the VCF file 
    
    param params: dictionary, Prediction parameter file
    param predictions: pandas dataframe, dicast predictions
    param workdir: Output directory """

    for cohort in params:
        for sample in params[cohort]['samples']:
            for tech in params[cohort]['vcf']:
                for method in params[cohort]['vcf'][tech]:

                    fname_vcf_in = replace_filename(params[cohort]['vcf'][tech][method], sample, params[cohort]['ref'])
                    if not os.path.exists(fname_vcf_in):
                        fname_vcf_in = fname_vcf_in.replace(sample, sample.replace('-', '_'))
                    vcf_in = VariantFile(fname_vcf_in, 'r')
                    vcf_in.header.info.add("DQ", number="1", type="String", description="Dicast Quality Score")

                    fname_vcf_out = fname_vcf_in.replace('.vcf', '.dicast.vcf').replace('.gz', '')
                    vcf_out = VariantFile(fname_vcf_out, 'w', header=vcf_in.header)
                    
                    variants_subset = predictions.loc[(predictions['sample'] == sample) & (predictions['tech'] == tech) & (predictions['method'] == method)].copy().reset_index(drop=True)
                    variants_subset = variants_subset[['id', 'qual_dicast']].set_index('id').T.to_dict('list')
                    for rec in vcf_in.fetch():
                        if rec.id in variants_subset.keys():
                            qual_dicast = variants_subset[rec.id][0]
                        else:
                            qual_dicast = -1

                        rec.info['DQ'] = str(qual_dicast)
                        vcf_out.write(rec)
                    vcf_out.close()
                    vcf_in.close()
                    logging.info(f'Added DQ tag to {method} VCF file for sample {sample}')



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
        
        VP = VariantPrep(arguments.sample, arguments.ref, params, arguments.workdir, CHROMS)
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
        logging.info(f'CURATION: {arguments.cur}')
        print('')

        params = read_parameters(arguments.params)
        clfparams = read_parameters(arguments.clfparams)
        model_out = clfparams[arguments.clfname]['directory'] + '/' + arguments.clfname + '_' + arguments.svtype.upper() + '.pkl'

        dicast = Dicast('train', arguments.svtype, params, pkl=model_out, clf=arguments.clfname, clfparams=clfparams, chr_excl=arguments.chr_excl, incl_cur=arguments.cur)
        logging.info('# Load Training Data')
        dicast.load_data()
        logging.info('# Train Model')
        dicast.train()
        logging.info('# Save Model')
        dicast.save_model()

    elif arguments.command == 'test':

        logging.info('MODE: test')
        logging.info(f'TYPE: {arguments.svtype}')
        logging.info(f'PARAMS: {arguments.params}')
        logging.info(f'INCLUDED CHROMS: {arguments.chr_incl}')
        logging.info(f'CURATION: {arguments.cur}')
        print('')

        params = read_parameters(arguments.params)
        clfparams = read_parameters(arguments.clfparams)
        model_in = clfparams[arguments.clfname]['directory'] + '/' + arguments.clfname + '_' + arguments.svtype.upper() + '.pkl'

        dicast = Dicast('test', arguments.svtype, params, pkl=model_in, clf=arguments.clfname, clfparams=clfparams, chr_incl=arguments.chr_incl, incl_cur=arguments.cur)
        logging.info('# Load Test Data')
        dicast.load_data()
        logging.info('# Load Model')
        dicast.load_model()
        logging.info('# Predict')
        dicast.predict()
        logging.info('# Compute Quality Scores')
        dicast.compute_qual()
        logging.info('# Save Test Results')
        dicast.save_test()


    elif arguments.command == 'predict':

        logging.info('MODE: test')
        logging.info(f'PARAMS: {arguments.params}')
        print('')

        params = read_parameters(arguments.params)
        clfparams = read_parameters(arguments.clfparams)
        model_dir = clfparams[arguments.clfname]['directory']

        predictions = []
        for file in glob(model_dir + '/*.pkl'):
            svtype = file.split('_')[-1].split('.')[0]
            logging.info(f'Predicting: {svtype}')
            dicast = Dicast('predict', svtype, params, pkl=file, clf=arguments.clfname, clfparams=clfparams)
            dicast.load_data()
            dicast.load_model()
            dicast.predict()
            predictions.append(dicast.variants)
        logging.info('# Save Predictions')
        save_predictions(arguments.workdir, pd.concat(predictions, ignore_index=True))

        if arguments.vcf:
            add_predictions_vcf(params, pd.concat(predictions, ignore_index=True), arguments.workdir)

    elif arguments.command == 'curate':
        logging.info('MODE: curate')
        logging.info(f'TYPE: {arguments.svtype}')
        logging.info(f'CLASSIFIER: {arguments.clfname}')
        logging.info(f'CLASSIFIER PARAMS: {arguments.clfparams}')
        logging.info(f'PARAMS: {arguments.params}')
        print('')

        params = read_parameters(arguments.params)
        clfparams = read_parameters(arguments.clfparams)

        curation_set = Parallel(n_jobs=len(CHROMS))(delayed(create_manual_curation_set)(chrom) for chrom in CHROMS)
        curation_set = pd.concat(curation_set, ignore_index=True)
        save_manual_curation_set(curation_set, params)
    else:
        raise ValueError('Invalid command')

    print('')
    logging.info('############### Finished DICAST ###############\n')


    