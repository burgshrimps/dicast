import pandas as pd 
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
import pickle
import logging
import bioframe as bf
import numpy as np
from datetime import datetime
import os
import pysam
import pathlib
from lib.utils import replace_filename
import sys

class Dicast:
    """ Class to train a model. """

    def __init__(self, mode, svtype, params, pkl=None, clf=None, clfparams=None, chr_excl=[], chr_incl=['all']):
        """ Initialize class. 
        
        param mode: string, 'train', 'test' or 'predict'
        param svtype: string, SV type
        param pkl: string, model input/output file
        param params: dictionary, parameters including location of data 
        param clf: string, classifier name
        param chr_excl: list, chromosomes to exclude
        param chr_incl: list, chromosomes to include """

        self.mode = mode
        self.type = svtype.upper()
        self.params = params
        self.chr_excl = chr_excl
        self.chr_incl = chr_incl

        if pkl != None:
            self.pkl = pkl
        if clf != None:
            self.clf = clf
        if clfparams != None:
            self.clfparams = clfparams[clf]


    def load_sample(self, sample, ref, variant_features, variant_labels=None):
        """ Load sample data. 
        
        param sample: string, sample name 
        param ref: string, reference genome name 
        param variant_features: string, variant feature file
        param variant_labels: string, variant label file
        
        return pandas dataframe with sample data """

        # Load variant features
        filename_features = replace_filename(variant_features, sample, ref)
        df_features = pd.read_csv(filename_features, sep='\t', low_memory=False)

        # Load variant labels
        if variant_labels != None:
            filename_labels = replace_filename(variant_labels, sample, ref)
            df_labels = pd.read_csv(filename_labels, low_memory=False)

            # Filter based on confirmation status
            df_labels = df_labels[(df_labels['StatusSimple'] == 'Confirmed') | (df_labels['StatusSimple'] == 'ConfirmedPublic')].copy().reset_index(drop=True)

            # Extract IDs of confirmed variants
            confirmed_ids = defaultdict(list)
            methods = df_features['method'].unique()
            for i in range(len(df_labels)):
                sub_graph = df_labels.loc[i, 'sub_graph'][1:-1].split(', ')
                for entry in sub_graph:
                    for method in methods:
                        if entry[1:-1].startswith(method):
                            confirmed_ids[method].append(entry[1:-1].split('_')[-1])

            # Add confirmation status to feature dataframe
            df_features['confirmed'] = 0
            for i in range(len(df_features)):
                if df_features.loc[i, 'id'] in confirmed_ids[df_features.loc[i, 'method']]:
                    df_features.loc[i, 'confirmed'] = 1 

        return df_features


    def load_cohort(self, cohort):
        """ Load cohort data. 
        
        param cohort: string, cohort name """
        
        samples = self.params[cohort]['samples']
        ref = self.params[cohort]['ref']
        variant_features = self.params[cohort]['variant_features']
        sample_dfs = []

        if self.mode == 'train' or self.mode == 'test' or self.mode == 'curate':
            variant_labels = self.params[cohort]['variant_labels']
            for sample in samples:
                sample_df = self.load_sample(sample, ref, variant_features, variant_labels=variant_labels)
                sample_dfs.append(sample_df)

        elif self.mode == 'predict':
            for sample in samples:
                sample_df = self.load_sample(sample, ref, variant_features)
                sample_dfs.append(sample_df)

        else:
            raise ValueError(f'Invalid mode: {self.mode}')

        return pd.concat(sample_dfs, ignore_index=True)


    def load_data(self):
        """ Load data. """

        cohort_dfs = []
        for cohort in self.params:
            cohort_df = self.load_cohort(cohort)
            cohort_dfs.append(cohort_df)

        self.variants = pd.concat(cohort_dfs, ignore_index=True)

        # Filter based on SV type
        self.variants = self.variants[self.variants['type'] == self.type].copy().reset_index(drop=True)

        # Order based on confirmation status
        if self.mode == 'train' or self.mode == 'test' or self.mode == 'curate':
            self.variants = self.variants.sort_values(by=['confirmed'], ascending=True).reset_index(drop=True)


    def check_caller_support(self, row, min_num_callers):
        """ Adds caller support to dataframe. 
        
        param row: pandas dataframe row
        param min_num_callers: int, minimum number of callers required to call a variant
        
        return int, caller support """

        caller_count = 0
        for qual in row:
            if qual != 0:
                caller_count += 1
        if caller_count >= min_num_callers:
            return row.sum()
        else:
            return 0


    def extract_overlap_ids(self, df1, df2):
        """ Extracts SV IDs of overlapping variants
        
        param df1: pandas dataframe, first dataframe
        param df2: pandas dataframe, second dataframe 
        
        return pandas dataframe with overlapping SV IDs """

        closest_intervals = bf.closest(df1, df2, suffixes=('_1','_2'))
        closest_intervals = closest_intervals.dropna(subset=['id_1', 'id_2']).reset_index(drop=True)
        closest_intervals['diff_start'] = abs(closest_intervals['start_1'] - closest_intervals['start_2'])
        closest_intervals['diff_end'] = abs(closest_intervals['end_1'] - closest_intervals['end_2'])
        closest_intervals['diff_size'] = closest_intervals.apply(lambda x: min([x['size_1'], x['size_2']]) / max([x['size_1'], x['size_2']]), axis=1)
        overlapping_svs = closest_intervals[(closest_intervals['diff_start'] < 50) & (closest_intervals['diff_end'] < 50) & (closest_intervals['diff_size'] > 0.7)].copy()
        
        return overlapping_svs[['id_1', 'id_2']].reset_index(drop=True)


    def compute_qual(self):
        """ Checks which SV is called by multiple methods and adds the quality score of the other method to the result dataframe. """

        variants_samples = []
        for sample in self.variants['sample'].unique():
            variants_sample = self.variants[self.variants['sample'] == sample].copy().reset_index(drop=True)
            methods = variants_sample['method'].unique()

            # Add quality score for each method to dataframe
            for method in methods:
                variants_sample['qual_' + method] = 0
                variants_sample.loc[variants_sample['method'] == method, 'qual_' + method] = variants_sample.loc[variants_sample['method'] == method, 'qual']

            # Check for overlapping SVs
            for i in range(len(methods)):
                for j in range(i+1, len(methods)):
                    df1 = variants_sample[variants_sample['method'] == methods[i]].copy().reset_index(drop=True)
                    df2 = variants_sample[variants_sample['method'] == methods[j]].copy().reset_index(drop=True)    
                    df1 = df1[['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)
                    df2 = df2[['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)
                    overlap_ids = self.extract_overlap_ids(df1, df2)

                    # Add quality scores of other method to result dataframe
                    for k in range(len(overlap_ids)):
                        variants_sample.loc[variants_sample['id'] == overlap_ids['id_1'][k], 'qual_' + methods[j]] = variants_sample.loc[variants_sample['id'] == overlap_ids['id_2'][k], 'qual_' + methods[j]].values
                        variants_sample.loc[variants_sample['id'] == overlap_ids['id_2'][k], 'qual_' + methods[i]] = variants_sample.loc[variants_sample['id'] == overlap_ids['id_1'][k], 'qual_' + methods[i]].values

            # Normalize Qualities and Round
            for method in methods:
                variants_sample['qual_' + method] = variants_sample['qual_' + method] / variants_sample['qual_' + method].max()
                variants_sample['qual_' + method] = variants_sample['qual_' + method].apply(lambda x: np.round(x, 2))

            variants_samples.append(variants_sample)

        self.variants = pd.concat(variants_samples, ignore_index=True)
        self.variants.drop('qual', axis=1, inplace=True)

        # For each variant add caller support
        for i in range(len(self.variants['method'].unique())):
            self.variants['qual_' + str(i+1) + '_caller_support'] = self.variants.apply(lambda x: self.check_caller_support(x[['qual_' + method for method in methods]], i+1), axis=1)


    def train(self):
        """ Train model. """

        # Filter based on excluded chromosomes
        variants_training = self.variants[~self.variants['chrom'].isin(self.chr_excl)].copy().reset_index(drop=True)

        features = list(variants_training.columns[12:-1]) + ['size']

        # Check how many variants are excluded due to missing feature values
        num_svs = len(variants_training)
        variants_training = variants_training.dropna(subset=features).copy().reset_index(drop=True)
        logging.info(f'# Dropped {num_svs - len(variants_training)}/{num_svs} SVs due to missing feature values')

        X = variants_training[features]
        y = variants_training['confirmed']

        if self.clfparams['classifier'] == 'RandomForestClassifier':
            self.model = RandomForestClassifier(**self.clfparams['parameters'])
        else:
            raise ValueError(f'Invalid classifier: {self.clfparams["classifier"]}')

        self.model.fit(X, y)


    def predict(self):
        """ Predict probability of being TP for set of variants. """

        # Filter based on included chromosomes
        if self.chr_incl[0] != 'all':
            self.variants = self.variants[self.variants['chrom'].isin(self.chr_incl)].copy().reset_index(drop=True)

        features = list(self.variants.columns[12:]) + ['size']

        # Remove confirmed column if present, this is the case during model curation
        if 'confirmed' in features:
            features.remove('confirmed')

        # Check how many variants are excluded due to missing feature values
        num_svs = len(self.variants)
        self.variants = self.variants.dropna(subset=features).copy().reset_index(drop=True)
        logging.info(f'# Dropped {num_svs - len(self.variants)}/{num_svs} SVs due to missing feature values')

        X = self.variants[features]

        self.variants['pred_dicast'] = self.model.predict(X)
        self.variants['qual_dicast'] = self.model.predict_proba(X)[:, 1]


    def save_curation(self):
        """ Save curation SVs. """

        qual_cols = [col for col in self.variants.columns if col.startswith('qual_')]
        variants_curation = self.variants[['sample', 'tech', 'method', 'id', 'type', 'chrom', 'chrom2', 
                                           'start', 'end', 'size', 'filter', 'confirmed', 
                                           'pred_dicast'] + qual_cols].copy().reset_index(drop=True)

        # Determine which variants are sent for curation
        variants_curation_fp = variants_curation[(variants_curation['confirmed'] == 0) & ((variants_curation['qual_dicast'] > 0.4) | (variants_curation['qual_1_caller_support'] > 0.4) | (variants_curation['qual_2_caller_support'] > 0.1))].copy().reset_index(drop=True)
        variants_curation_fn = variants_curation[(variants_curation['confirmed'] == 1) & (variants_curation['qual_dicast'] < 0.4)].copy().reset_index(drop=True)
        variants_curation = pd.concat([variants_curation_fp, variants_curation_fn], ignore_index=True)
        logging.info(f'# Sending {len(variants_curation)} SVs from {self.chr_incl[0]} for manual curation')
        
        for cohort in self.params:
            ref = self.params[cohort]['ref']
            for sample in self.params[cohort]['samples']:
                variants_curation_sample = variants_curation[variants_curation['sample'] == sample].copy().reset_index(drop=True)

                workdir_root_sample = replace_filename(self.params[cohort]['workdir'], sample, ref)
                workdir_sample = workdir_root_sample + '/' + self.type + '/' + self.chr_incl[0]
                if not os.path.exists(workdir_sample):
                    os.makedirs(workdir_sample)

                filename_sample = '_'.join([datetime.today().strftime('%Y%m%d'), sample, self.type, self.chr_incl[0]]) + '.tsv'
                variants_curation_sample.to_csv('/'.join([workdir_sample, filename_sample]), sep='\t', index=False, na_rep='NA')

        
    def save_predictions_tsv(self, output_file, concatinated_dfs):
        """ Save predictions. """
        abs_path = os.path.abspath(output_file)
        logging.info(f'# Saving predictions to {abs_path}')
        concatinated_dfs[['id', 'sample', 'tech', 'method', 'type', 'chrom', 'chrom2', 'start', 'end', 'size', 
                       'filter', 'qual', 'pred_dicast', 'qual_dicast']].to_csv( output_file , sep='\t', index=False, na_rep='NA')


    def add_predictions_to_vcf(self, df_variants):
        """ Add predictions to VCF. """
        l_vcfsnames = []
        for cohort in self.params:
            for sample in self.params[cohort]['samples']:
                for techs in self.params[cohort]['vcf']:
                    for vcf in self.params[cohort]['vcf'][techs]:
                        filename = replace_filename(self.params[cohort]['vcf'][techs][vcf], sample, self.params[cohort]['ref'])
                        if not os.path.exists(filename): 
                            filename = replace_filename(self.params[cohort]['vcf'][techs][vcf], sample.replace('-','_'), self.params[cohort]['ref'])
                        if not os.path.exists(filename): # Try to find file in workdir, with a misspellings in sample name
                            logging.error(f'VCF file {filename} does not exist')
                            sys.exit(1)
                        l_vcfsnames += [(vcf,filename)]
        for _, path in l_vcfsnames:
            input = pysam.VariantFile(path, 'r')
            input.header.info.add("DICAST", number="1", type="String", description="Dicast prediction score")
            new_vcf = path.replace('.vcf.gz','.dicast.vcf.gz')
            output = pysam.VariantFile(new_vcf, 'w', header=input.header)
            for record in input.fetch():
                if record.id in df_variants['id'].values:
                    record.info['DICAST'] = str(df_variants[df_variants['id'] == record.id]['qual_dicast'].values[0])
                else:
                    record.info['DICAST'] = '-1'
                output.write(record)

            input.close()
            logging.info(f'VCF file {new_vcf} has been updated with DICAST predictions')
            output.close()

    def save_test(self, out_root):
        """ Save test result. """

        self.variants.to_csv(out_root + '/test.tsv', sep='\t', index=False, na_rep='NA')


    def load_model(self):
        """ Load model. """

        with open(self.pkl, 'rb') as f:
            self.model = pickle.load(f)
        
        
    def save_model(self):
        """ Save model. """

        with open(self.pkl, 'wb') as f:
            pickle.dump(self.model, f)