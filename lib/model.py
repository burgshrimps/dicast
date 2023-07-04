import pandas as pd 
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import resample
import pickle
import logging
import bioframe as bf
import numpy as np
from datetime import datetime
import os
import pysam
import pathlib
import sys
import json
import xgboost as xgb

from lib.utils import replace_filename


class Dicast:
    """ Class to train a model. """

    def __init__(self, mode, svtype, params, pkl=None, clf=None, clfparams=None, chr_excl=[], chr_incl=['all'], incl_cur=False, balance=False):
        """ Initialize class. 
        
        param mode: string, 'train', 'test' or 'predict'
        param svtype: string, SV type
        param pkl: string, model input/output file
        param params: dictionary, parameters including location of data 
        param clf: string, classifier name
        param chr_excl: list, chromosomes to exclude
        param chr_incl: list, chromosomes to include 
        param incl_cur: boolean, correct curated variants in training data """

        self.mode = mode
        self.type = svtype.upper()
        self.params = params
        self.chr_excl = chr_excl
        self.chr_incl = chr_incl
        self.incl_cur = incl_cur
        self.balance = balance

        if pkl != None:
            self.pkl = pkl
        if clf != None:
            self.clf = clf
        if clfparams != None:
            self.clfparams = clfparams[clf]

        self.cov_thr = 6 # Threshold for log2 change in coverage to be considered for feature extraction, otherwise jump


    def load_sample(self, sample, ref, variant_features, variant_labels=None, variant_curated=None):
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
            df_labels = pd.read_csv(filename_labels, low_memory=False, index_col=0)

            # Filter based on confirmation status
            df_labels = df_labels[(df_labels['FinalStatus'] == 'Confirmed') | (df_labels['FinalStatus'] == 'ConfirmedPublic')].copy().reset_index(drop=True)

            # Extract IDs of confirmed variants
            confirmed_ids = defaultdict(list)
            methods = df_features['method'].unique()
            for i in range(len(df_labels)):
                sub_graph = df_labels.loc[i, 'sub_graph'][1:-1].split(', ')
                for entry in sub_graph:
                    for method in methods:
                        if entry[1:-1].startswith(method):
                            confirmed_ids[method].append(entry[1:-1])

            # Add confirmation status to feature dataframe
            df_features['confirmed'] = 0
            for i in range(len(df_features)):
                if df_features.loc[i, 'id'] in confirmed_ids[df_features.loc[i, 'method']]:
                    df_features.loc[i, 'confirmed'] = 1 

            # Load variant curation
            if variant_curated != None:
                filename_curated = replace_filename(variant_curated, sample, ref)
                df_curated = pd.read_csv(filename_curated, low_memory=False, sep='\t') 

                confirmed_ids_cur = defaultdict(list)
                unconfirmed_ids_cur = defaultdict(list)
                for method in methods:
                    df_curated_method = df_curated[df_curated['method'] == method].copy().reset_index(drop=True)
                    confirmed_ids_cur[method] = df_curated_method[df_curated_method['Confirmed (Consensus)'] == 1]['id'].tolist()
                    unconfirmed_ids_cur[method] = df_curated_method[df_curated_method['Confirmed (Consensus)'] == 0]['id'].tolist()

                for i in range(len(df_features)):
                    if df_features.loc[i, 'id'] in confirmed_ids_cur[df_features.loc[i, 'method']]:
                        df_features.loc[i, 'confirmed'] = 1
                    elif df_features.loc[i, 'id'] in unconfirmed_ids_cur[df_features.loc[i, 'method']]:
                        df_features.loc[i, 'confirmed'] = 0

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
                
                if self.incl_cur:
                    variant_curated = self.params[cohort]['variant_curation']
                    sample_df = self.load_sample(sample, ref, variant_features, variant_labels=variant_labels, variant_curated=variant_curated)
                else:
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
        
        #### QUICK FIX RM LATER
        self.variants['size'].fillna('', inplace=True)
        self.variants.loc[self.variants['size'].str.startswith('('), 'size'] = self.variants.loc[self.variants['size'].str.startswith('('), 'size'].str.extract('\(([^,]*),')
        self.variants['size'].replace('', np.nan, inplace=True)
        self.variants['size'] = self.variants['size'].astype(float)
        ############

        # Filter based on SV type
        self.variants = self.variants[self.variants['type'] == self.type].copy().reset_index(drop=True)

        # Order based on confirmation status
        if self.mode == 'train' or self.mode == 'test' or self.mode == 'curate':
            self.variants = self.variants.sort_values(by=['confirmed'], ascending=True).reset_index(drop=True)

        # Resample to balance classes, because in all cases we have more negative than positive examples
        if self.mode == 'train':
            if self.balance:
                variants_confirmed = self.variants[self.variants['confirmed'] == 1].copy().reset_index(drop=True)
                variants_unconfirmed = self.variants[self.variants['confirmed'] == 0].copy().reset_index(drop=True)

                variants_confirmed_upsampled = resample(variants_confirmed, replace=True, n_samples=len(variants_confirmed)*2)
                variants_unconfirmed_downsampled = resample(variants_unconfirmed, replace=False, n_samples=len(variants_confirmed)*2)

                self.variants = pd.concat([variants_confirmed_upsampled, variants_unconfirmed_downsampled], ignore_index=True)


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
            self.variants['qual_' + str(i+1) + '_caller_support'] = self.variants['qual_' + str(i+1) + '_caller_support'].apply(lambda x: np.round(x, 2))

    def set_features(self, columns):
        """ Set features. 
        
        param columns: list of strings, column names """

        features = []
        if 'variant' in self.clfparams['features']:
            features += ['size']
        if 'alignment' in self.clfparams['features']:
            features += [f for f in columns if f.startswith('ill')]
        if 'reference' in self.clfparams['features']:
            features += [f for f in columns if f.startswith('rep')] + ['cpg_islands', 'centromeres', 'asmb_gaps', 'alt_haps', 'GC_content_left', 'GC_content_right']
        
        return features


    def prepare_data_for_model(self, variants, features, labels=False):

        if self.clfparams['classifier'] == 'RandomForestClassifier':
            
            # Drop variants with missing values
            variants = variants.dropna(subset=features).copy().reset_index(drop=True)
            self.num_svs_rm_na = len(variants)

        elif self.clfparams['classifier'] == 'XGBoostClassifier':

            self.num_svs_rm_na = 0

        X = variants[features]
        if labels:
            y = variants['confirmed']
            return X, y
        else:
            return X


    def train(self):
        """ Train model. """

        # Filter based on excluded chromosomes
        variants_training = self.variants[~self.variants['chrom'].isin(self.chr_excl)].copy().reset_index(drop=True)

        columns = list(variants_training.columns[12:-1]) + ['size']
        features = self.set_features(columns)

        if self.clfparams['classifier'] == 'RandomForestClassifier':

            X, y = self.prepare_data_for_model(variants_training, features, labels=True)
            self.model = RandomForestClassifier(**self.clfparams['parameters'])
            self.model.fit(X, y)

        elif self.clfparams['classifier'] == 'XGBoostClassifier':
            
            X, y = self.prepare_data_for_model(variants_training, features, labels=True)
            self.model = xgb.XGBClassifier(**self.clfparams['parameters'])
            self.model.fit(X, y)

        else:

            raise ValueError(f'Invalid classifier: {self.clfparams["classifier"]}')


    def impute(self):

        # Impute size, for now only relevant for unresolved insertions by manta
        manta_ins_with_size = self.variants[(self.variants['method'] == 'manta') & (self.variants['type'] == 'INS') & (self.variants['size'].notna())].copy().reset_index(drop=True)
        manta_ins_without_size = self.variants[(self.variants['method'] == 'manta') & (self.variants['type'] == 'INS') & (self.variants['size'].isna())].copy().reset_index(drop=True)
        manta_ins_without_size['size'] = manta_ins_with_size['size'].median()
        self.variants.loc[self.variants['id'].isin(manta_ins_without_size['id']), 'size'] = manta_ins_without_size['size']

        # Impute GC Content
        self.variants['GC_content_left'] = self.variants['GC_content_left'].fillna(self.variants['GC_content_left'].median())
        self.variants['GC_content_right'] = self.variants['GC_content_right'].fillna(self.variants['GC_content_right'].median())

        # Impute coverage, was set to NA during feature extraction process if coverage exceeded threshold
        self.variants['ill_cov_mean_I'] = self.variants['ill_cov_mean_I'].fillna(self.cov_thr)
        self.variants['ill_cov_mean_II'] = self.variants['ill_cov_mean_II'].fillna(self.cov_thr)
        self.variants['ill_cov_mean_III'] = self.variants['ill_cov_mean_III'].fillna(self.cov_thr)
        self.variants['ill_cov_mean_IV'] = self.variants['ill_cov_mean_IV'].fillna(self.cov_thr)
        self.variants['ill_cov_std_I'] = self.variants['ill_cov_std_I'].fillna(1)
        self.variants['ill_cov_std_II'] = self.variants['ill_cov_std_II'].fillna(1)
        self.variants['ill_cov_std_III'] = self.variants['ill_cov_std_III'].fillna(1)
        self.variants['ill_cov_std_IV'] = self.variants['ill_cov_std_IV'].fillna(1)


    def predict(self):
        """ Predict probability of being TP for set of variants. """

        # Filter based on included chromosomes
        if self.chr_incl[0] != 'all':
            variants_predicting = self.variants[self.variants['chrom'].isin(self.chr_incl)].copy().reset_index(drop=True)
        else:
            variants_predicting = self.variants.copy().reset_index(drop=True)

        columns = list(variants_predicting.columns[12:]) + ['size']
        features = self.set_features(columns)
        

        # Remove confirmed column if present, this is the case during model curation
        if 'confirmed' in features:
            features.remove('confirmed')

        if self.clfparams['classifier'] == 'RandomForestClassifier':

            X = self.prepare_data_for_model(variants_predicting, features)
            
            variants_na = variants_predicting[variants_predicting[features].isna().any(axis=1)].copy().reset_index(drop=True)
            variants_predicting = variants_predicting.dropna(subset=features).copy().reset_index(drop=True)

            variants_predicting['pred_dicast'] = self.model.predict(X)
            variants_predicting['qual_dicast'] = self.model.predict_proba(X)[:, 1]
            variants_predicting['qual_dicast'] = variants_predicting['qual_dicast'].apply(lambda x: np.round(x, 2))

            variants_na = variants_predicting[variants_predicting[features].isna().any(axis=1)].copy().reset_index(drop=True)
            variants_na['pred_dicast'] = 0
            variants_na['qual_dicast'] = 0

            self.variants = pd.concat([variants_predicting, variants_na], ignore_index=True)

        elif self.clfparams['classifier'] == 'XGBoostClassifier':

            X = self.prepare_data_for_model(variants_predicting, features)

            variants_predicting['pred_dicast'] = self.model.predict(X)
            variants_predicting['qual_dicast'] = self.model.predict_proba(X)[:, 1]
            variants_predicting['qual_dicast'] = variants_predicting['qual_dicast'].apply(lambda x: np.round(x, 2))

            self.variants = variants_predicting.copy().reset_index(drop=True)

        
    def get_curation_set(self):
        """ Save curation SVs. """

        qual_cols = [col for col in self.variants.columns if col.startswith('qual_')]
        variants_curation = self.variants[['sample', 'tech', 'method', 'id', 'type', 'chrom', 'chrom2', 
                                           'start', 'end', 'size', 'filter', 'confirmed', 
                                           'pred_dicast'] + qual_cols].copy().reset_index(drop=True)

        # Determine which variants are sent for curation
        variants_curation_fp = variants_curation[(variants_curation['confirmed'] == 0) & (variants_curation['qual_dicast'] > 0.4)].copy().reset_index(drop=True)
        variants_curation_fn = variants_curation[(variants_curation['confirmed'] == 1) & (variants_curation['qual_dicast'] < 0.4)].copy().reset_index(drop=True)
        variants_curation_fp['err_type'] = 'FP'
        variants_curation_fn['err_type'] = 'FN'
        variants_curation = pd.concat([variants_curation_fp, variants_curation_fn], ignore_index=True)

        logging.info(f'# Sending {len(variants_curation_fp)} FP and {len(variants_curation_fn)} FN and SVs from {self.chr_incl[0]} for manual curation')
        
        return variants_curation


    def save_test(self):
        """ Save test result. """

        model_dir = '/'.join(self.pkl.split('/')[:-1])
        eval_dir = model_dir + '/eval'

        if not os.path.exists(eval_dir):
            os.makedirs(eval_dir)

        test_file = eval_dir + '/' + self.clf + '_' + self.type + '_eval.tsv'
        self.variants['model_dicast'] = self.clf
        self.variants.to_csv(test_file, sep='\t', index=False, na_rep='NA')


    def load_model(self):
        """ Load model. """

        with open(self.pkl, 'rb') as f:
            self.model = pickle.load(f)


    def save_model(self):
        """ Save model. """
        
        print(self.pkl)

        # Create model directory if it does not exist
        model_dir = '/'.join(self.pkl.split('/')[:-1])
        if not os.path.exists(model_dir):
            os.makedirs(model_dir)

        # Save metadata in JSON file
        json_file = model_dir + '/' + self.clf + '_metadata.json'
        if os.path.exists(json_file):
            meta = json.load(open(json_file))
        else:
            meta = {}

        meta[self.type] = {}
        meta[self.type]['samples'] = self.variants['sample'].unique().tolist()
        meta[self.type]['techs'] = self.variants['tech'].unique().tolist()
        meta[self.type]['methods'] = self.variants['method'].unique().tolist()
        meta[self.type]['num_svs'] = len(self.variants)
        meta[self.type]['num_nas'] = len(self.variants) - self.num_svs_rm_na
        meta[self.type]['num_pos'] = len(self.variants[self.variants['confirmed'] == 1])
        meta[self.type]['num_neg'] = len(self.variants[self.variants['confirmed'] == 0])
        meta[self.type]['curation'] = self.incl_cur
        meta[self.type]['balance'] = self.balance
        meta[self.type]['chr_excl'] = self.chr_excl
        meta[self.type]['features'] = self.clfparams['features']
        meta[self.type]['params'] = self.clfparams['parameters']
        json_object = json.dumps(meta, indent=4)
        with open(json_file, 'w') as f:
            f.write(json_object)

        # Save model
        with open(self.pkl, 'wb') as f:
            pickle.dump(self.model, f)