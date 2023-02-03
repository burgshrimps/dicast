import pandas as pd 
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
import pickle
import logging

from lib.utils import replace_filename

class Dicast:
    """ Class to train a model. """

    def __init__(self, mode, svtype, pkl, params, clf=None, chr_excl=[], chr_incl=['all']):
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
        self.clf = clf
        self.pkl = pkl
        self.params = params
        self.chr_excl = chr_excl
        self.chr_incl = chr_incl


    def load_sample(self, sample, ref, variant_features, variant_labels=None):
        """ Load sample data. 
        
        param sample: string, sample name 
        param ref: string, reference genome name 
        param variant_features: string, variant feature file
        param variant_labels: string, variant label file
        
        return pandas dataframe with sample data """

        # Load variant features
        filename_features = replace_filename(variant_features, sample, ref)
        df_features = pd.read_csv(filename_features, sep='\t')

        # Load variant labels
        if variant_labels != None:
            filename_labels = replace_filename(variant_labels, sample, ref)
            df_labels = pd.read_csv(filename_labels)

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


    def filter_data(self):
        """ Filter data. """
        
        # Filter based on excluded chromosomes
        self.variants = self.variants[~self.variants['chrom'].isin(self.chr_excl)].copy().reset_index(drop=True)

        # Filter based on included chromosomes
        if self.chr_incl[0] != 'all':
            self.variants = self.variants[self.variants['chrom'].isin(self.chr_incl)].copy().reset_index(drop=True)

        # Filter based on SV type
        self.variants = self.variants[self.variants['type'] == self.type].copy().reset_index(drop=True)


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


    def merge_overlapping_svs(result, method_dfs, sample, methods):
        """ Checks which SV is called by multiple methods and adds the quality score of the other method to the result dataframe. 
        
        param result: pandas dataframe, result dataframe
        param method_dfs: dictionary, method dataframes
        param sample: string, sample name
        param methods: list, method names
        
        return pandas dataframe with updated quality scores """

        result_sample = result[result['sample'] == sample].copy().reset_index(drop=True)
        for i in range(len(methods)):
            for j in range(i+1, len(methods)):
                method_df_sample_1 = method_dfs[methods[i]][method_dfs[methods[i]]['sample'] == sample].copy().reset_index(drop=True)
                method_df_sample_2 = method_dfs[methods[j]][method_dfs[methods[j]]['sample'] == sample].copy().reset_index(drop=True)
                overlap_ids = extract_overlap_ids(method_df_sample_1, method_df_sample_2)
                for k in range(len(overlap_ids)):
                    result_sample.loc[result_sample['id'] == overlap_ids['id_1'][k], 'qual_' + methods[j]] = result_sample.loc[result_sample['id'] == overlap_ids['id_2'][k], 'qual_' + methods[j]].values
                    result_sample.loc[result_sample['id'] == overlap_ids['id_2'][k], 'qual_' + methods[i]] = result_sample.loc[result_sample['id'] == overlap_ids['id_1'][k], 'qual_' + methods[i]].values
                
        # Normalize Qualities
        for method in methods:
            result_sample['qual_' + method] = result_sample['qual_' + method] / result_sample['qual_' + method].max()

        return result_sample
    

    def train(self):
        """ Train model. """

        features = list(self.variants.columns[12:-1]) + ['size']
        num_svs = len(self.variants)
        self.variants = self.variants.dropna(subset=features).copy().reset_index(drop=True)
        logging.info(f'Dropped {num_svs - len(self.variants)}/{num_svs} SVs due to missing feature values')

        X = self.variants[features]
        y = self.variants['confirmed']

        if self.clf == 'rf':
            self.model = RandomForestClassifier(n_estimators=100)
        else:
            raise ValueError(f'Invalid classifier: {self.clf}')

        self.model.fit(X, y)


    def predict(self):
        """ Predict probability of being TP for set of variants. """

        features = list(self.variants.columns[12:-1]) + ['size']
        num_svs = len(self.variants)
        self.variants = self.variants.dropna(subset=features).copy().reset_index(drop=True)
        logging.info(f'Dropped {num_svs - len(self.variants)}/{num_svs} SVs due to missing feature values')

        X = self.variants[features]
        y = self.variants['confirmed']

        self.variants['dicast_pred'] = self.model.predict(X)
        self.variants['dicast_qual'] = self.model.predict_proba(X)[:, 1]

    
    def save_predictions_tsv(self, out_file):
        """ Save predictions. """

        self.variants[['id', 'sample', 'tech', 'method', 'type', 'chrom', 'chrom2', 'start', 'end', 'size', 
                       'filter', 'qual', 'dicast_pred', 'dicast_qual']].to_csv(out_file, sep='\t', index=False, na_rep='NA')


    def load_model(self):
        """ Load model. """

        with open(self.pkl, 'rb') as f:
            self.model = pickle.load(f)
        
        
    def save_model(self):
        """ Save model. """

        with open(self.pkl, 'wb') as f:
            pickle.dump(self.model, f)