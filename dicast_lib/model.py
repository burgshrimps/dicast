import pandas as pd 
import xgboost as xgb
import pickle
import json
import datetime
import numpy as np


class Dicast:
    """ Structural variant detection from short-read sequencing data. """

    def __init__(self, sv_type: str):
        """ Initialize Dicast object. """
        
        self.sv_type = sv_type
        self.cov_thr = 6 # Log2 threshold were feature collection was aborted
        
        # Determine set of features to use
        self.features_var = ['sv_len']
        self.features_ref = ['rep_LINE', 'rep_SINE', 'rep_LTR', 'rep_DNA', 'rep_Simple_repeat', 'rep_Satellite', 'rep_Low_complexity',
                             'rep_Retroposon', 'rep_snRNA', 'rep_tRNA', 'rep_srpRNA', 'rep_rRNA','rep_RC', 'rep_scRNA', 'rep_RNA', 'rep_VNTR', 
                             'rep_STR', 'cpg_islands', 'centromeres', 'asmb_gaps', 'alt_haps', 'GC_content_left', 'GC_content_right']
        self.features_aln = {'DEL': ['ill_cov_mean', 'ill_cov_std', 'ill_isize_mean', 'ill_isize_std', 'ill_mapq_mean', 'ill_mapq_std',
                                     'ill_clipreads', 'ill_splitreads', 'ill_disco_ff', 'ill_disco_rr'],
                             'INS': ['ill_cov_mean', 'ill_cov_std', 'ill_isize_mean', 'ill_isize_std', 'ill_mapq_mean', 'ill_mapq_std',
                                     'ill_clipreads', 'ill_splitreads', 'ill_disco_ff', 'ill_disco_rr'],
                             'INV': ['ill_cov_mean', 'ill_cov_std', 'ill_isize_mean', 'ill_isize_std', 'ill_mapq_mean', 'ill_mapq_std',
                                     'ill_clipreads', 'ill_splitreads', 'ill_disco_ff', 'ill_disco_rr'],
                             'DUP': ['ill_cov_mean', 'ill_cov_std', 'ill_isize_mean', 'ill_isize_std', 'ill_mapq_mean', 'ill_mapq_std',
                                     'ill_clipreads', 'ill_splitreads', 'ill_disco_ff', 'ill_disco_rr'],
                             'BND': ['ill_cov_mean', 'ill_cov_std', 'ill_isize_mean', 'ill_isize_std', 'ill_mapq_mean', 'ill_mapq_std',
                                     'ill_clipreads', 'ill_splitreads', 'ill_disco_ff', 'ill_disco_rr']}
        suffices = ['I', 'II'] if self.sv_type == 'INS' else ['I', 'II', 'III', 'IV']
        self.features = self.features_var + self.features_ref + [aln_feature + '_' + suffix for aln_feature in self.features_aln[self.sv_type] for suffix in suffices]
        
    
    def load_from_df(self, df_variants: pd.DataFrame):
        """ Load variants from dataframe.

        Args:
            df_variants (pd.DataFrame): Dataframe containing variants with features
        """        
        
        self.variants = df_variants
        
        
    def load_from_csv(self, variant_filename: str):
    
        self.variants = pd.read_csv(variant_filename, sep='\t', low_memory=False)
        
        
    def impute_missing_values(self):
        """ Impute missing values. """
        
        # Impute GC Content
        self.variants['GC_content_left'] = self.variants['GC_content_left'].fillna(self.variants['GC_content_left'].median())
        self.variants['GC_content_right'] = self.variants['GC_content_right'].fillna(self.variants['GC_content_right'].median())

        # Impute coverage, was set to NA during feature extraction process if coverage exceeded threshold
        self.variants['ill_cov_mean_I'] = self.variants['ill_cov_mean_I'].fillna(self.cov_thr)
        self.variants['ill_cov_mean_II'] = self.variants['ill_cov_mean_II'].fillna(self.cov_thr)
        self.variants['ill_cov_std_I'] = self.variants['ill_cov_std_I'].fillna(1)
        self.variants['ill_cov_std_II'] = self.variants['ill_cov_std_II'].fillna(1)
        
        if self.sv_type != 'INS':
            self.variants['ill_cov_mean_III'] = self.variants['ill_cov_mean_III'].fillna(self.cov_thr)
            self.variants['ill_cov_mean_IV'] = self.variants['ill_cov_mean_IV'].fillna(self.cov_thr)
            self.variants['ill_cov_std_III'] = self.variants['ill_cov_std_III'].fillna(1)
            self.variants['ill_cov_std_IV'] = self.variants['ill_cov_std_IV'].fillna(1)
            
            
    def train(self, model_type: str, model_params: dict, chroms: list=[]):
        """ Train model.

        Args:
            model_type (str): Type of model to train
            model_params (dict): Parameters for model
            chroms (list, optional): List of chromosomes to use for training. Defaults to [].
        """           
        
        self.chroms_train = chroms
        self.model_type = model_type
        self.model_params = model_params
        
        # Subset variants to chromosomes used for training
        if len(self.chroms_train) > 0:
            self.variants_train = self.variants[self.variants['chrom'].isin(self.chroms_train)].copy().reset_index(drop=True)
        else:
            self.variants_train = self.variants.copy()
        
        # Train model
        X = self.variants_train[self.features]
        y = self.variants_train['confirmation_status']
        
        if model_type == 'XGBoost':
            self.model = xgb.XGBClassifier(**self.model_params)
            self.model.fit(X, y)
            
            
    def predict(self, chroms: list=[]):
        """ Make predictions for variants.

        Args:
            chroms (list, optional): List of chromosomes to use for predicting. Defaults to [].
        """        
        
        self.chroms_predict = chroms
        
        # Subset variants to chromosomes used for prediction
        if len(self.chroms_predict) > 0:
            self.variants_predict = self.variants[self.variants['chrom'].isin(self.chroms_predict)].copy().reset_index(drop=True)
        else:
            self.variants_predict = self.variants.copy()
            
        # Predict
        X = self.variants_predict[self.features]
        
        self.variants_predict['dicast_qual'] = np.round(self.model.predict_proba(X)[:, 1], 3)
        
        
    def save_model(self, model_filename: str):
        """ Save model to file.

        Args:
            model_filename (str): Filename to save the model to, must end with .pkl
        """        
        
        # Save model
        with open(model_filename, 'wb') as f:
            pickle.dump(self.model, f)
            
        # Save model metadata
        meta_filename = model_filename.replace('.pkl', '_metadata.json')
        meta = {'date' : datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                'model_type': self.model_type,
                'model_params': self.model_params,
                'features': self.features,
                'sv_type': self.sv_type,
                'reference': self.variants['reference'].unique().tolist()[0],
                'cohorts': self.variants['cohort'].unique().tolist(),
                'samples': self.variants['sample'].unique().tolist(),
                'callers': self.variants['caller'].unique().tolist(),
                'chroms_train': self.chroms_train,
                'number_variants': self.variants.shape[0],
                'number_variants_positive': self.variants[self.variants['confirmation_status'] == 1].shape[0],
                'number_variants_negative': self.variants[self.variants['confirmation_status'] == 0].shape[0]}
        json_object = json.dumps(meta, indent=4)
        with open(meta_filename, 'w') as f:
            f.write(json_object)
            
            
    def load_model(self, model_filename: str):
        """ Load model from file.

        Args:
            model_filename (str): Filename to load the model from, must end with .pkl
        """        
        
        with open(model_filename, 'rb') as f:
            self.model = pickle.load(f)
            
            
    def to_database(self) -> pd.DataFrame:
        """ Get predictions for variants.

        Returns:
            pd.DataFrame: Dataframe containing variants with predictions
        """        
        
        columns_for_export = ['single_id', 'merged_id', 'caller_id', 'cohort', 'sample', 'reference', 'technology', 'caller', 'sv_type', 
                              'chrom', 'chrom_2', 'start', 'end', 'sv_len', 'filter', 'caller_qual', 'dicast_qual', 'genotype',
                              'performed_confirmation', 'confirmation_status', 'performed_curation', 'curation_status']
        return self.variants_predict[columns_for_export]
    
    
    def to_df(self) -> pd.DataFrame:
        """ Get predictions for variants.

        Returns:
            pd.DataFrame: Dataframe containing variants with predictions
        """        
        
        columns_for_export = ['id', 'cohort', 'sample', 'reference', 'technology', 'caller', 'sv_type', 
                              'chrom', 'chrom_2', 'start', 'end', 'sv_len', 'filter', 'qual', 'dicast_qual', 'genotype']
        return self.variants_predict[columns_for_export]