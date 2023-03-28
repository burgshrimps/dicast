import pysam 
import pandas as pd
import os
import numpy as np
import re
import bioframe as bf
from sklearn.metrics import roc_auc_score, precision_recall_curve, roc_curve
import plotly.express as px
from tqdm import tqdm


from lib.utils import replace_filename, parse_vcf

CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 
          'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
          'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']

class Eva:
    """ Class to evaluate dicast compared to other methods. """


    def __init__(self, params):
        self.sample = params['sample']
        self.ref = params['ref']
        self.fname_dicast = params['dicast']
        self.fname_benchmark = params['benchmark']
        self.fnames_methods = params['vcf']
        
        self.max_dist_overlap = 500
        self.min_size_overlap = 0.7
        self.svtypes = ['DEL']
        self.variants = dict()


    def read_vcf(self, filename):
        """ Reads VCF file as a pysam.VariantFile object and takes care of sample names. """

        if os.path.exists(filename):
            vcf = pysam.VariantFile(filename)
        else:
            # Smoove appears in lumpy filenames
            vcf = pysam.VariantFile(filename.replace('-', '_').replace('_smoove', '-smoove'))

        return vcf

    
    def read_method_variants(self):
        """ Reads all VCFs in params file and saves them in pandas dataframe. """

        vcf_dfs = []
        for tech in self.fnames_methods:
            for method in self.fnames_methods[tech]:
                filename_vcf = replace_filename(self.fnames_methods[tech][method], self.sample, self.ref)
                vcf = self.read_vcf(filename_vcf)
                df = parse_vcf(vcf, tech, method, self.sample)
                df = df[df['chrom'].isin(CHROMS)].copy().reset_index(drop=True)
                vcf_dfs.append(df)
        
        self.variants_methods = pd.concat(vcf_dfs, ignore_index=True)
        self.variants_methods = self.variants_methods[self.variants_methods['type'].isin(self.svtypes)].copy().reset_index(drop=True)

        # Set quality of cvnator to 1 because it does not have a quality score
        self.variants_methods.loc[self.variants_methods['method'] == 'cnvnator', 'qual'] = 1



    def read_benchmark_variants(self):
        """ Reads benchmark file and saves it in pandas dataframe. """

        self.variants_bench = pd.read_csv(self.fname_benchmark, sep='\t')
        self.variants_bench = self.variants_bench[['ID', '#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'MERGE_SAMPLES']].copy()
        self.variants_bench.columns = ['id', 'chrom', 'start', 'end', 'type', 'size', 'samples']
        self.variants_bench = self.variants_bench[self.variants_bench['samples'].str.contains(self.sample)].copy().reset_index(drop=True)
        self.variants_bench = self.variants_bench[self.variants_bench['chrom'].isin(CHROMS)].copy().reset_index(drop=True)
        self.variants_bench.drop('samples', axis=1, inplace=True)
        self.variants_bench['confirmed'] = 1

        self.variants_bench = self.variants_bench[self.variants_bench['type'].isin(self.svtypes)].copy().reset_index(drop=True)


    def read_dicast_variants(self):
        """ Reads dicast file and saves it in pandas dataframe. """

        self.variants_dicast = pd.read_csv(self.fname_dicast, sep='\t')
        self.variants_dicast = self.variants_dicast[self.variants_dicast['type'].isin(self.svtypes)].copy().reset_index(drop=True)
        self.variants_dicast.rename(columns={'id': 'id_org'}, inplace=True)

        # create unique IDs for dicast variants 
        dicast_variant_dfs = []
        for svtype in self.svtypes:

            # create unique id
            dicast_variants_type = self.variants_dicast.loc[self.variants_dicast['type'] == svtype].copy().reset_index(drop=True)
            dicast_variants_type['index'] = dicast_variants_type.index + 1
            dicast_variants_type['id'] = 'dicast.' + svtype + '.' + dicast_variants_type['index'].astype(str)
            dicast_variant_dfs.append(dicast_variants_type)

        self.variants_dicast = pd.concat(dicast_variant_dfs, ignore_index=True)
        self.variants_dicast.drop('index', axis=1, inplace=True)


    def extract_overlap_ids(self, df1, df2):
        """ Extracts SV IDs of overlapping variants
        
        param df1: pandas dataframe, first dataframe
        param df2: pandas dataframe, second dataframe 
        
        return pandas dataframe with overlapping SV IDs """

        df1['end'] = df1['end'].astype(int)
        df2['end'] = df2['end'].astype(int)

        closest_intervals = bf.closest(df1, df2, suffixes=('_1','_2'), k=10) # k=25 to get all overlapping variants
        closest_intervals = closest_intervals.dropna(subset=['id_1', 'id_2']).reset_index(drop=True)
        closest_intervals['diff_start'] = abs(closest_intervals['start_1'] - closest_intervals['start_2'])
        closest_intervals['diff_end'] = abs(closest_intervals['end_1'] - closest_intervals['end_2'])
        closest_intervals['diff_size'] = closest_intervals.apply(lambda x: min([x['size_1'], x['size_2']]) / max([x['size_1'], x['size_2']]), axis=1)

        overlapping_svs = closest_intervals[(closest_intervals['diff_start'] < self.max_dist_overlap) & (closest_intervals['diff_end'] < self.max_dist_overlap) & (closest_intervals['diff_size'] > self.min_size_overlap)].copy()
        
        #overlapping_svs = closest_intervals

        return overlapping_svs.reset_index(drop=True)


    def overlap_bench_methods(self):

        for tech in self.fnames_methods:
            for method in self.fnames_methods[tech]:

                dfs = []
                for svtype in self.svtypes:

                    curr_bench_df = self.variants_bench[self.variants_bench['type'] == svtype].copy().reset_index(drop=True)
                    curr_method_df = self.variants_methods[(self.variants_methods['type'] == svtype) & (self.variants_methods['method'] == method)].copy().reset_index(drop=True)

                    overlap_bench_method = self.extract_overlap_ids(curr_bench_df, curr_method_df)
                    ids_bench = overlap_bench_method['id_1'].unique()
                    ids_method = overlap_bench_method['id_2'].unique()
                    benchmark_method = overlap_bench_method.groupby('id_1').agg({'chrom_1' : 'first', 'start_1' : 'first', 'end_1' : 'first',
                                                                                 'type_1' : 'first', 'size_1' : 'first', 'confirmed_1' : 'first',
                                                                                 'qual_2' : max}).reset_index()
                    benchmark_method.columns = ['id', 'chrom', 'start', 'end', 'type', 'size', 'confirmed', 'qual']

                    adding_bench = curr_bench_df[~curr_bench_df['id'].isin(ids_bench)].copy()
                    adding_bench['qual'] = 0

                    adding_method = curr_method_df[~curr_method_df['id'].isin(ids_method)].copy()
                    adding_method['confirmed'] = 0
                    adding_method = adding_method[['id', 'chrom', 'start', 'end', 'type', 'size', 'confirmed', 'qual']]

                    dfs.append(pd.concat([benchmark_method, adding_bench, adding_method], ignore_index=True))

                self.variants[method] = pd.concat(dfs, axis=0, ignore_index=True)



    def overlap_bench_dicast(self):

        dfs = []
        for svtype in self.svtypes:

            curr_bench_df = self.variants_bench[self.variants_bench['type'] == svtype].copy().reset_index(drop=True)
            curr_dicast_df = self.variants_dicast[self.variants_dicast['type'] == svtype].copy().reset_index(drop=True)

            overlap_bench_dicast = self.extract_overlap_ids(curr_bench_df, curr_dicast_df)
            ids_bench = overlap_bench_dicast['id_1'].unique()
            ids_dicast = overlap_bench_dicast['id_2'].unique() 
            benchmark_dicast = overlap_bench_dicast.groupby('id_1').agg({'chrom_1' : 'first', 'start_1' : 'first', 'end_1' : 'first', 
                                                                        'type_1' : 'first', 'size_1' : 'first', 'confirmed_1' : 'first',
                                                                        'qual_dicast_2' : max}).reset_index()
            benchmark_dicast.columns = ['id', 'chrom', 'start', 'end', 'type', 'size', 'confirmed', 'qual_dicast']

            adding_bench = curr_bench_df[~curr_bench_df['id'].isin(ids_bench)].copy()
            adding_bench['qual_dicast'] = 0

            adding_dicast = curr_dicast_df[~curr_dicast_df['id'].isin(ids_dicast)].copy()
            adding_dicast['confirmed'] = 0
            adding_dicast = adding_dicast[['id', 'chrom', 'start', 'end', 'type', 'size', 'confirmed', 'qual_dicast']]

            dfs.append(pd.concat([benchmark_dicast, adding_bench, adding_dicast], axis=0, ignore_index=True))
        
        self.variants['dicast'] = pd.concat(dfs, axis=0, ignore_index=True).rename(columns={'qual_dicast' : 'qual'})


    def compute_precision_recall_df(self):
        """ Compute precision and recall for different thresholds. 
        
        param df_eval: pandas dataframe with evaluation info
        
        return: pandas dataframe with precision and recall for different thresholds """

        pr_rc_dict = {'method' : [], 'precision' : [], 'recall' : [], 'type' : []}

        for method in self.variants:
            for svtype in self.svtypes:
                df = self.variants[method][self.variants[method]['type'] == svtype].copy().reset_index(drop=True)
                precision, recall, _ = precision_recall_curve(df['confirmed'], df['qual'])
                pr_rc_dict['method'].extend([method] * (len(precision) - 1))
                pr_rc_dict['type'].extend([svtype] * (len(precision) - 1))
                pr_rc_dict['precision'].extend(precision[1:])
                pr_rc_dict['recall'].extend(recall[1:])
        
        self.precision_recall_df = pd.DataFrame(pr_rc_dict)

        

                        

    