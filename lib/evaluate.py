import pysam 
import pandas as pd
import os
import numpy as np
import re
import bioframe as bf
from sklearn.metrics import roc_auc_score, precision_recall_curve, roc_curve
import plotly.express as px


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
        self.svtypes = ['DEL', 'INS']


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
                vcf_dfs.append(df)
        
        self.variants_methods = pd.concat(vcf_dfs, ignore_index=True)


    def read_benchmark_variants(self):
        """ Reads benchmark file and saves it in pandas dataframe. """

        self.variants_bench = pd.read_csv(self.fname_benchmark, sep='\t')
        self.variants_bench = self.variants_bench[['ID', '#CHROM', 'POS', 'END', 'SVTYPE', 'SVLEN', 'MERGE_SAMPLES']].copy()
        self.variants_bench.columns = ['id', 'chrom', 'start', 'end', 'type', 'size', 'samples']
        self.variants_bench = self.variants_bench[self.variants_bench['samples'].str.contains(self.sample)].copy().reset_index(drop=True)
        self.variants_bench = self.variants_bench[self.variants_bench['chrom'].isin(CHROMS)].copy().reset_index(drop=True)
        self.variants_bench.drop('samples', axis=1, inplace=True)


    def read_dicast_variants(self):
        """ Reads dicast file and saves it in pandas dataframe. """

        self.variants_dicast = pd.read_csv(self.fname_dicast, sep='\t')
        self.variants_dicast = self.variants_dicast[self.variants_dicast['type'].isin(self.svtypes)].copy().reset_index(drop=True)
        self.variants_dicast['tech'] = 'mgi'
        self.variants_dicast.rename(columns={'id': 'id_org'}, inplace=True)

        dicast_variant_dfs = []
        for svtype in self.svtypes:
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


    def overlap_benchmark(self):
        
        overlap_ids_methods = []
        overlap_ids_dicast = []
        for tech in self.fnames_methods:
            for method in self.fnames_methods[tech]:
                for svtype in self.svtypes:
                    curr_method_df = self.variants_methods.loc[(self.variants_methods['tech'] == tech) & (self.variants_methods['method'] == method) & (self.variants_methods['type'] == svtype), ['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)
                    curr_dicast_df = self.variants_dicast.loc[(self.variants_dicast['tech'] == tech) & (self.variants_dicast['method'] == method) & (self.variants_dicast['type'] == svtype), ['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)
                    curr_benchmark_df = self.variants_bench.loc[self.variants_bench['type'] == svtype, ['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)


                    if not curr_method_df.empty and not curr_benchmark_df.empty:
                        curr_method_df['end'] = curr_method_df['end'].astype(int)
                        curr_dicast_df['end'] = curr_dicast_df['end'].astype(int)
                        curr_benchmark_df['end'] = curr_benchmark_df['end'].astype(int)

                        curr_overlap_ids_method = self.extract_overlap_ids(curr_benchmark_df, curr_method_df)
                        curr_overlap_ids_dicast = self.extract_overlap_ids(curr_benchmark_df, curr_dicast_df)

                        curr_overlap_ids_method['tech'] = tech
                        curr_overlap_ids_method['method'] = method
                        curr_overlap_ids_method['type'] = svtype
                        curr_overlap_ids_dicast['tech'] = tech
                        curr_overlap_ids_dicast['method'] = method
                        curr_overlap_ids_dicast['type'] = svtype

                        overlap_ids_methods.append(curr_overlap_ids_method)
                        overlap_ids_dicast.append(curr_overlap_ids_dicast)
        self.overlap_ids_methods = pd.concat(overlap_ids_methods, ignore_index=True)
        self.overlap_ids_dicast = pd.concat(overlap_ids_dicast, ignore_index=True)


    def confirm_variants(self):
        
        method_dfs = []
        dicast_dfs = []
        for tech in self.fnames_methods:
            for method in self.fnames_methods[tech]:
                for svtype in self.svtypes:
                    curr_method_df = self.variants_methods[(self.variants_methods['tech'] == tech) & (self.variants_methods['method'] == method) & (self.variants_methods['type'] == svtype)].copy().reset_index(drop=True)
                    curr_dicast_df = self.variants_dicast[(self.variants_dicast['tech'] == tech) & (self.variants_dicast['method'] == method) & (self.variants_dicast['type'] == svtype)].copy().reset_index(drop=True)
                    curr_overlap_ids_methods = self.overlap_ids_methods[(self.overlap_ids_methods['tech'] == tech) & (self.overlap_ids_methods['method'] == method) & (self.overlap_ids_methods['type'] == svtype)].copy().reset_index(drop=True)
                    curr_overlap_ids_dicast = self.overlap_ids_dicast[(self.overlap_ids_dicast['tech'] == tech) & (self.overlap_ids_dicast['method'] == method) & (self.overlap_ids_dicast['type'] == svtype)].copy().reset_index(drop=True)

                    curr_method_df['confirmed'] = 0
                    curr_method_df.loc[curr_method_df['id'].isin(curr_overlap_ids_methods['id_2']), 'confirmed'] = 1
                    method_dfs.append(curr_method_df)

                    curr_dicast_df['confirmed'] = 0
                    curr_dicast_df.loc[curr_dicast_df['id'].isin(curr_overlap_ids_dicast['id_2']), 'confirmed'] = 1
                    dicast_dfs.append(curr_dicast_df)

        self.variants_methods = pd.concat(method_dfs, ignore_index=True)
        self.variants_dicast = pd.concat(dicast_dfs, ignore_index=True)

        self.variants = dict()

        # Dicast
        self.variants['dicast'] = self.variants_dicast[['id', 'qual_dicast', 'confirmed', 'type']].rename(columns={'qual_dicast': 'qual'})
        variants_missed_dicast = self.variants_bench.loc[~self.variants_bench['id'].isin(self.overlap_ids_dicast['id_1']), ['id', 'type']]
        variants_missed_dicast['confirmed'] = 1
        variants_missed_dicast['qual'] = 0
        self.variants['dicast'] = pd.concat([self.variants['dicast'], variants_missed_dicast], ignore_index=True)

        # Methods
        for tech in self.fnames_methods:
            for method in self.fnames_methods[tech]:
                curr_variants = self.variants_methods.loc[(self.variants_methods['tech'] == tech) & (self.variants_methods['method'] == method), ['id', 'qual', 'confirmed', 'type']].copy()
                self.variants[tech + '_' + method] = curr_variants

                variants_missed = self.variants_bench.loc[~self.variants_bench['id'].isin(self.overlap_ids_methods.loc[(self.overlap_ids_methods['tech'] == tech) & (self.overlap_ids_methods['method'] == method), 'id_1']), ['id', 'type']]
                variants_missed['confirmed'] = 1
                variants_missed['qual'] = 0
                self.variants[tech + '_' + method] = pd.concat([self.variants[tech + '_' + method], variants_missed], ignore_index=True)


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


    def plot_precision_recall_curve(self, svtype):
        """ Plot precision recall curve. 

        param pr_rc_df: pandas dataframe with precision and recall for different thresholds """

        pr_rc_df = self.precision_recall_df[self.precision_recall_df['type'] == svtype].copy().reset_index(drop=True)

        colors = ['black', '#1f77b4', '#ff7f0e', 'darkred', '#1f77b4', '#ff7f0e', 'darkred']
        dash = ['solid', 'solid', 'solid', 'solid', 'dot', 'dot', 'dot']
        circle_bg_white = [0, 0, 0, 1, 1, 1, 1]
        fig = px.line(x='recall', y='precision', color='method', 
                    data_frame=pr_rc_df, line_dash='method', 
                    line_dash_sequence=dash,
                    color_discrete_sequence=colors)


        for i, method in enumerate(self.variants.keys()):
            x = pr_rc_df[pr_rc_df['method'] == method].reset_index(drop=True).loc[0, 'recall']
            y = pr_rc_df[pr_rc_df['method'] == method].reset_index(drop=True).loc[0, 'precision']
            if circle_bg_white[i] == 1:
                fig.add_shape(type='circle', xref='x', yref='y', x0=x-0.002, y0=y-0.015, x1=x+0.003, y1=y+0.005, line_color=colors[i], line_width=2, opacity=1, fillcolor='white')
            else:
                fig.add_shape(type='circle', xref='x', yref='y', x0=x-0.002, y0=y-0.015, x1=x+0.003, y1=y+0.005, line_color=colors[i], line_width=2, opacity=1, fillcolor=colors[i])


        fig.update_layout(plot_bgcolor='white', xaxis_title='Recall', yaxis_title='Precision', xaxis_linecolor='black', yaxis_linecolor='black')
        fig.update_traces(line=dict(width=2))
        fig.update_xaxes(ticks='outside', tickcolor='black', tickwidth=1, ticklen=5, gridcolor='lightgray', gridwidth=0.5, range=[0, 1.1])
        fig.update_yaxes(ticks='outside', tickcolor='black', tickwidth=1, ticklen=5, gridcolor='lightgray', gridwidth=0.5, range=[0, 1.1])
        fig.show()

        

                        

    