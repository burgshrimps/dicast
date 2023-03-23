import pysam 
import pandas as pd
import os
import numpy as np
import re
import bioframe as bf


from lib.utils import replace_filename, parse_vcf

class Eva:
    """ Class to evaluate dicast compared to other methods. """


    def __init__(self, params):
        self.sample = params['sample']
        self.ref = params['ref']
        self.dicast = params['dicast']
        self.benchmark = params['benchmark']
        self.vcf = params['vcf']
        
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
        for tech in self.vcf:
            for method in self.vcf[tech]:
                filename_vcf = replace_filename(self.vcf[tech][method], self.sample, self.ref)
                vcf = self.read_vcf(filename_vcf)
                df = parse_vcf(vcf, tech, method, self.sample)
                vcf_dfs.append(df)
        
        self.method_variants = pd.concat(vcf_dfs, ignore_index=True)


    def read_benchmark_variants(self):
        """ Reads benchmark file and saves it in pandas dataframe. """

        vcf = self.read_vcf(self.benchmark)
        self.benchmark_variants = self.parse_vcf(vcf, 'benchmark', 'benchmark', self.sample, check_gt=True)


    def read_dicast_variants(self):
        """ Reads dicast file and saves it in pandas dataframe. """

        self.dicast_variants = pd.read_csv(self.dicast, sep='\t')
        self.dicast_variants = self.dicast_variants[self.dicast_variants['type'].isin(self.svtypes)].copy().reset_index(drop=True)
        self.dicast_variants['tech'] = 'mgi'


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

        overlapping_svs = closest_intervals[(closest_intervals['diff_start'] < self.max_dist_overlap) & (closest_intervals['diff_end'] < self.max_dist_overlap) & (closest_intervals['diff_size'] > self.min_size_overlap)].copy()
        
        #overlapping_svs = closest_intervals

        return overlapping_svs.reset_index(drop=True)


    def overlap_benchmark(self):
        
        overlap_ids_methods = []
        overlap_ids_dicast = []
        for tech in self.vcf:
            for method in self.vcf[tech]:
                for svtype in self.svtypes:
                    curr_method_df = self.method_variants.loc[(self.method_variants['tech'] == tech) & (self.method_variants['method'] == method) & (self.method_variants['type'] == svtype), ['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)
                    curr_dicast_df = self.dicast_variants.loc[(self.dicast_variants['tech'] == tech) & (self.dicast_variants['method'] == method) & (self.dicast_variants['type'] == svtype), ['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)
                    curr_benchmark_df = self.benchmark_variants.loc[self.benchmark_variants['type'] == svtype, ['id', 'chrom', 'start', 'end', 'size']].copy().reset_index(drop=True)


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
        for tech in self.vcf:
            for method in self.vcf[tech]:
                for svtype in self.svtypes:
                    curr_method_df = self.method_variants[(self.method_variants['tech'] == tech) & (self.method_variants['method'] == method) & (self.method_variants['type'] == svtype)].copy().reset_index(drop=True)
                    curr_dicast_df = self.dicast_variants[(self.dicast_variants['tech'] == tech) & (self.dicast_variants['method'] == method) & (self.dicast_variants['type'] == svtype)].copy().reset_index(drop=True)
                    curr_overlap_ids_methods = self.overlap_ids_methods[(self.overlap_ids_methods['tech'] == tech) & (self.overlap_ids_methods['method'] == method) & (self.overlap_ids_methods['type'] == svtype)].copy().reset_index(drop=True)
                    curr_overlap_ids_dicast = self.overlap_ids_dicast[(self.overlap_ids_dicast['tech'] == tech) & (self.overlap_ids_dicast['method'] == method) & (self.overlap_ids_dicast['type'] == svtype)].copy().reset_index(drop=True)

                    curr_method_df['confirmed'] = 0
                    curr_method_df.loc[curr_method_df['id'].isin(curr_overlap_ids_methods['id_2']), 'confirmed'] = 1
                    method_dfs.append(curr_method_df)

                    curr_dicast_df['confirmed'] = 0
                    curr_dicast_df.loc[curr_dicast_df['id'].isin(curr_overlap_ids_dicast['id_2']), 'confirmed'] = 1
                    dicast_dfs.append(curr_dicast_df)

        self.method_variants = pd.concat(method_dfs, ignore_index=True)
        self.dicast_variants = pd.concat(dicast_dfs, ignore_index=True)
                        

    