import pandas as pd
import pysam
import numpy as np
import re
import os

from lib.utils import replace_filename


class VariantPrep:
    """ Class to prepare raw variant calls for feature extraction. """

    def __init__(self, sample, ref, params, workdir):
        """ Initialize class. """

        # Input parameters
        self.sample = sample
        self.ref = ref
        self.params = params
        self.workdir = workdir

        # Auxiallary files for preparation
        self.ref_dir = self.params['ref'][self.ref]['directory']
        self.filename_chrom_sizes = self.params['ref'][self.ref]['filename_chrom_sizes']
        self.chrom_sizes = pd.read_csv('/'.join([self.ref_dir, self.filename_chrom_sizes]), sep='\t',
                                       header=None, names=['size', 'offset', 'linebases', 'linewidth'], index_col=0)

        # Output files
        self.out_dir = '/'.join([self.workdir, 'ensemble'])
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        self.out_filename = f'{self.sample}_{self.ref}.SVs.raw.tsv'
        self.out_file = '/'.join([self.out_dir, self.out_filename])

        # List of chromosomes to use
        self.chroms = params['chroms']


    def parse_vcf(self, vcf, tech, method, sample):
        """ Parses VCF and saves info in pandas dataframe.

        param vcf: pysam.VariantFile object 
        param tech: string, technology used for SV calling
        param method: string, method used for SV calling
        param sample: string, sample name 
        
        return: pandas dataframe with SV info """  

        vcf_dict = {'id': [], 'sample': [], 'tech' : [], 'method' : [], 'type': [], 'chrom': [], 'start' : [], 'chrom2' : [], 'end': [], 'size' : [], 'filter': [], 'qual' : []}

        for rec in vcf.fetch():
            # Information that is always present
            vcf_dict['id'].append(rec.id)
            vcf_dict['sample'].append(sample)
            vcf_dict['tech'].append(tech)
            vcf_dict['method'].append(method)
            vcf_dict['type'].append(rec.info['SVTYPE'])
            vcf_dict['chrom'].append(rec.chrom)
            vcf_dict['start'].append(rec.start)
            vcf_dict['filter'].append(', '.join(rec.filter.keys()))
            vcf_dict['qual'].append(rec.qual)


            # Deletions
            if rec.info['SVTYPE'] == 'DEL':
                vcf_dict['end'].append(rec.stop)
                vcf_dict['size'].append(rec.stop - rec.start)
                vcf_dict['chrom2'].append(np.nan)
            
            # Insertions
            elif rec.info['SVTYPE'] == 'INS':
                if rec.stop == rec.start:
                    vcf_dict['end'].append(rec.stop + 1)
                else:
                    vcf_dict['end'].append(rec.stop)
                try:
                    if method == 'manta':
                        vcf_dict['size'].append(rec.info['SVLEN'][0])
                    else:
                        vcf_dict['size'].append(rec.info['SVLEN'])
                except KeyError:
                    vcf_dict['size'].append(np.nan)
                vcf_dict['chrom2'].append(np.nan)

            # Inversions
            elif rec.info['SVTYPE'] == 'INV':
                vcf_dict['end'].append(rec.stop)
                vcf_dict['size'].append(rec.stop - rec.start)
                vcf_dict['chrom2'].append(np.nan)

            # Duplications
            elif rec.info['SVTYPE'] == 'DUP':
                vcf_dict['end'].append(rec.stop)
                vcf_dict['size'].append(rec.stop - rec.start)
                vcf_dict['chrom2'].append(np.nan)

            # Breakends (Translocations)
            elif rec.info['SVTYPE'] == 'BND':
                vcf_dict['chrom2'].append(re.search(r'chr.*:', rec.alts[0]).group(0)[:-1])
                vcf_dict['end'].append(re.search(r':[0-9]*', rec.alts[0]).group(0)[1:])
                vcf_dict['size'].append(np.nan)

        return pd.DataFrame(vcf_dict)


    def check_out_of_bounds(self, svtype, chrom, chrom2, start, end, chrom_sizes, padding=50):
        """ Checks if SV is out of chromosome bounds. 

        param svtype: string, SV type
        param chrom: string, chromosome name
        param chrom2: string, chromosome name for translocations
        param start: int, SV start position
        param end: int, SV end position
        param chrom_sizes: pandas dataframe, chromosome sizes 
        param padding: int, padding to add to SV start and end positions for feature collection
        
        return: boolean, True if SV is out of bounds, False otherwise """

        if svtype != 'BND':
            return start - padding < 0 or end + padding > chrom_sizes.loc[chrom, 'size']
        else:
            # For translocations, check both chromosomes
            outbounds_chrom1 = start - padding < 0 or start + padding > chrom_sizes.loc[chrom, 'size']
            outbounds_chrom2 = end - padding < 0 or end + padding > chrom_sizes.loc[chrom2, 'size']
            return outbounds_chrom1 or outbounds_chrom2


    def read_vcf(self, filename):
        """ Reads VCF file as a pysam.VariantFile object and takes care of sample names. """

        if os.path.exists(filename):
            vcf = pysam.VariantFile(filename)
        else:
            # Smoove appears in lumpy filenames
            vcf = pysam.VariantFile(filename.replace('-', '_').replace('_smoove', '-smoove'))

        return vcf

    
    def read_variants(self):
        """ Reads all VCFs in params file and saves them in pandas dataframe. """

        vcf_dfs = []
        for tech in self.params['vcf']:
            for method in self.params['vcf'][tech]:
                filename_vcf = replace_filename(self.params['vcf'][tech][method], self.sample, self.ref)
                vcf = self.read_vcf(filename_vcf)
                df = self.parse_vcf(vcf, tech, method, self.sample)
                vcf_dfs.append(df)
        
        self.df_variants = pd.concat(vcf_dfs, ignore_index=True)


    def filter_variants(self):
        """ Removes variants that are out of chromosomes bounds or have other problems. """

        # Remove that are on non-canonical chromosomes
        self.df_variants = self.df_variants[(self.df_variants['type'] == 'BND') | (self.df_variants['chrom'].isin(self.chroms))].copy().reset_index(drop=True)
        self.df_variants.drop(self.df_variants[(self.df_variants['type'] == 'BND') & ((~self.df_variants['chrom'].isin(self.chroms)) | (~self.df_variants['chrom2'].isin(self.chroms)))].index, inplace=True)
        self.df_variants.reset_index(drop=True, inplace=True)

        # Remove calls that are out of chromosome bounds
        self.df_variants['start'] = self.df_variants['start'].astype(int)
        self.df_variants['end'] = self.df_variants['end'].astype(int)
        self.df_variants['outbounds'] = self.df_variants.apply(lambda x: self.check_out_of_bounds(x['type'], x['chrom'], x['chrom2'], x['start'], x['end'], self.chrom_sizes), axis=1)
        self.df_variants = self.df_variants[~self.df_variants['outbounds']].copy().drop('outbounds', axis=1).reset_index(drop=True)


    def save_variants(self):
        """ Saves variants dataframe to file. """

        self.df_variants.to_csv(self.out_file, index=False, sep='\t', na_rep='NA')
                


    



