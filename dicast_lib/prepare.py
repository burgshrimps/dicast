import pandas as pd
import pysam
import numpy as np
import re
import os

from dicast_lib.utils import replace_filename, caller_vcf_to_dataframe


class VariantPrep:
    """ Class to prepare raw variant calls for feature extraction. """

    def __init__(self, cohort: str, sample: str, ref: str, workdir: str, technology: str, vcfs: list, chroms: list, chrom_sizes: str, sv_types: list):
        """ Constructor for VariantPrep class.

        Args:
            cohort (str): Cohort name
            sample (str): Sample name
            ref (str): Reference genome name
            workdir (str): Working and output directory
            chroms (list): Chromosomes to use
            chrom_sizes (str): FAI file containing chromosome sizes
            sv_types (list): SV types supported by dicast
        """        

        # Input parameters
        self.cohort = cohort
        self.sample = sample
        self.ref = ref
        self.technology = technology
        self.workdir = workdir
        self.vcfs = vcfs
        self.sv_types = sv_types

        # Auxiallary files for preparation
        self.chrom_sizes = pd.read_csv(chrom_sizes, sep='\t', header=None, 
                                       names=['size', 'offset', 'linebases', 'linewidth'], index_col=0)

        # List of chromosomes to use
        self.chroms = chroms


    def read_variants(self):
        """ Reads all VCF files and stores them in pandas dataframe.
        """        

        vcf_dfs = []
        for caller, vcf_file in self.vcfs:
            
            # Open VCF file
            save = pysam.set_verbosity(0)
            if os.path.exists(vcf_file):
                vcf = pysam.VariantFile(vcf_file)
            else:
                # Smoove appears in lumpy filenames
                vcf = pysam.VariantFile(vcf_file.replace('-', '_'))
            pysam.set_verbosity(save)
            
            # Parse VCF file
            df = caller_vcf_to_dataframe(vcf, self.cohort, self.sample, self.ref, self.technology, caller, self.chroms)
            vcf_dfs.append(df)
            
        # Merge all VCF files
        self.df_variants = pd.concat(vcf_dfs, ignore_index=True)
            


    def check_out_of_bounds(self, svtype: str, chrom: str, chrom_2: str, start: int, end: int, chrom_sizes: pd.DataFrame, padding: int=50) -> bool:
        """ Checks if variant is out of chromosome bounds.

        Args:
            svtype (str): SV type
            chrom (str): Chromosome
            chrom2 (str): Second chromosome for translocations
            start (int): Start position
            end (int): End position
            chrom_sizes (pd.DataFrame): Dataframe with chromosome sizes
            padding (int, optional): Padding around SV borders. Defaults to 50.

        Returns:
            bool: True if variant is out of bounds, False otherwise
        """        

        if svtype != 'BND':
            return start - padding < 0 or end + padding > chrom_sizes.loc[chrom, 'size']
        else:
            # For translocations, check both chromosomes
            outbounds_chrom1 = start - padding < 0 or start + padding > chrom_sizes.loc[chrom, 'size']
            outbounds_chrom2 = end - padding < 0 or end + padding > chrom_sizes.loc[chrom_2, 'size']
            return outbounds_chrom1 or outbounds_chrom2


    def filter_variants(self):
        """ Removes variants that are out of chromosomes bounds or have other problems. """

        # Remove calls that are out of chromosome bounds
        self.df_variants['start'] = self.df_variants['start'].astype(int)
        self.df_variants['end'] = self.df_variants['end'].astype(int)
        self.df_variants['outbounds'] = self.df_variants.apply(lambda x: self.check_out_of_bounds(x['sv_type'], x['chrom'], x['chrom_2'], x['start'], x['end'], self.chrom_sizes), axis=1)
        self.df_variants = self.df_variants[~self.df_variants['outbounds']].copy().drop('outbounds', axis=1).reset_index(drop=True)
        
        # Remove SV types that are currently not supported by dicast
        self.df_variants = self.df_variants[self.df_variants['sv_type'].isin(self.sv_types)].copy().reset_index(drop=True)


    def save_variants(self):
        """ Saves variants dataframe to file. """
        
        filename = self.sample + '_' + self.ref + '.SVs.raw.tsv'
        self.df_variants.to_csv('/'.join([self.workdir, filename]), index=False, sep='\t', na_rep='NA')
                


    



