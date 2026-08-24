import pandas as pd
import pyBigWig
import bioframe as bf
import numpy as np
from tqdm import tqdm
import sys

from dicast_lib.utils import replace_filename

class ReferenceAnnotator:
    """ Object to annotate a set of SV calls based on features obtained from a reference genome. """


    def __init__(self, reference_filenames: dict):
        """ Constructor for ReferenceAnnotator class.

        Args:
            reference_filenames (dict): Dictionary containing filenames of reference files
        """        
        
        # Reference files
        self.repeats_filename = reference_filenames['repeats_filename']
        self.vntrs_filename = reference_filenames['vntrs_filename']
        self.strs_filename = reference_filenames['strs_filename']
        self.gc_filename = reference_filenames['gc_filename']
        self.cpgislands_filename = reference_filenames['cpgislands_filename']
        self.centromeres_filename = reference_filenames['centromeres_filename']
        self.asmb_gaps_filename = reference_filenames['asmb_gaps_filename']
        self.alt_haps_filename = reference_filenames['alt_haps_filename']
        
        # Load reference data
        self.df_repeats = pd.read_csv(self.repeats_filename, index_col=0, sep='\t')
        self.df_vntrs = pd.read_csv(self.vntrs_filename, sep='\t', header=None, names=['chrom', 'start', 'stop', 'class'])
        self.df_strs = pd.read_csv(self.strs_filename, sep='\t', header=None, names=['chrom', 'start', 'stop', 'len_unit', 'seq_unit', 'unknown'])
        self.bw_gc = pyBigWig.open(self.gc_filename)
        self.df_cpgislands = pd.read_csv(self.cpgislands_filename, sep='\t')
        self.df_centromeres = pd.read_csv(self.centromeres_filename, sep='\t')
        self.df_asmb_gaps = pd.read_csv(self.asmb_gaps_filename, sep='\t')
        self.df_alt_haps = pd.read_csv(self.alt_haps_filename, sep='\t')

            
    def load_from_csv(self, filename: str):
        """ Load SV calls from TSV file.

        Args:
            filename (str): Path to TSV file with raw SV calls
        """        
        
        self.df_calls = pd.read_csv(filename, sep='\t')
        self.df_calls_annot = self.df_calls[['id', 'cohort', 'sample', 'reference', 'technology', 'caller', 'sv_type', 'chrom', 'chrom_2', 'start', 'end']].copy()
        
        
    def load_from_df(self, df: pd.DataFrame):
        """ Load SV calls from pandas dataframe.

        Args:
            df (pd.DataFrame): Pandas dataframe with raw SV calls
        """        
        
        self.df_calls = df.copy()
        self.df_calls_annot = self.df_calls[['id', 'cohort', 'sample', 'reference', 'technology', 'caller', 'sv_type', 'chrom', 'chrom_2', 'start', 'end']].copy()
        
        
    def split_bnd(self):
        """ Split annotation dataframes into chromosomal and interchromosomal.
        """        
        
        if 'BND' in self.df_calls_annot['sv_type'].unique():
            self.split_bnd = True
            df_bnd_chrom1 = self.df_calls_annot.loc[self.df_calls_annot['sv_type'] == 'BND'].copy().reset_index(drop=True)
            df_bnd_chrom_2 = self.df_calls_annot.loc[self.df_calls_annot['sv_type'] == 'BND'].copy().reset_index(drop=True)
            self.df_calls_annot = self.df_calls_annot.loc[self.df_calls_annot['sv_type'] != 'BND'].copy().reset_index(drop=True)

            df_bnd_chrom1['end'] = df_bnd_chrom1['start'] + 50
            df_bnd_chrom1['start'] = df_bnd_chrom1['start'] - 50
            df_bnd_chrom1.drop('chrom_2', axis=1, inplace=True)
            self.df_calls_annot_bnd1 = df_bnd_chrom1

            df_bnd_chrom_2['chrom'] = df_bnd_chrom_2['chrom_2']
            df_bnd_chrom_2['start'] = df_bnd_chrom_2['end'] - 50
            df_bnd_chrom_2['end'] = df_bnd_chrom_2['end'] + 50
            df_bnd_chrom_2.drop('chrom_2', axis=1, inplace=True)
            self.df_calls_annot_bnd2 = df_bnd_chrom_2

        else:
            self.split_bnd = False


    def annotate_repeats(self):
        """ Annotate the variants with repeats called by repeatmasker. 
        """

        # Drop ambiguous RepeatMasker classes (trailing-? labels and "Unknown")
        self.df_repeats = self.df_repeats.loc[~self.df_repeats['repClass'].isin(['LTR?', 'Unknown', 'DNA?', 'RC?', 'SINE?'])].reset_index(drop=True)

        # Compute distance between variants and repeats, for translocation we need to do it twice; once for each chromosome
        df_ref = self.df_repeats[['genoName', 'genoStart', 'genoEnd', 'repClass']].rename(columns={'genoName': 'chrom', 'genoStart': 'start', 'genoEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        if not df_svs.empty:
            for rep_class in self.df_repeats['repClass'].unique():
                df_closest = bf.closest(df_svs, df_ref[df_ref['repClass'] == rep_class], return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_' + rep_class})
                self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'rep_' + rep_class]], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annoate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
            for rep_class in self.df_repeats['repClass'].unique():
                df_closest = bf.closest(df_svs, df_ref[df_ref['repClass'] == rep_class], return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_' + rep_class})
                self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'rep_' + rep_class]], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annoate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
            for rep_class in self.df_repeats['repClass'].unique():
                df_closest = bf.closest(df_svs, df_ref[df_ref['repClass'] == rep_class], return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_' + rep_class})
                self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'rep_' + rep_class]], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_vntrs(self):
        """ Annotate the variants with VNTRs from the Chaisson paper. 
        """

        # Compute distance between variants and VNTRs
        df_ref = self.df_vntrs.rename(columns={'stop': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        if not df_svs.empty:
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_VNTR'})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'rep_VNTR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annotate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_VNTR'})
            self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'rep_VNTR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annotate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_VNTR'})
            self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'rep_VNTR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_strs(self):
        """ Annotate the variants with STRs from the Chaisson paper. 
        """

        # Compute distance between variants and STRs
        df_ref = self.df_strs.rename(columns={'stop': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        if not df_svs.empty:
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_STR'})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'rep_STR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annotate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_STR'})
            self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'rep_STR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annotate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_STR'})
            self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'rep_STR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_cpg_islands(self):
        """ Annotate the variants with CpG islands from the UCSC Genome Browser Track. 
        """

        # Compute distance between variants and CpG islands
        df_ref = self.df_cpgislands[['chrom', 'chromStart', 'chromEnd']].rename(columns={'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        if not df_svs.empty:
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'cpg_islands'})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'cpg_islands']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annotate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'cpg_islands'})
            self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'cpg_islands']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annotate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'cpg_islands'})
            self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'cpg_islands']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_centromeres(self):
        """ Annotate the variants with distance to centromeres from the UCSC Genome Browser Track. 
        """

        # Compute distance between variants and centromeres
        df_ref = self.df_centromeres[['chrom', 'chromStart', 'chromEnd']].rename(columns={'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        if not df_svs.empty:
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'centromeres'})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'centromeres']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])
        
        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annotate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'centromeres'})
            self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'centromeres']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annotate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'centromeres'})
            self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'centromeres']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_asmb_gaps(self):
        """ Annotate the variants with distance to assembly gaps from the UCSC Genome Browser Track. 
        """

        # Compute distance between variants and assembly gaps
        df_ref = self.df_asmb_gaps[['chrom', 'chromStart', 'chromEnd']].rename(columns={'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        if not df_svs.empty:
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'asmb_gaps'})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'asmb_gaps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annotate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'asmb_gaps'})
            self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'asmb_gaps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annotate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'asmb_gaps'})
            self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'asmb_gaps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_alt_haps(self):
        """ Annotate the variants with distance to alternative haplotypes from the UCSC Genome Browser Track. 
        """

        # Compute distance between variants and alternative haplotypes
        df_ref = self.df_alt_haps[['tName', 'tStart', 'tEnd']].rename(columns={'tName': 'chrom', 'tStart' : 'start', 'tEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        if not df_svs.empty:
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'alt_haps'})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'alt_haps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            # Annotate BNDs on first chromosome
            df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'alt_haps'})
            self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'alt_haps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

            # Annotate BNDs on second chromosome
            df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
            df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'alt_haps'})
            self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'alt_haps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def calculate_gc_content(self, row: pd.Series, region: str) -> float:
        """ Calculate GC content of a region.

        Args:
            row (pd.Series): Row from dataframe df_calls_annot
            region (str): Either 'start' or 'end' of the SV

        Returns:
            float: GC content of the region
        """        

        try:
            return self.bw_gc.stats(row['chrom'], row[region] - 50, row[region] + 50)[0]
        except Exception:
            return np.nan

        
    def annotate_gc_content(self):
        """ Annotate the variants with GC content. 
        """
        
        if not self.df_calls_annot.empty:
            self.df_calls_annot['GC_content_left'] = self.df_calls_annot.apply(lambda x: self.calculate_gc_content(x, 'start'), axis=1)
            self.df_calls_annot['GC_content_right'] = self.df_calls_annot.apply(lambda x: self.calculate_gc_content(x, 'end'), axis=1)

        if self.split_bnd and not self.df_calls_annot_bnd1.empty:

            self.df_calls_annot_bnd1['GC_content_left'] = self.df_calls_annot_bnd1.apply(lambda x: self.calculate_gc_content(x, 'start'), axis=1)
            self.df_calls_annot_bnd1['GC_content_right'] = self.df_calls_annot_bnd1.apply(lambda x: self.calculate_gc_content(x, 'end'), axis=1)

            self.df_calls_annot_bnd2['GC_content_left'] = self.df_calls_annot_bnd2.apply(lambda x: self.calculate_gc_content(x, 'start'), axis=1)
            self.df_calls_annot_bnd2['GC_content_right'] = self.df_calls_annot_bnd2.apply(lambda x: self.calculate_gc_content(x, 'end'), axis=1)


    def aggregate_results(self):
        """ Save the dataframe to a csv file. 
        """

        if self.split_bnd:

            # Add BNDs to main frame
            self.df_calls_annot_bnd2.rename(columns={'chrom': 'chrom_2'}, inplace=True)
            self.df_calls_annot_bnd1['start'] = self.df_calls_annot_bnd1['start'] + 50
            self.df_calls_annot_bnd2['end'] = self.df_calls_annot_bnd2['end'] - 50

            self.df_calls_annot_bnd1['GC_content_left'] = (self.df_calls_annot_bnd1['GC_content_left'] + self.df_calls_annot_bnd1['GC_content_right']) / 2
            self.df_calls_annot_bnd1['GC_content_right'] = np.nan
            self.df_calls_annot_bnd2['GC_content_right'] = (self.df_calls_annot_bnd2['GC_content_left'] + self.df_calls_annot_bnd2['GC_content_right']) / 2
            self.df_calls_annot_bnd2['GC_content_left'] = np.nan

            df_calls_annot_bnd = pd.concat([self.df_calls_annot_bnd1, self.df_calls_annot_bnd2], ignore_index=True)
            df_calls_annot_bnd = df_calls_annot_bnd.groupby('id').agg({'cohort' : 'first', 'sample' : 'first', 'reference' : 'first', 'technology' : 'first', 'caller' : 'first', 'sv_type' : 'first', 'chrom' : 'first', 'chrom_2' : 'last', 'start' : 'first', 'end' : 'last',
                                                                    'rep_LINE' : 'min', 'rep_SINE' : 'min', 'rep_LTR' : 'min', 'rep_DNA' : 'min', 'rep_Simple_repeat' : 'min', 'rep_Satellite' : 'min', 'rep_Low_complexity' : 'min', 
                                                                    'rep_Retroposon' : 'min', 'rep_snRNA' : 'min', 'rep_tRNA' : 'min', 'rep_srpRNA' : 'min', 'rep_rRNA' : 'min', 'rep_RC' : 'min', 'rep_scRNA' : 'min',
                                                                    'rep_RNA' : 'min', 'rep_VNTR' : 'min', 'rep_STR' : 'min', 'cpg_islands' : 'min', 'centromeres' : 'min', 'asmb_gaps' : 'min', 'alt_haps' : 'min',
                                                                    'GC_content_left' : 'first', 'GC_content_right' : 'first'}).reset_index()
            if not self.df_calls_annot.empty:
                self.df_calls_annot = pd.concat([self.df_calls_annot, df_calls_annot_bnd[self.df_calls_annot.columns]], ignore_index=True)
            
            else:
                self.df_calls_annot = df_calls_annot_bnd.copy()

    def to_csv(self, filename: str):
        """ Save the dataframe to a csv file.

        Args:
            filename (str): Path to output TSV file
        """        
        
        self.df_calls_annot.to_csv(filename, index=False, na_rep='NA', sep='\t')
        
        
    def to_df(self):
        """ Return the dataframe. """
        
        return self.df_calls_annot