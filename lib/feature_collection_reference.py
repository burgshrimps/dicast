import pandas as pd
import pyBigWig
import bioframe as bf
import numpy as np
from tqdm import tqdm

from lib.utils import replace_filename


class ReferenceAnnotator:
    """ Object to annotate a set of SV calls based on features obtained from a reference genome. """


    def __init__(self, sample, ref, workdir, params):

        # Meta data
        self.sample = sample
        self.ref = ref
        self.workdir = workdir
        
        # Reference files
        ref_dir_root = params['ref'][self.ref]['directory']
        ref_dir_annot = params['ref'][self.ref]['subdirectory_annotation']
        self.ref_dir = '/'.join([ref_dir_root, ref_dir_annot])
        self.filename_repeats = params['ref'][self.ref]['filename_repeats']
        self.filename_vntrs = params['ref'][self.ref]['filename_vntrs']
        self.filename_strs = params['ref'][self.ref]['filename_strs']
        self.filename_gc = params['ref'][self.ref]['filename_gc']
        self.filename_cpgislands = params['ref'][self.ref]['filename_cpgislands']
        self.filename_centromeres = params['ref'][self.ref]['filename_centromeres']
        self.filename_asmb_gaps = params['ref'][self.ref]['filename_asmb_gaps']
        self.filename_alt_haps = params['ref'][self.ref]['filename_alt_haps']

        # Variant Files
        self.filename_variants = self.workdir + '/ensemble/' + self.sample + '_' + self.ref + '.SVs.raw.tsv'
        self.filename_variants_ref_annot = self.workdir + '/ensemble/' + self.sample + '_' + self.ref + '.SVs.ref.tsv'

        # Load data
        self.df_calls = pd.read_csv(self.filename_variants, sep='\t', low_memory=False)
        self.df_calls_annot = self.df_calls[['sample', 'id', 'type', 'chrom', 'chrom2', 'start', 'end']].copy()
        self.df_repeats = pd.read_csv('/'.join([self.ref_dir, self.filename_repeats]), index_col=0, sep='\t')
        self.df_vntrs = pd.read_csv('/'.join([self.ref_dir, self.filename_vntrs]), sep='\t', header=None, names=['chrom', 'start', 'stop', 'class'])
        self.df_strs = pd.read_csv('/'.join([self.ref_dir, self.filename_strs]), sep='\t', header=None, names=['chrom', 'start', 'stop', 'len_unit', 'seq_unit', 'unknown'])
        self.bw_gc = pyBigWig.open('/'.join([self.ref_dir, self.filename_gc]))
        self.df_cpgislands = pd.read_csv('/'.join([self.ref_dir, self.filename_cpgislands]), sep='\t')
        self.df_centromeres = pd.read_csv('/'.join([self.ref_dir, self.filename_centromeres]), sep='\t')
        self.df_asmb_gaps = pd.read_csv('/'.join([self.ref_dir, self.filename_asmb_gaps]), sep='\t')
        self.df_alt_haps = pd.read_csv('/'.join([self.ref_dir, self.filename_alt_haps]), sep='\t')

        # Split annotation dataframes into chromosomal and interchromosomal
        df_bnd_chrom1 = self.df_calls_annot.loc[self.df_calls_annot['type'] == 'BND'].copy().reset_index(drop=True)
        df_bnd_chrom2 = self.df_calls_annot.loc[self.df_calls_annot['type'] == 'BND'].copy().reset_index(drop=True)
        self.df_calls_annot = self.df_calls_annot.loc[self.df_calls_annot['type'] != 'BND'].copy().reset_index(drop=True)

        df_bnd_chrom1['end'] = df_bnd_chrom1['start'] + 50
        df_bnd_chrom1['start'] = df_bnd_chrom1['start'] - 50
        df_bnd_chrom1.drop('chrom2', axis=1, inplace=True)
        self.df_calls_annot_bnd1 = df_bnd_chrom1

        df_bnd_chrom2['chrom'] = df_bnd_chrom2['chrom2']
        df_bnd_chrom2['start'] = df_bnd_chrom2['end'] - 50
        df_bnd_chrom2['end'] = df_bnd_chrom2['end'] + 50
        df_bnd_chrom2.drop('chrom2', axis=1, inplace=True)
        self.df_calls_annot_bnd2 = df_bnd_chrom2


    def annotate_repeats(self):
        """ Annotate the variants with repeats. """

        # Remove weird repeat types
        self.df_repeats = self.df_repeats.loc[~self.df_repeats['repClass'].isin(['LTR?', 'Unknown', 'DNA?', 'RC?', 'SINE?'])].reset_index(drop=True)

        # Compute distance between variants and repeats, for translocation we need to do it twice; once for each chromosome
        df_ref = self.df_repeats[['genoName', 'genoStart', 'genoEnd', 'repClass']].rename(columns={'genoName': 'chrom', 'genoStart': 'start', 'genoEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        for rep_class in self.df_repeats['repClass'].unique():
            df_closest = bf.closest(df_svs, df_ref[df_ref['repClass'] == rep_class], return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_' + rep_class})
            self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'rep_' + rep_class]], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

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
        """ Annotate the variants with VNTRs from the Chaisson paper. """

        # Compute distance between variants and VNTRs
        df_ref = self.df_vntrs.rename(columns={'stop': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_VNTR'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'rep_VNTR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_VNTR'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'rep_VNTR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_VNTR'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'rep_VNTR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_strs(self):
        """ Annotate the variants with STRs from the Chaisson paper. """

        # Compute distance between variants and STRs
        df_ref = self.df_strs.rename(columns={'stop': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_STR'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'rep_STR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_STR'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'rep_STR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy().reset_index(drop=True)
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'rep_STR'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'rep_STR']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_cpg_islands(self):
        """ Annotate the variants with CpG islands from the UCSC Genome Browser Track. """

        # Compute distance between variants and CpG islands
        df_ref = self.df_cpgislands[['chrom', 'chromStart', 'chromEnd']].rename(columns={'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'cpg_islands'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'cpg_islands']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'cpg_islands'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'cpg_islands']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'cpg_islands'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'cpg_islands']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_centromeres(self):
        """ Annotate the variants with distance to centromeres from the UCSC Genome Browser Track. """

        # Compute distance between variants and centromeres
        df_ref = self.df_centromeres[['chrom', 'chromStart', 'chromEnd']].rename(columns={'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'centromeres'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'centromeres']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'centromeres'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'centromeres']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'centromeres'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'centromeres']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_asmb_gaps(self):
        """ Annotate the variants with distance to assembly gaps from the UCSC Genome Browser Track. """

        # Compute distance between variants and assembly gaps
        df_ref = self.df_asmb_gaps[['chrom', 'chromStart', 'chromEnd']].rename(columns={'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'asmb_gaps'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'asmb_gaps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'asmb_gaps'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'asmb_gaps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'asmb_gaps'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'asmb_gaps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_alt_haps(self):
        """ Annotate the variants with distance to alternative haplotypes from the UCSC Genome Browser Track. """

        # Compute distance between variants and alternative haplotypes
        df_ref = self.df_alt_haps[['tName', 'tStart', 'tEnd']].rename(columns={'tName': 'chrom', 'tStart' : 'start', 'tEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'alt_haps'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'alt_haps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'alt_haps'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'alt_haps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'alt_haps'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'alt_haps']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_genes(self):
        """ Annotate the variants with distance to genes (NCBI RefSeq) from the UCSC Genome Browser Track. """

        # Compute distance between variants and genes
        df_ref = self.df_genes[['chrom', 'txStart', 'txEnd']].rename(columns={'txStart' : 'start', 'txEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'genes'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'genes']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'genes'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'genes']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'genes'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'genes']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def annotate_orphanet(self):
        """ Annotate the variants with distance to orphanet diseases gene associations from the UCSC Genome Browser Track. """

        # Compute distance between variants and orphanet disease gene associations
        df_ref = self.df_orphanet[['#chrom', 'chromStart', 'chromEnd']].rename(columns={'#chrom': 'chrom', 'chromStart' : 'start', 'chromEnd': 'end'}).copy()

        # Annotate non BNDs
        df_svs = self.df_calls_annot[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'orphanet'})
        self.df_calls_annot = pd.merge(self.df_calls_annot, df_closest[['id_1', 'orphanet']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on first chromosome
        df_svs = self.df_calls_annot_bnd1[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'orphanet'})
        self.df_calls_annot_bnd1 = pd.merge(self.df_calls_annot_bnd1, df_closest[['id_1', 'orphanet']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])

        # Annotate BNDs on second chromosome
        df_svs = self.df_calls_annot_bnd2[['chrom', 'start', 'end', 'id']].copy()
        df_closest = bf.closest(df_svs, df_ref, return_overlap=True, suffixes=['_1', '_2']).rename(columns={'distance': 'orphanet'})
        self.df_calls_annot_bnd2 = pd.merge(self.df_calls_annot_bnd2, df_closest[['id_1', 'orphanet']], how='left', left_on='id', right_on='id_1').drop(columns=['id_1'])


    def calculate_gc_content(self, row, region):
        """ Calculate the GC content of the region around the variant. """

        try:
            return self.bw_gc.stats(row['chrom'], row[region] - 50, row[region] + 50)[0]
        except:
            return np.nan

        
    def annotate_gc_content(self):

        """ Annotate the variants with GC content. """
        self.df_calls_annot['GC_content_left'] = self.df_calls_annot.apply(lambda x: self.calculate_gc_content(x, 'start'), axis=1)
        self.df_calls_annot['GC_content_right'] = self.df_calls_annot.apply(lambda x: self.calculate_gc_content(x, 'end'), axis=1)

        self.df_calls_annot_bnd1['GC_content_left'] = self.df_calls_annot_bnd1.apply(lambda x: self.calculate_gc_content(x, 'start'), axis=1)
        self.df_calls_annot_bnd1['GC_content_right'] = self.df_calls_annot_bnd1.apply(lambda x: self.calculate_gc_content(x, 'end'), axis=1)

        self.df_calls_annot_bnd2['GC_content_left'] = self.df_calls_annot_bnd2.apply(lambda x: self.calculate_gc_content(x, 'start'), axis=1)
        self.df_calls_annot_bnd2['GC_content_right'] = self.df_calls_annot_bnd2.apply(lambda x: self.calculate_gc_content(x, 'end'), axis=1)


    def to_csv(self):
        """ Save the dataframe to a csv file. """

        # Add BNDs to main frame
        self.df_calls_annot_bnd2.rename(columns={'chrom': 'chrom2'}, inplace=True)
        self.df_calls_annot_bnd1['start'] = self.df_calls_annot_bnd1['start'] + 50
        self.df_calls_annot_bnd2['end'] = self.df_calls_annot_bnd2['end'] - 50

        self.df_calls_annot_bnd1['GC_content_left'] = (self.df_calls_annot_bnd1['GC_content_left'] + self.df_calls_annot_bnd1['GC_content_right']) / 2
        self.df_calls_annot_bnd1['GC_content_right'] = np.nan
        self.df_calls_annot_bnd2['GC_content_right'] = (self.df_calls_annot_bnd2['GC_content_left'] + self.df_calls_annot_bnd2['GC_content_right']) / 2
        self.df_calls_annot_bnd2['GC_content_left'] = np.nan

        df_calls_annot_bnd = pd.concat([self.df_calls_annot_bnd1, self.df_calls_annot_bnd2], ignore_index=True)
        df_calls_annot_bnd = df_calls_annot_bnd.groupby('id').agg({'sample': 'first', 'type' : 'first', 'chrom' : 'first', 'chrom2' : 'last', 'start' : 'first', 'end' : 'last',
                                                                   'rep_LINE' : min, 'rep_SINE' : min, 'rep_LTR' : min, 'rep_DNA' : min, 'rep_Simple_repeat' : min, 'rep_Satellite' : min, 'rep_Low_complexity' : min, 
                                                                   'rep_Retroposon' : min, 'rep_snRNA' : min, 'rep_tRNA' : min, 'rep_srpRNA' : min, 'rep_rRNA' : min, 'rep_RC' : min, 'rep_scRNA' : min,
                                                                   'rep_RNA' : min, 'rep_VNTR' : min, 'rep_STR' : min, 'cpg_islands' : min, 'centromeres' : min, 'asmb_gaps' : min, 'alt_haps' : min,
                                                                   'GC_content_left' : 'first', 'GC_content_right' : 'first'}).reset_index()

        self.df_calls_annot = pd.concat([self.df_calls_annot, df_calls_annot_bnd[self.df_calls_annot.columns]], ignore_index=True)
        self.df_calls_annot.to_csv(self.filename_variants_ref_annot, index=False, na_rep='NA', sep='\t')