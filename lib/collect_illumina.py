import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import logging

from lib.utils import replace_filename, mad

tqdm.pandas()


class AlignmentAnnotatorIllumina:
    """ Object to annotate a set of SV calls based on features obtained from an Illumina alignment file. """

    def __init__(self, sample, ref, workdir, params, chrom=None):
        """ Initialize class. """

        # Meta data
        self.sample = sample
        self.ref = ref
        self.workdir = workdir
        self.params = params
        self.chrom = chrom
        self.cov_thr = 6 # Threshold for log2 change in coverage to be considered for feature extraction, otherwise jump
        
        # Alignment file
        self.alignment_file = replace_filename(params['bam']['mgi'], sample, ref)
        if not os.path.exists(self.alignment_file):
            self.alignment_file = self.alignment_file.replace('-', '_')
        if not os.path.exists(self.alignment_file):
            raise FileNotFoundError(f'Alignment file {self.alignment_file} does not exist.')

        # Variant files
        self.filename_variants = self.workdir + '/ensemble/' + self.sample + '_' + self.ref + '.SVs.raw.tsv'
        if chrom != None:
            self.filename_variants_ill_annot = self.workdir + '/ensemble/' + self.sample + '_' + self.ref + '.SVs.aln.ill.' + self.chrom + '.tsv'
        else:
            self.filename_variants_ill_annot = self.workdir + '/ensemble/' + self.sample + '_' + self.ref + '.SVs.aln.ill.tsv'

        # Load data
        self.df_calls = pd.read_csv(self.filename_variants, sep='\t', low_memory=False)
        self.df_calls_annot = self.df_calls[['sample', 'id', 'type', 'chrom', 'chrom2', 'start', 'end']].copy()
        self.df_calls_annot = self.df_calls_annot.loc[self.df_calls_annot['chrom'] == chrom].reset_index(drop=True)
        self.bam = pysam.AlignmentFile(self.alignment_file, 'rb')

        # Split annotation dataframes into chromosomal and interchromosomal
        self.df_calls_annot_bnd = self.df_calls_annot.loc[self.df_calls_annot['type'] == 'BND'].copy().reset_index(drop=True)
        self.df_calls_annot = self.df_calls_annot.loc[self.df_calls_annot['type'] != 'BND'].copy().reset_index(drop=True)


    def calculate_coverage_baseline(self, s=100000, n=100):
        """ Calculates baseline mean and std coverage based on n sampled regions from the respective chromosome. """

        if self.chrom == 'chrX':
            chrom_idx = 22
        elif self.chrom == 'chrY':
            chrom_idx = 23
        else:
            chrom_idx = int(self.chrom[3:]) - 1
        self.baseline_coverage_mean = 0
        self.baseline_coverage_std = 0
        for i in range(n):
            start = np.random.randint(0, self.bam.lengths[chrom_idx])
            stop = start + s
            region = self.chrom + ':' + str(start) + '-' + str(stop)
            df = pd.DataFrame([x.split('\t') for x in pysam.depth(self.alignment_file, 
                                                                  '-r', region, '-a', '-g', 'SECONDARY,SUPPLEMENTARY').split('\n')[:-1]])
            df.rename({0 : 'chrom', 1 : 'pos', 2 : 'coverage'}, axis=1, inplace=True)
            df['coverage'] = df['coverage'].astype(int)
            self.baseline_coverage_mean += df['coverage'].mean()
            self.baseline_coverage_std += df['coverage'].std()

        self.baseline_coverage_mean /= n
        self.baseline_coverage_std /= n


    def calculate_insertsize_baseline(self, s=1000, n=1000):
        """ Calculates baseline mean and std insertsize based on n sampled regions from the respective chromosome. """

        if self.chrom == 'chrX':
            chrom_idx = 22
        elif self.chrom == 'chrY':
            chrom_idx = 23
        else:
            chrom_idx = int(self.chrom[3:]) - 1
        insert_sizes = []
        for i in range(n):
            start = np.random.randint(0, self.bam.lengths[chrom_idx] - s)
            stop = start + s
            for read in self.bam.fetch(self.chrom, start, stop):
                if not read.is_unmapped and not read.mate_is_unmapped:
                    insert_sizes.append(abs(read.template_length))
            
        self.baseline_insertsize_median = np.median(insert_sizes)
        self.baseline_insertsize_mad = mad(insert_sizes)


    def calculate_mapping_quality_baseline(self, s=1000, n=1000):
        """ Calculates baseline mean and std mapping quality based on n sampled regions from the respective chromosome. """

        if self.chrom == 'chrX':
            chrom_idx = 22
        elif self.chrom == 'chrY':
            chrom_idx = 23
        else:
            chrom_idx = int(self.chrom[3:]) - 1
        mapqs = []

        for i in range(n):
            start = np.random.randint(0, self.bam.lengths[chrom_idx] - s)
            stop = start + s
            for read in self.bam.fetch(self.chrom, start, stop):
                if not read.is_unmapped:
                    mapqs.append(read.mapping_quality)
            
        self.baseline_mapq_mean = np.mean(mapqs)
        self.baseline_mapq_std = np.std(mapqs)

    
    def calculate_coverage_region(self, chrom, start, stop, suffix):
        """ Calculates mean and std coverage of a region. """

        df = pd.DataFrame([x.split('\t') for x in pysam.depth(self.alignment_file, '-r', chrom + ':' + str(start) + '-' + str(stop), '-a', '-g', 'SECONDARY,SUPPLEMENTARY').split('\n')[:-1]])
        df.rename({0 : 'chrom', 1 : 'pos', 2 : 'coverage'}, axis=1, inplace=True)
        try:
            df['coverage'] = df['coverage'].astype(int)

             # add 0.1 to avoid -inf values
            coverage_mean = np.round(np.log2((df['coverage'].mean() + 0.1) / (self.baseline_coverage_mean + 0.1)), 3)
            coverage_std = np.round(np.log2((df['coverage'].std() + 0.1) / (self.baseline_coverage_std + 0.1)), 3)
        except KeyError:
            coverage_mean = np.nan
            coverage_std = np.nan

        return pd.Series([coverage_mean, coverage_std], index =['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix])


    def calculate_coverage(self, df, chrom_col, pos_col, offset_left, offset_right, suffix):
        """ Calculates mean and std coverage for all regions for a bin suffix. """
        df = df.copy()
        i = 0
        pbar = tqdm(total=len(df))
        all_exclude_idx = []
        while i < len(df):
            df.loc[i, ['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix]] = self.calculate_coverage_region(df.loc[i, chrom_col], df.loc[i, pos_col] - offset_left, df.loc[i, pos_col] + offset_right, suffix)

            # Mechanism to jump regions with extremly high coverage
            if df.loc[i, 'ill_cov_mean_' + suffix] > self.cov_thr:
                exclude_idx = df[(df[pos_col] >= df.loc[i, pos_col] - offset_left) & (df[pos_col] <= df.loc[i, pos_col] + offset_right)].index
                i += len(exclude_idx)
                pbar.update(len(exclude_idx))
                all_exclude_idx += list(exclude_idx)
            else:
                i += 1
                pbar.update(1)
        pbar.close()

        df.drop(all_exclude_idx, inplace=True)
        df.reset_index(drop=True, inplace=True)

        return df


    def annotate_coverage(self):
        """ Annotate the variants with coverage information. """
        print(self.chrom + ': Annotation Coverage')

        # Annotate non-BNDs
        self.df_calls_annot = self.calculate_coverage(self.df_calls_annot, 'chrom', 'start', 50, 0, 'I')
        self.df_calls_annot = self.calculate_coverage(self.df_calls_annot, 'chrom', 'start', 0, 50, 'II')
        self.df_calls_annot = self.calculate_coverage(self.df_calls_annot, 'chrom', 'end', 50, 0, 'III')
        self.df_calls_annot = self.calculate_coverage(self.df_calls_annot, 'chrom', 'end', 0, 50, 'IV')

        # Annotate BNDs
        self.df_calls_annot_bnd = self.calculate_coverage(self.df_calls_annot_bnd, 'chrom', 'start', 50, 0, 'I')
        self.df_calls_annot_bnd = self.calculate_coverage(self.df_calls_annot_bnd, 'chrom', 'start', 0, 50, 'II')
        self.df_calls_annot_bnd = self.calculate_coverage(self.df_calls_annot_bnd, 'chrom2', 'end', 50, 0, 'III')
        self.df_calls_annot_bnd = self.calculate_coverage(self.df_calls_annot_bnd, 'chrom2', 'end', 0, 50, 'IV')


    def get_overlap(self, a, b):
        """ Returns the overlap between two intervals. """
        return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1 )


    def get_clipped_span(self, read):
        """ Returns reference positions in which a read is clipped. """
        ref_start = read.reference_start
        ref_end = read.reference_end
        cigar = read.cigartuples
        if cigar[0][0] == 4 or cigar[0][0] == 5:
            # beginnging of read is clipped
            return (ref_start - cigar[0][1], ref_start - 1)
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            # end of read is clipped
            return (ref_end + 1, ref_end + cigar[-1][1])
        return None

        
    def calculate_read_based_features(self, chrom, start, stop, suffix):
        """ Collects features regarding insert size, mapping quality and split reads. """
        insert_sizes = []
        mapqs = []
        all_reads = 0
        all_reads_extended = 0
        clipped_reads = 0
        split_reads = 0
        disco_ff_reads = 0
        disco_rr_reads = 0

        for read in self.bam.fetch(chrom, start-5, stop+5):
            # add 5 bp padding because clipped bases cannot be used to fetch reads from a BAM file
            if not read.is_unmapped:
                if not read.reference_end <= start and not read.reference_start >= stop:
                    # only consider reads that overlap with the region for which we want to calculate the features
                    insert_sizes.append(abs(read.template_length))
                    mapqs.append(read.mapping_quality)
                    all_reads += 1
                    if read.has_tag('SA'):
                        split_reads += 1
                    if read.is_reverse and read.mate_is_reverse:
                        disco_rr_reads += 1
                    if not read.is_reverse and not read.mate_is_reverse:
                        disco_ff_reads += 1

                # for clipped reads, we need to extend the region by 5 bp
                all_reads_extended += 1
                clip_span = self.get_clipped_span(read)
                if clip_span != None:
                    overlap = self.get_overlap(clip_span, [start, stop])
                    if overlap > 0:
                        clipped_reads += 1
            
        if all_reads > 0:
             # add 0.1 to avoid -inf values
            insertsize_mean = np.round(np.log2((np.mean(insert_sizes) + 0.1) / (self.baseline_insertsize_median + 0.1)), 3)
            insertsize_std = np.round(np.log2((np.std(insert_sizes) + 0.1) / (self.baseline_insertsize_mad + 0.1)), 3)

            mapping_quality_mean = np.round(np.log2((np.mean(mapqs) + 0.1) / (self.baseline_mapq_mean + 0.1)), 3)
            mapping_quality_std = np.round(np.log2((np.std(mapqs) + 0.1) / (self.baseline_mapq_std + 0.1)), 3)

            splitreads_proportion = np.round(split_reads / all_reads, 3)
            clippedreads_proportion = np.round(clipped_reads / all_reads_extended, 3)
            disco_ff_proportion = np.round(disco_ff_reads / all_reads, 3)
            disco_rr_proportion = np.round(disco_rr_reads / all_reads, 3)
        else:
            # special case: no reads in region
            insertsize_mean = np.round(np.log2(0.1 / (self.baseline_insertsize_median + 0.1)), 3)
            insertsize_std = np.round(np.log2(0.1 / (self.baseline_insertsize_mad + 0.1)), 3)

            mapping_quality_mean = np.round(np.log2(0.1 / (self.baseline_mapq_mean + 0.1)), 3)
            mapping_quality_std = np.round(np.log2(0.1 / (self.baseline_mapq_std + 0.1)), 3)

            splitreads_proportion = 0
            clippedreads_proportion = 0
            disco_ff_proportion = 0
            disco_rr_proportion = 0

        return pd.Series([insertsize_mean, insertsize_std, mapping_quality_mean, mapping_quality_std, splitreads_proportion, clippedreads_proportion, disco_ff_proportion, disco_rr_proportion], index =['ill_isize_mean_' + suffix, 'ill_isize_std_' + suffix, 'ill_mapq_mean_' + suffix, 'ill_mapq_std_' + suffix, 'ill_splitreads_' + suffix, 'ill_clipreads_' + suffix, 'ill_disco_ff_' + suffix, 'ill_disco_rr_' + suffix])
        
        
    def annotate_read_based_features(self):
        """ Annotate the variants with read-based features. """
        print(self.chrom + ': Annotation Read-Based Features')
        
        # Annotate non-BNDs
        self.df_calls_annot.loc[:, ['ill_isize_mean_I', 'ill_isize_std_I', 'ill_mapq_mean_I', 'ill_mapq_std_I', 'ill_splitreads_I', 'ill_clipreads_I', 'ill_disco_ff_I', 'ill_disco_rr_I']] = self.df_calls_annot.progress_apply(lambda x: 
                                                                          self.calculate_read_based_features(x['chrom'], x['start'] - 50, 
                                                                          x['start'], 'I'), axis=1, result_type ='expand')
       
        self.df_calls_annot.loc[:, ['ill_isize_mean_II', 'ill_isize_std_II', 'ill_mapq_mean_II', 'ill_mapq_std_II', 'ill_splitreads_II', 'ill_clipreads_II', 'ill_disco_ff_II', 'ill_disco_rr_II']] = self.df_calls_annot.progress_apply(lambda x: 
                                                                            self.calculate_read_based_features(x['chrom'], x['start'], 
                                                                            x['start'] + 50, 'II'), axis=1, result_type ='expand')
       
        self.df_calls_annot.loc[:, ['ill_isize_mean_III', 'ill_isize_std_III', 'ill_mapq_mean_III', 'ill_mapq_std_III', 'ill_splitreads_III', 'ill_clipreads_III', 'ill_disco_ff_III', 'ill_disco_rr_III']] = self.df_calls_annot.progress_apply(lambda x: 
                                                                              self.calculate_read_based_features(x['chrom'], x['end'] - 50, 
                                                                              x['end'], 'III'), axis=1, result_type ='expand')

        self.df_calls_annot.loc[:, ['ill_isize_mean_IV', 'ill_isize_std_IV', 'ill_mapq_mean_IV', 'ill_mapq_std_IV', 'ill_splitreads_IV', 'ill_clipreads_IV', 'ill_disco_ff_IV', 'ill_disco_rr_IV']] = self.df_calls_annot.progress_apply(lambda x: 
                                                                              self.calculate_read_based_features(x['chrom'], x['end'], 
                                                                              x['end'] + 50, 'IV'), axis=1, result_type ='expand')

        # Annotate BNDs
        self.df_calls_annot_bnd.loc[:, ['ill_isize_mean_I', 'ill_isize_std_I', 'ill_mapq_mean_I', 'ill_mapq_std_I', 'ill_splitreads_I', 'ill_clipreads_I', 'ill_disco_ff_I', 'ill_disco_rr_I']] = self.df_calls_annot_bnd.progress_apply(lambda x: 
                                                                          self.calculate_read_based_features(x['chrom'], x['start'] - 50, 
                                                                          x['start'], 'I'), axis=1, result_type ='expand')
       
        self.df_calls_annot_bnd.loc[:, ['ill_isize_mean_II', 'ill_isize_std_II', 'ill_mapq_mean_II', 'ill_mapq_std_II', 'ill_splitreads_II', 'ill_clipreads_II', 'ill_disco_ff_II', 'ill_disco_rr_II']] = self.df_calls_annot_bnd.progress_apply(lambda x: 
                                                                            self.calculate_read_based_features(x['chrom'], x['start'], 
                                                                            x['start'] + 50, 'II'), axis=1, result_type ='expand')

        self.df_calls_annot_bnd.loc[:, ['ill_isize_mean_III', 'ill_isize_std_III', 'ill_mapq_mean_III', 'ill_mapq_std_III', 'ill_splitreads_III', 'ill_clipreads_III', 'ill_disco_ff_III', 'ill_disco_rr_III']] = self.df_calls_annot_bnd.progress_apply(lambda x: 
                                                                              self.calculate_read_based_features(x['chrom2'], x['end'] - 50, 
                                                                              x['end'], 'III'), axis=1, result_type ='expand')

        
        self.df_calls_annot_bnd.loc[:, ['ill_isize_mean_IV', 'ill_isize_std_IV', 'ill_mapq_mean_IV', 'ill_mapq_std_IV', 'ill_splitreads_IV', 'ill_clipreads_IV', 'ill_disco_ff_IV', 'ill_disco_rr_IV']] = self.df_calls_annot_bnd.progress_apply(lambda x: 
                                                                              self.calculate_read_based_features(x['chrom2'], x['end'], 
                                                                              x['end'] + 50, 'IV'), axis=1, result_type ='expand')


    def to_csv(self):
        """ Writes the annotated variants to a csv file. """

        print(self.chrom + ': Writing to CSV')

        self.df_calls_annot = pd.concat([self.df_calls_annot, self.df_calls_annot_bnd], ignore_index=True)

        alignment_columns = []
        for feature in ['ill_cov_mean_', 'ill_cov_std_', 'ill_isize_mean_', 'ill_isize_std_', 'ill_mapq_mean_', 'ill_mapq_std_', 
                        'ill_clipreads_', 'ill_splitreads_', 'ill_disco_ff_', 'ill_disco_rr_']:
            for suffix in ['I', 'II', 'III', 'IV']:
                alignment_columns.append(feature + suffix)
        columns_reordered = ['sample', 'id', 'type', 'chrom', 'chrom2', 'start', 'end'] + alignment_columns
        self.df_calls_annot = self.df_calls_annot[columns_reordered]
        self.df_calls_annot.to_csv(self.filename_variants_ill_annot, index=False, na_rep='NA', sep='\t')
        self.bam.close()