import pysam
import pandas as pd
import numpy as np
from tqdm import tqdm
import re
import time

from dicast_lib.utils import mad


class AlignmentAnnotatorIllumina:
    """ Class for annotating SVs with features derived from the alignment of short reads. """    

    def __init__(self, bam_filename: str, chrom: str, sv_type: str, sample: str, log_file: str = None, job_id: str = None):
        """ Initializes the class.

        Args:
            bam_filename (str): Path to the BAM file containing the alignment of short reads
            chrom (str): Chromosome name
            sv_type (str): Structural variant type (DEL, DUP, INV, INS, BND)
            log_file (str): Path to log file for writing
            job_id (str): Job identifier for logging
        """        

        # Meta data
        self.sample = str(sample)
        self.chrom = chrom
        self.sv_type = sv_type
        self.log_file = log_file
        self.job_id = job_id or f'{sample}_{chrom}_{sv_type}'
        self.features_breakpoints = ['ill_cov_mean_', 'ill_cov_std_', 'ill_isize_mean_', 'ill_isize_std_', 'ill_mapq_mean_', 'ill_mapq_std_', 
                                     'ill_clipreads_', 'ill_splitreads_', 'ill_disco_ff_', 'ill_disco_rr_', 'ill_disco_rf_', 'ill_disco_tx_']
        self.features_body = ['ill_cov_mean_', 'ill_cov_std_']
        self.features_connection = ['ill_disco_ff_', 'ill_disco_rr_', 'ill_disco_rf_', 'ill_splitreads_']
        self.cov_thr = 3 # Threshold for log2 change in coverage to be considered for feature extraction, otherwise jump
        self.max_fetch_reads = 10000  # Hard cap on reads fetched per bin
        
        # Alignment file'
        self.alignment_file = bam_filename
        self.bam = pysam.AlignmentFile(self.alignment_file, 'rb')
        
        # Log message for tqdm
        self.log_message = ' '.join([self.sample, self.chrom, self.sv_type])
        
        # Define bins around SV to calculate features for
        # Values are: [chrom, pos_col_left, pos_col_right, offset_left, offset_right, check_len]
        if self.sv_type == 'DEL' or self.sv_type == 'DUP' or self.sv_type == 'INV':
            # Bins around the breakpoints
            self.bin_dict_bps = {'I' : ['chrom', 'start', 'start', -52, -2],
                                 'II' : ['chrom', 'start', 'start', 2, 52],
                                 'III' : ['chrom', 'end', 'end', -52, -2],
                                 'IV' : ['chrom', 'end', 'end', 2, 52]}
            
            # Bins inside the SV body
            self.bin_dict_body = {'IIa' : ['chrom', 'start', 'body_I', 52, 0],
                                  'IIb' : ['chrom', 'body_I', 'body_II', 0, 0],
                                  'IIIb' : ['chrom', 'body_II', 'body_III', 0, 0],
                                  'IIIa' : ['chrom', 'body_III', 'end', 0, -52]}
            
        elif self.sv_type == 'INS':
            # Bins around the breakpoints
            self.bin_dict_bps = {'I' : ['chrom', 'start', 'start', -52, -2],
                                 'II' : ['chrom', 'start', 'start', 2, 52]}
            
        elif self.sv_type == 'BND':
            # Bins around the breakpoints
            self.bin_dict_bps = {'I' : ['chrom', 'start', 'start', -52, -2],
                                'II' : ['chrom', 'start', 'start', 2, 52],
                                'III' : ['chrom_2', 'end', 'end', -52, -2],
                                'IV' : ['chrom_2', 'end', 'end', 2, 52]}
        
        # Define connections between bins    
        self.bin_connections = {'I' : ['II', 'III', 'IV'],
                                'II' : ['III', 'IV'],
                                'III' : ['IV'],
                                'IV' : []}
        
        
    def _log(self, msg: str, level: str = 'INFO'):
        """Write log message to file and stdout."""
        timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        log_msg = f'{timestamp} [{self.job_id}] {level}: {msg}'
        if self.log_file:
            with open(self.log_file, 'a') as f:
                f.write(log_msg + '\n')
                f.flush()
        print(log_msg, flush=True)

    def prepare_dataframe(self):
        """ Prepares DataFrame for SV annotation. """
        
        self.df_calls_annot = self.df_calls[['id', 'cohort', 'sample', 'reference', 'technology', 'caller', 'sv_type', 'chrom', 'chrom_2', 'start', 'end', 'sv_len']].copy()
        
        # Check if correct subset of SVs is being processed
        self.df_calls_annot = self.df_calls_annot.loc[self.df_calls_annot['chrom'] == self.chrom].reset_index(drop=True)
        self.df_calls_annot = self.df_calls_annot.loc[self.df_calls_annot['sv_type'] == self.sv_type].reset_index(drop=True)
        
        # Add columns such that DataFrame has all the necessary columns
        for feature in self.features_breakpoints:
            for suffix in ['I', 'II', 'III', 'IV']:
                self.df_calls_annot[feature + suffix] = np.nan
                
        for feature in self.features_body:
            for suffix in ['IIa', 'IIb', 'IIIb', 'IIIa']:
                self.df_calls_annot[feature + suffix] = np.nan
                
        for feature in self.features_connection:
            suffices = ['I', 'II', 'III', 'IV']
            for i in range(4):
                for j in range(i+1, 4):
                    self.df_calls_annot[feature + suffices[i] + '_' + suffices[j]] = np.nan
    
    
    def load_from_csv(self, variants_filename: str):
        """ Loads SVs from file.

        Args:
            variants_filename (str): Path to the file containing the SVs
        """        
        
        self.df_calls = pd.read_csv(variants_filename, sep='\t', low_memory=False)
        self.prepare_dataframe()
        
    
    def load_from_df(self, df: pd.DataFrame):
        """ Loads SVs from a DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing the SVs
        """        
        
        self.df_calls = df.copy()
        self.prepare_dataframe()
        
        
    def calculate_coverage_baseline(self, s: int=1000, n: int=1000):
        """ Calculates baseline mean and std coverage based on n sampled regions of size s from the respective chromosome.

        Args:
            s (int, optional): Size of sampled region. Defaults to 100000.
            n (int, optional): Number of sampled regions. Defaults to 100.
        """        

        if self.chrom == 'chrX':
            chrom_idx = 22
        elif self.chrom == 'chrY':
            chrom_idx = 23
        elif self.chrom == 'chrM':
            chrom_idx = 24
        else:
            chrom_idx = int(self.chrom[3:]) - 1
            
        self.baseline_coverage_mean = 0
        self.baseline_coverage_std = 0
        baseline_coverage_mean = []
        baseline_coverage_std = []
        
        for i in range(n):
            if len(self.bam.lengths)>1:
                start = np.random.randint(0, self.bam.lengths[chrom_idx])
            else:
                start = np.random.randint(0, self.bam.lengths[0])            
            stop = start + s
            region = self.chrom + ':' + str(start) + '-' + str(stop)
            
            coverage = [0] * (stop - start + 1)
            for pileupcolumn in self.bam.pileup(self.chrom, start, stop, min_mapping_quality=20, flag_filter=1540, stepper='samtools', ignore_orphans=False, ignore_overlaps=False):
                if start <= pileupcolumn.pos <= stop:
                    coverage[pileupcolumn.pos - start] = pileupcolumn.nsegments
            
            baseline_coverage_mean.append(np.mean(coverage))
            baseline_coverage_std.append(np.std(coverage))

        self.baseline_coverage_mean = np.median(baseline_coverage_mean)
        self.baseline_coverage_std = np.median(baseline_coverage_std)
        
        
    def calculate_insertsize_baseline(self, s: int=1000, n: int=1000):
        """ Calculates baseline median and MAD insert size based on n sampled regions of size s from the respective chromosome.

        Args:
            s (int, optional): Size of sampled region. Defaults to 1000.
            n (int, optional): Number of sampled regions. Defaults to 1000.
        """        

        if self.chrom == 'chrX':
            chrom_idx = 22
        elif self.chrom == 'chrY':
            chrom_idx = 23
        elif self.chrom == 'chrM':
            chrom_idx = 24
        else:
            chrom_idx = int(self.chrom[3:]) - 1
        insert_sizes = []
        for i in range(n):
            if len(self.bam.lengths)>1:
                start = np.random.randint(0, self.bam.lengths[chrom_idx])
            else:
                start = np.random.randint(0, self.bam.lengths[0]) 
            stop = start + s
            for read in self.bam.fetch(self.chrom, start, stop):
                if not read.is_unmapped and not read.mate_is_unmapped:
                    insert_sizes.append(abs(read.template_length))
            
        self.baseline_insertsize_median = np.median(insert_sizes)
        self.baseline_insertsize_mad = mad(insert_sizes)
        
        
    def calculate_mapping_quality_baseline(self, s=1000, n=1000):
        """ Calculates baseline mean and std mapping quality based on n sampled regions of size s from the respective chromosome.

        Args:
            s (int, optional): Size of sampled region. Defaults to 1000.
            n (int, optional): Number of sampled regions. Defaults to 1000.
        """      

        if self.chrom == 'chrX':
            chrom_idx = 22
        elif self.chrom == 'chrY':
            chrom_idx = 23
        elif self.chrom == 'chrM':
            chrom_idx = 24
        else:
            chrom_idx = int(self.chrom[3:]) - 1
        mapqs = []

        for i in range(n):
            if len(self.bam.lengths)>1:
                start = np.random.randint(0, self.bam.lengths[chrom_idx])
            else:
                start = np.random.randint(0, self.bam.lengths[0]) 
            stop = start + s
            for read in self.bam.fetch(self.chrom, start, stop):
                if not read.is_unmapped:
                    mapqs.append(read.mapping_quality)
            
        self.baseline_mapq_mean = np.mean(mapqs)
        self.baseline_mapq_std = np.std(mapqs)
        
        
    def jump_to_next_variant_for_coverage_calculation(self, df: pd.DataFrame, i: int, pos_col_left: str, pos_col_right: str, offset_left: int, offset_right: int, suffix: str) -> pd.DataFrame:
        """ Checks if the mean coverage in a region exceeds the threshold. If yes, all variants overlapping that bin are excluded from further analysis.

        Args:
            df (pd.DataFrame): DataFrame containing the SVs
            i (int): Index of current SV
            pos_col_left (str): Column name for the left position
            pos_col_right (str): Column name for the right position
            offset_left (int):  Position +- Offset Left = Start position of the bin
            offset_right (int): Position +- Offset Right = End position of the bin
            suffix (str):  Suffix for column names

        Returns:
            pd.DataFrame: _description_
        """        
        
        if df.loc[i, 'ill_cov_mean_' + suffix] > self.cov_thr:
            mask_left = (df[(df[pos_col_left] + offset_left >= df.loc[i, pos_col_left] + offset_left) & (df[pos_col_left] + offset_left <= df.loc[i, pos_col_right] + offset_right)])
            mask_right = (df[(df[pos_col_right] + offset_right >= df.loc[i, pos_col_left] + offset_left) & (df[pos_col_right] + offset_right <= df.loc[i, pos_col_right] + offset_right)])
            exclude_idx = list(set(mask_left.index) | set(mask_right.index))
            
            return len(exclude_idx), list(exclude_idx)
        
        else:
            return 1, []         
        
        
    def calculate_coverage_region(self, chrom: str, start: int, stop: int, suffix: str) -> pd.Series:
        """ Calculates mean and std coverage in a region.

        Args:
            chrom (str): Chromosome name
            start (int): Start position
            stop (int): End position
            suffix (str): Suffix for column names

        Returns:
            pd.Series: Series containing mean and std coverage
        """        

        # Fetch + cumulative sum replacement for pileup
        # We basically retrieve all reads overlapping the region
        # and clip their start and end position to the regions exact interval
        # Then count read "starts" and "ends" and use the cumulative sum to get the coverage
        self._log(f'Calculating coverage for type {self.sv_type}_{suffix}:{chrom}:{start}-{stop}')
        n = stop - start + 1
        rs_list: list = []
        re_list: list = []
        for _read in self.bam.fetch(chrom, start, stop):
            if not (_read.flag & 1540) and _read.mapping_quality >= 20:
                rs_list.append(_read.reference_start)
                re_list.append(_read.reference_end)
        if rs_list:
            _rs = np.clip(np.array(rs_list, dtype=np.int32) - start, 0, n)
            _re = np.clip(np.array(re_list, dtype=np.int32) - start, 0, n)
            _delta = np.bincount(_rs, minlength=n + 1).astype(np.int32)
            _delta -= np.bincount(_re, minlength=n + 1).astype(np.int32)
            coverage = np.cumsum(_delta)[:n].tolist()
        else:
            coverage = [0] * n
                
        coverage_mean = np.round(np.log2((np.mean(coverage) + 0.1) / (self.baseline_coverage_mean + 0.1)), 3)
        coverage_std = np.round(np.log2((np.std(coverage) + 0.1) / (self.baseline_coverage_std + 0.1)), 3)
        
        return pd.Series([coverage_mean, coverage_std], index =['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix])
    
    
    def calculate_coverage(self, df: pd.DataFrame, chrom_col: str, pos_col_left: str, pos_col_right: str, offset_left: int, offset_right: int, suffix: str, check_len: bool = False) -> pd.DataFrame:
        """ Calculates mean and std coverage for one bin for every SV in the DataFrame.

        Args:
            df (pd.DataFrame): DataFrame containing the SVs
            chrom_col (str): Column name for the chromosome
            pos_col (str): Column name for the position
            offset_left (int): Position - Offset Left = Start position of the bin
            offset_right (int): Position + Offset Right = End position of the bin
            suffix (str): Suffix for column names

        Returns:
            pd.DataFrame: DataFrame containing the SVs with the added columns for mean and std coverage
        """        
        
        df = df.copy()
        i = 0
        all_exclude_idx = []
        #pbar = tqdm(total=len(df), desc=' '.join(['Annotation Coverage:', self.log_message, suffix]), position=0, leave=True)
        coverage_dict = {}
        
        while i < len(df):
            
            if check_len:
                sv_len = df.loc[i, 'sv_len']
                
                if sv_len < 150:
                    breakpoint_suffix = suffix[:-1]
                    df.loc[i, ['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix]] = df.loc[i, ['ill_cov_mean_' + breakpoint_suffix, 'ill_cov_std_' + breakpoint_suffix]]
                    i += 1
                    #pbar.update(1)
                    continue
                     
                elif sv_len < 5000:
                    left_border = df.loc[i, pos_col_left] + offset_left
                    right_border = df.loc[i, pos_col_right] + offset_right
                
                else:
                    left_border = df.loc[i, pos_col_left] + offset_left
                    right_border = left_border + 50
                    
                df.loc[i, ['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix]] = self.calculate_coverage_region(df.loc[i, chrom_col], left_border, right_border, suffix)
                i += 1
                #pbar.update(1)
                
            else:       
                left_border = df.loc[i, pos_col_left] + offset_left
                right_border = df.loc[i, pos_col_right] + offset_right

                # the use of a coverage dictionary is mainly for BNDs that share the exact same breakpoint hotspots
                if f'{df.loc[i, chrom_col]}_{left_border}_{right_border}' in coverage_dict:
                    df.loc[i, ['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix]] = coverage_dict[f'{df.loc[i, chrom_col]}_{left_border}_{right_border}']
                else:
                    df.loc[i, ['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix]] = self.calculate_coverage_region(df.loc[i, chrom_col], left_border, right_border, suffix)
                    coverage_dict[f'{df.loc[i, chrom_col]}_{left_border}_{right_border}'] = df.loc[i, ['ill_cov_mean_' + suffix, 'ill_cov_std_' + suffix]]
                step_size, exclude_idx = self.jump_to_next_variant_for_coverage_calculation(df, i, pos_col_left, pos_col_right, offset_left, offset_right, suffix)
                i += step_size
                all_exclude_idx += exclude_idx
                #pbar.update(step_size)

        df.drop(all_exclude_idx, inplace=True)
        df.reset_index(drop=True, inplace=True)

        return df
    
    
    def divide_sv_body(self, start: int, end: int) -> pd.Series:
        
        bin_edges = np.linspace(start, end, 5).astype(int)
        
        return pd.Series(bin_edges[1:-1], index=['body_I', 'body_II', 'body_III'])
            
        
    def annotate_coverage(self):
        """ Annotates SVs with mean and std coverage for all bins around the SV breakpoints. """        
        
        self._log(f'Starting coverage annotation: {len(self.df_calls_annot)} variants')
        start_coverage_time = time.time()

        # Calculate coverage around breakpoints
        for suffix, values in self.bin_dict_bps.items():
            variants_before = len(self.df_calls_annot)
            self.df_calls_annot = self.calculate_coverage(self.df_calls_annot, values[0], values[1], values[2], values[3], values[4], suffix)
            variants_after = len(self.df_calls_annot)
            excluded = variants_before - variants_after
            elapsed = time.time() - start_coverage_time
            #self._log(f'  Bin {suffix}: {variants_before:>4d} → {variants_after:>4d} vars ({excluded:>3d} excl), elapsed {elapsed:>6.1f}s', level='BIN')
            
        if self.sv_type == 'DEL' or self.sv_type == 'DUP' or self.sv_type == 'INV':    
            
            # Generate borders for bins inside the SV body
            self.df_calls_annot.loc[:, ['body_I', 'body_II', 'body_III']] = self.df_calls_annot.apply(lambda x: self.divide_sv_body(x['start'] + 52, x['end'] - 52), axis=1, result_type ='expand')
            
            # Calculate coverage inside the SV body
            for suffix, values in self.bin_dict_body.items():
                variants_before = len(self.df_calls_annot)
                self.df_calls_annot = self.calculate_coverage(self.df_calls_annot, values[0], values[1], values[2], values[3], values[4], suffix, check_len=True)
                variants_after = len(self.df_calls_annot)
                excluded = variants_before - variants_after
                elapsed = time.time() - start_coverage_time
                #self._log(f'  Body {suffix}: {variants_before:>4d} → {variants_after:>4d} vars ({excluded:>3d} excl), elapsed {elapsed:>6.1f}s', level='BIN')
                
            self.df_calls_annot.drop(['body_I', 'body_II', 'body_III'], axis=1, inplace=True)

        total_elapsed = time.time() - start_coverage_time
        self._log(f'DONE: {len(self.df_calls_annot)} variants remaining after coverage filtering (total {total_elapsed:.1f}s)', level='DONE')
            
            
    def get_overlap(self, a: tuple, b: tuple) -> int:
        """ Calculates the overlap between two intervals.

        Args:
            a (tuple): Interval a = (start, end)
            b (tuple): Interval b = (start, end)

        Returns:
            int: Overlap between a and b
        """        
        
        return max(0, min(a[1], b[1]) - max(a[0], b[0]) + 1)
    
    
    def get_clipped_span(self, read: pysam.AlignedSegment) -> tuple:
        """ Calculates the span of a read that is clipped.

        Args:
            read (pysam.AlignedSegment): Sequencing Read

        Returns:
            tuple: Reference start and end of the clipped part of the read
        """        
        
        ref_start = read.reference_start
        ref_end = read.reference_end
        cigar = read.cigartuples
        
        # Beginnging of read is clipped
        if cigar[0][0] == 4 or cigar[0][0] == 5:
            return (ref_start - cigar[0][1], ref_start - 1)
        
        # End of read is clipped
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            return (ref_end + 1, ref_end + cigar[-1][1])
        
        return None
    
    
    def get_end_position(self, start: int, cigar: str) -> int:
        """ Calculates the end position of a read based on its start position and CIGAR string.

        Args:
            start (int): Start position
            cigar (str): CIGAR string

        Returns:
            int: End position
        """        
    
        cigar_tuples = [(int(length), op) for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
        
        # Because positions are inclusive
        end = start - 1
        for length, op in cigar_tuples:
            if op in ["M", "D", "N", "=", "X"]:
                end += length

        return end
    
    
    def get_overlap_mate_bins(self, read: pysam.AlignedSegment, bins: list) -> np.array:
        """ Calculates the overlap between the mate of a read and all bins.

        Args:
            read (pysam.AlignedSegment): Sequencing Read
            bins (list): List of tuples containing the start and end position of the bins

        Returns:
            np.array: For each bin, 1 if there is an overlap, 0 otherwise
        """        
    
        conns = np.zeros(len(bins))
        
        if read.is_mapped:
            for i, bin in enumerate(bins):
                conns[i] += int(bool(self.get_overlap((read.next_reference_start, read.next_reference_start + 150), bin)))
        
        return conns
    
    
    def suffix_to_bin_idx(self, suffix: str) -> int:
        """ Converts a suffix to a bin index.

        Args:
            suffix (str): Suffix for column names

        Returns:
            int: Bin index
        """        
        
        suffix_dict = {'I' : 1, 'II' : 2, 'III' : 3, 'IV' : 4}
        
        return suffix_dict[suffix]
    
    
    def calculate_read_based_features(self, chrom: str, start: int, end: int, sv_start: int, sv_end: int, suffix: str) -> pd.Series:
        """ Calculates read-based features for a region.

        Args:
            bam (pysam.AlignmentFile): BAM file
            chrom (str): Chromosome name
            bin (int): bin number
            bin_start (int): Bin Start position
            bin_stop (int): Bin end position
            start (int): SV start position
            end (int): SV end position
            suffix (str): Suffix for column names

        Returns:
            pd.Series: Series containing the read-based features
        """  
        
        # Set bins for overlap computations
        bin_idx = self.suffix_to_bin_idx(suffix)
        bins = [(sv_start + 2, sv_start + 52), (sv_end - 52, sv_end - 2), (sv_end + 2, sv_end + 52)][bin_idx-1:]      

        # Initialize features
        insert_sizes = []
        mapqs = []
        all_reads = 0
        all_reads_extended_ids = set()
        clipped_reads = 0
        split_reads = 0
        split_reads_conn = np.zeros(4-bin_idx)
        disco_ff_reads = 0
        disco_ff_conn = np.zeros(4-bin_idx)
        disco_rr_reads = 0
        disco_rr_conn = np.zeros(4-bin_idx)
        disco_rf_reads = 0
        disco_rf_conn = np.zeros(4-bin_idx)
        disco_tx_reads = 0

        # Define labels for binned split read features
        labels_split_reads_conn = [f'ill_splitreads_{suffix}_II', f'ill_splitreads_{suffix}_III', f'ill_splitreads_{suffix}_IV'][bin_idx-1:]
        
        # Define labels for binned discordant read pair features
        # CHANGE THIS TO THREE ARRAYS
        labels_disco_conn = np.array([[f'ill_disco_ff_{suffix}_II', f'ill_disco_ff_{suffix}_III', f'ill_disco_ff_{suffix}_IV'], 
                                      [f'ill_disco_rr_{suffix}_II', f'ill_disco_rr_{suffix}_III', f'ill_disco_rr_{suffix}_IV'], 
                                      [f'ill_disco_rf_{suffix}_II',f'ill_disco_rf_{suffix}_III', f'ill_disco_rf_{suffix}_IV']])[:,bin_idx-1:]
        labels_disco_conn = list(labels_disco_conn.flatten())
        
        # Iterate over reads in region
        reads_fetched = 0
        for read in self.bam.fetch(chrom, start - 20, end + 20):
            
            reads_fetched += 1
            # for regions exceeding the max fetch reads, reset all counts and break
            # these regions are highly repetitive and cannot be properly scored anyways
            if reads_fetched > self.max_fetch_reads:
                all_reads = 0
                all_reads_extended_ids = set()
                split_reads_conn[:] = 0
                disco_ff_conn[:] = 0
                disco_rr_conn[:] = 0
                disco_rf_conn[:] = 0
                break

            if not read.is_unmapped and not read.is_duplicate and not read.is_qcfail:
                
                # Only consider reads that overlap with the region for which we want to calculate the features
                if not read.reference_end <= start and not read.reference_start >= end:
                    
                    # Insert size
                    insert_sizes.append(abs(read.template_length))
                    
                    # Mapping quality
                    mapqs.append(read.mapping_quality)
                        
                    # Read orientation for inversions
                    if read.is_reverse and read.mate_is_reverse:
                        disco_rr_reads += 1
                        disco_rr_conn += self.get_overlap_mate_bins(read, bins)
                        
                    elif not read.is_reverse and not read.mate_is_reverse:
                        disco_ff_reads += 1
                        disco_ff_conn += self.get_overlap_mate_bins(read, bins)
                    
                    # Read orientation for duplications
                    elif read.is_read1 and read.is_reverse and not read.mate_is_reverse and read.template_length>0:
                        disco_rf_reads += 1
                        disco_rf_conn += self.get_overlap_mate_bins(read, bins)
                        
                    elif read.is_read1 and not read.is_reverse and read.mate_is_reverse and read.template_length<0:
                        disco_rf_reads += 1
                        disco_rf_conn += self.get_overlap_mate_bins(read, bins)
                        
                    elif read.is_read2 and not read.is_reverse and read.mate_is_reverse and read.template_length<0:
                        disco_rf_reads += 1
                        disco_rf_conn += self.get_overlap_mate_bins(read, bins)
                        
                    elif read.is_read2 and read.is_reverse and not read.mate_is_reverse and read.template_length>0:
                        disco_rf_reads += 1
                        disco_rf_conn += self.get_overlap_mate_bins(read, bins)  
                        
                    # Read orientation for translocations
                    elif read.reference_id != read.next_reference_id:
                        disco_tx_reads += 1
                    
                    # Count number of reads for normalization 
                    read_id = read.query_name + '_' + str(int(read.is_read1)) + str(int(read.is_supplementary)) + str(int(read.is_secondary))
                    all_reads_extended_ids.add(read_id)
                    all_reads += 1

                # Clipped reads
                clip_span = self.get_clipped_span(read)
                
                if clip_span != None:
                    overlap = self.get_overlap(clip_span, [start, end])
                    
                    if overlap > 0:
                        clipped_reads += 1
                        read_id = read.query_name + '_' + str(int(read.is_read1)) + str(int(read.is_supplementary)) + str(int(read.is_secondary))
                        all_reads_extended_ids.add(read_id)
                        
                        # Split-reads
                        if read.has_tag('SA'):
                            split_reads += 1
                            
                            # Split-read connections
                            supplementary_alignment = read.get_tag('SA').split(',')
                            sa_cigar = supplementary_alignment[3]
                            sa_chrom = supplementary_alignment[0]
                            sa_start = int(supplementary_alignment[1])
                            sa_end = self.get_end_position(sa_start, sa_cigar)
                            
                            if sa_chrom == chrom:
                                for i, bin in enumerate(bins):
                                    split_reads_conn[i] += int(bool(self.get_overlap((sa_start, sa_end), bin)))
        
        # Number of reads including clipped reads whose aligned segment falls outside current bin
        all_reads_extended = len(all_reads_extended_ids)
         
        if all_reads > 0:
            # Case that we have aligned reads in bin, add 0.1 to avoid -inf values
            
            # Normalize insert size
            insertsize_mean = np.round(np.log2((np.mean(insert_sizes) + 0.1) / (self.baseline_insertsize_median + 0.1)), 3)
            insertsize_std = np.round(np.log2((np.std(insert_sizes) + 0.1) / (self.baseline_insertsize_mad + 0.1)), 3)

            # Normalize mapping quality
            mapping_quality_mean = np.round(np.log2((np.mean(mapqs) + 0.1) / (self.baseline_mapq_mean + 0.1)), 3)
            mapping_quality_std = np.round(np.log2((np.std(mapqs) + 0.1) / (self.baseline_mapq_std + 0.1)), 3)
            
            # Normalize discordant read pair features
            disco_ff_proportion = np.round(disco_ff_reads / all_reads, 3)
            disco_ff_conn_proportion = np.round(disco_ff_conn / disco_ff_reads, 3) if disco_ff_reads > 0 else disco_ff_conn
            disco_rr_proportion = np.round(disco_rr_reads / all_reads, 3)
            disco_rr_conn_proportion = np.round(disco_rr_conn / disco_rr_reads, 3) if disco_rr_reads > 0 else disco_rr_conn
            disco_rf_proportion = np.round(disco_rf_reads / all_reads, 3)
            disco_rf_conn_proportion = np.round(disco_rf_conn / disco_rf_reads, 3) if disco_rf_reads > 0 else disco_rf_conn
            disco_tx_proportion = np.round(disco_tx_reads / all_reads, 3)

        else:
            # Case that we have no aligned reads in bin
            
            # Set insert size to 0
            insertsize_mean = np.round(np.log2(0.1 / (self.baseline_insertsize_median + 0.1)), 3)
            insertsize_std = np.round(np.log2(0.1 / (self.baseline_insertsize_mad + 0.1)), 3)

            # Set mapping quality to 0
            mapping_quality_mean = np.round(np.log2(0.1 / (self.baseline_mapq_mean + 0.1)), 3)
            mapping_quality_std = np.round(np.log2(0.1 / (self.baseline_mapq_std + 0.1)), 3)
            
            # Set discordant read pair features to 0
            disco_ff_proportion = 0
            disco_ff_conn_proportion = disco_ff_conn
            disco_rr_proportion = 0
            disco_rr_conn_proportion = disco_rr_conn
            disco_rf_proportion = 0
            disco_rf_conn_proportion = disco_rf_conn
            disco_tx_proportion = 0
            
        if all_reads_extended > 0:
            # Case that we have reads in the extended bin, possibly corresponding to clipped reads
            
            # Normalize split-reads
            splitreads_proportion = np.round(split_reads / all_reads_extended, 3)
            if split_reads > 0:
                split_reads_conn_proportion = np.round(split_reads_conn / split_reads, 3)
            else:
                split_reads_conn_proportion = split_reads_conn
            
            # Normalize clipped reads
            clippedreads_proportion = np.round(clipped_reads / all_reads_extended, 3)
            
        else:
            # Case that we have no reads in the extended bin
            
            # Set split-reads to 0
            splitreads_proportion = 0
            split_reads_conn_proportion = split_reads_conn
            
            # Set clipped reads to 0
            clippedreads_proportion = 0
        
        # Create output pandas Series
        values = [insertsize_mean, insertsize_std, mapping_quality_mean, mapping_quality_std, 
                splitreads_proportion, clippedreads_proportion, disco_ff_proportion, disco_rr_proportion, 
                disco_rf_proportion, disco_tx_proportion] + list(split_reads_conn_proportion) + list(disco_ff_conn_proportion) + list(disco_rr_conn_proportion) + list(disco_rf_conn_proportion)
        index = ['ill_isize_mean_' + suffix, 'ill_isize_std_' + suffix, 'ill_mapq_mean_' + suffix, 'ill_mapq_std_' + suffix, 
                'ill_splitreads_' + suffix, 'ill_clipreads_' + suffix, 'ill_disco_ff_' + suffix, 'ill_disco_rr_' + suffix, 
                'ill_disco_rf_' + suffix, 'ill_disco_tx_' + suffix] + labels_split_reads_conn + labels_disco_conn
        
        return pd.Series(values, index=index)
        
        
    def annotate_read_based_features(self):
        """ Annotates SVs with read-based features for all bins around the SV breakpoints. """
        
        self._log(f'Starting read-based feature annotation: {len(self.df_calls_annot)} variants')
        start_total = time.time()

        for suffix, values in self.bin_dict_bps.items():
            start_bin = time.time()
            read_based_features = [feature + suffix for feature in self.features_breakpoints[2:]]
            connection_features = [feature + suffix + '_' + conn for feature in self.features_connection for conn in self.bin_connections[suffix]]
            self.df_calls_annot.loc[:, read_based_features + connection_features] = self.df_calls_annot.apply(lambda x: self.calculate_read_based_features(x[values[0]], 
                                                                                                                                x[values[1]] + values[3], 
                                                                                                                                x[values[2]] + values[4],
                                                                                                                                x['start'],
                                                                                                                                x['end'],
                                                                                                                                suffix), axis=1, result_type ='expand')
            self._log(f'  Bin {suffix} DONE: {len(self.df_calls_annot)} variants in {time.time() - start_bin:>7.1f}s', level='BIN')

        self._log(f'DONE: Read-based feature annotation completed in {time.time() - start_total:.1f}s', level='DONE')

     
    def aggregate_results(self):
        """ Cleans up result DataFrame. """        
            
        alignment_columns = []
        
        for feature in self.features_breakpoints:
            for suffix in ['I', 'II', 'III', 'IV']:
                alignment_columns.append(feature + suffix)
                
        for feature in self.features_body:
            for suffix in ['IIa', 'IIb', 'IIIb', 'IIIa']:
                alignment_columns.append(feature + suffix)
                
        for feature in self.features_connection:
            suffices = ['I', 'II', 'III', 'IV']
            for i in range(4):
                for j in range(i+1, 4):
                    alignment_columns.append(feature + suffices[i] + '_' + suffices[j])
        
        columns_reordered = ['id', 'cohort', 'sample', 'reference', 'technology', 'caller', 'sv_type', 'chrom', 'chrom_2', 'start', 'end'] + alignment_columns
        self.df_calls_annot = self.df_calls_annot[columns_reordered]
        self.bam.close()
        
        
    def to_csv(self, csv_filename: str):
        """ Writes the annotated SVs to a CSV file.

        Args:
            csv_filename (str): Path to the output file
        """        

        self.df_calls_annot.to_csv(csv_filename, index=False, na_rep='NA', sep='\t')
        
        
    def to_df(self) -> pd.DataFrame:
        """ Returns the annotated SVs as a DataFrame.

        Returns:
            pd.DataFrame: Result DataFrame
        """         
        
        return self.df_calls_annot
        

