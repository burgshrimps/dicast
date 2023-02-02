import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.colors import ListedColormap
from matplotlib import gridspec
import pandas as pd
import pysam
import os


def read_parameters(file):
    """ Reads parameter file and saves them to dictionary. """

    with open(file,'r') as f:
        params = json.load(f)
    return params


def replace_filename(filename, sample, ref):
    """ Replaces sample name in filename. """

    return filename.replace('SAMPLE', sample).replace('REF', ref)


def read_vcf(filename):
    """ Reads VCF file as a pysam.VariantFile object and takes care of sample names. """

    if os.path.exists(filename):
        vcf = pysam.VariantFile(filename)
    else:
        # Smoove appears in lumpy filenames
        vcf = pysam.VariantFile(filename.replace('-', '_').replace('_smoove', '-smoove'))

    return vcf


def compute_overlap(s1, s2, e1, e2):
    """ Computes overlap between two segments. """
    return max(0, min(e1, e2) - max(s1, s2))


def cigartuples_to_array(cigartuples):
    """ Converts list of cigar tuples to numpy array. """
    cigar_array = []
    for cigartuple in cigartuples:
        cigar_array += [cigartuple[0]] * cigartuple[1]
    return np.array(cigar_array)


def compute_aln_matrix(bam, chrom, start, stop, size=100):
    """ Computes alignment matrix consisting of CIGAR integers for a given region. """
    reads = []
    for read in bam.fetch(chrom, start-5, stop+5):
        if not read.is_unmapped:
            reads.append(read)

    aln_matrix = -1 * np.ones((len(reads), size))
    split_reads_idx = []
    low_mapq_idx = []

    for idx, read in enumerate(reads):
        if read.has_tag('SA'):
            split_reads_idx.append(idx)
        if read.mapping_quality < 5:
           low_mapq_idx.append(idx)

        cigararray = cigartuples_to_array(read.cigartuples)

        if read.cigartuples[0][0] == 4 or read.cigartuples[0][0] == 5:
            read_start = read.reference_start - read.cigartuples[0][1]
        else:
            read_start = read.reference_start
        if read.cigartuples[-1][0] == 4 or read.cigartuples[-1][0] == 5:
            read_end = read.reference_end + read.cigartuples[-1][1]
        else:
            read_end = read.reference_end

        del_idx = []
        for i in range(len(cigararray) -1):
            if np.all(np.equal(cigararray[i:i+2], np.array([0,1]))):
                cigararray[i] = 1
            elif np.all(np.equal(cigararray[i:i+2], np.array([1,1]))):
                del_idx.append(i)
            elif np.all(np.equal(cigararray[i:i+2], np.array([1,0]))):
                del_idx.append(i)

        cigararray = np.delete(cigararray, del_idx)
        
        if read_start < start:
            start_idx_read = start - read_start
            start_idx_read = start_idx_read
            start_idx_aln = 0
        else:
            start_idx_read = 0
            start_idx_aln = read_start - start
        if read_end > stop:
            end_idx_read = stop - read_start
            end_idx_aln = size
        else:
            end_idx_read = len(cigararray)
            end_idx_aln = read_end - start
        
        if end_idx_aln > 0 and end_idx_read > 0:
            try:
                aln_matrix[idx, start_idx_aln:end_idx_aln] = cigararray[start_idx_read:end_idx_read]
            except:
                print(start_idx_aln, end_idx_aln, start_idx_read, end_idx_read)
    
    return aln_matrix, split_reads_idx, low_mapq_idx


def pad_alignment_matrices(aln_matrix_left, aln_matrix_right):
    """ Pads alignment matrices with zeros to make them the same height. """

    if len(aln_matrix_right) < len(aln_matrix_left):
        aln_matrix_right = np.pad(aln_matrix_right, ((0, len(aln_matrix_left) - len(aln_matrix_right)), (0, 0)), 'constant', constant_values=(-1, -1))
    elif len(aln_matrix_right) > len(aln_matrix_left):
        aln_matrix_left = np.pad(aln_matrix_left, ((0, len(aln_matrix_right) - len(aln_matrix_left)), (0, 0)), 'constant', constant_values=(-1, -1))
    return aln_matrix_left, aln_matrix_right


def compute_cov_df(params, chrom, start, stop, minq=30):
    """ Computes coverage for given region. """

    alignment_dir = replace_filename(params['alignments']['illumina_directory'], params)
    alignment_filename = replace_filename(params['alignments']['illumina_filename'], params).replace('-', '_')
    cov = pd.DataFrame([x.split('\t') for x in pysam.depth('/'.join([alignment_dir, alignment_filename]), '-r', chrom + ':' + str(start) + '-' + str(stop), '-a').split('\n')[:-1]])
    cov_minq = pd.DataFrame([x.split('\t') for x in pysam.depth('/'.join([alignment_dir, alignment_filename]), '-r', chrom + ':' + str(start) + '-' + str(stop), '-a', '-Q', str(minq)).split('\n')[:-1]])
    cov[1] = cov[1].astype(int)
    cov[2] = cov[2].astype(int)
    cov[1] = cov[1] - cov.loc[0, 1]
    cov_minq[1] = cov_minq[1].astype(int)
    cov_minq[2] = cov_minq[2].astype(int)
    cov_minq[1] = cov_minq[1] - cov_minq.loc[0, 1]
    return cov, cov_minq


def compute_rep_df(params, chrom, start, stop, padding=500):
    """ Computes repeat overlap for a given region. """

    ref_dir_root = replace_filename(params['reference']['directory'], params)
    ref_dir_annot = params['reference']['subdirectory_annotation']
    ref_dir = '/'.join([ref_dir_root, ref_dir_annot])
    rep_filename = replace_filename(params['reference']['filename_repeats'], params).replace('-', '_')
    rep_df = pd.read_csv('/'.join([ref_dir, rep_filename]), index_col=0, sep='\t')
    rep_df = rep_df[(rep_df['genoName'] == chrom) & (rep_df['genoStart'] > start) & (rep_df['genoEnd'] < stop)].copy()
    rep_df['genoStart'] = rep_df['genoStart'] - start + padding
    rep_df['genoEnd'] = rep_df['genoEnd'] - start + padding
    rep_df = rep_df[(rep_df['genoEnd'] > 0) & (rep_df['genoStart'] < stop - start + (2*padding))]
    rep_df.reset_index(drop=True, inplace=True)
    return rep_df


def mad(arr):
    """ Median absolute deviation. """
    
    med = np.median(arr)
    return np.median(np.abs(arr - med))