import pandas as pd
import pysam
import numpy as np
import re
import sys
from collections import defaultdict

from lib.utils import read_parameters, replace_filename

# PARAMETERS
SAMPLE = sys.argv[1]
REF = sys.argv[2]
INFILE = sys.argv[3]
OUTFILE = sys.argv[4]
PARAMS_FILE = sys.argv[5]

CHROMS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 
          'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
METHODS = ['delly', 'manta', 'lumpy']

params = read_parameters(PARAMS_FILE, SAMPLE, REF)


## FUNCTIONS
def parse_vcf(vcf, tech, method, sample):
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

def check_out_of_bounds(svtype, chrom, chrom2, start, end, chrom_sizes, padding=50):
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
        return start - padding < 0 or end + padding > df_chrom_sizes.loc[chrom, 'size']
    else:
        # For translocations, check both chromosomes
        outbounds_chrom1 = start - padding < 0 or start + padding > df_chrom_sizes.loc[chrom, 'size']
        outbounds_chrom2 = end - padding < 0 or end + padding > df_chrom_sizes.loc[chrom2, 'size']
        return outbounds_chrom1 or outbounds_chrom2


# SCRIPT
# Load data
reference_dir = replace_filename(params['reference']['directory'], params)

df_tgenvar = pd.read_csv(INFILE, low_memory=False, index_col=0).reset_index(drop=True)
df_chrom_sizes = pd.read_csv('/'.join([reference_dir, params['reference']['filename_chrom_sizes']]), sep='\t',
                             header=None, names=['size', 'offset', 'linebases', 'linewidth'], index_col=0)

filename_delly = replace_filename(params['variant_calls']['filename_calls_raw_delly'], params)
filename_manta = replace_filename(params['variant_calls']['filename_calls_raw_manta'], params)
filename_lumpy = replace_filename(params['variant_calls']['filename_calls_raw_lumpy'], params)

vcf_delly = pysam.VariantFile(filename_delly.replace('-', '_'))
vcf_manta = pysam.VariantFile(filename_manta.replace('-', '_'))
vcf_lumpy = pysam.VariantFile(filename_lumpy.replace('-', '_').replace('_smoove', '-smoove'))


# Parse VCFs
delly = parse_vcf(vcf_delly, 'ILL', 'delly', SAMPLE)
manta = parse_vcf(vcf_manta, 'ILL', 'manta', SAMPLE)
lumpy = parse_vcf(vcf_lumpy, 'ILL', 'lumpy', SAMPLE)
df_raw = pd.concat([delly, manta, lumpy], ignore_index=True)

# Remove calls that are located on non-canonical chromosomes
df_raw = df_raw[(df_raw['type'] == 'BND') | (df_raw['chrom'].isin(CHROMS))].copy().reset_index(drop=True)
df_raw.drop(df_raw[(df_raw['type'] == 'BND') & ((~df_raw['chrom'].isin(CHROMS)) | (~df_raw['chrom2'].isin(CHROMS)))].index, inplace=True)

# Remove calls that are located on non-canonical chromosomes
df_raw = df_raw[(df_raw['type'] == 'BND') | (df_raw['chrom'].isin(CHROMS))].copy().reset_index(drop=True)
df_raw.drop(df_raw[(df_raw['type'] == 'BND') & ((~df_raw['chrom'].isin(CHROMS)) | (~df_raw['chrom2'].isin(CHROMS)))].index, inplace=True)

# Remove calls that are out of chromosome bounds
df_raw['start'] = df_raw['start'].astype(int)
df_raw['end'] = df_raw['end'].astype(int)
df_raw['outbounds'] = df_raw.apply(lambda x: check_out_of_bounds(x['type'], x['chrom'], x['chrom2'], x['start'], x['end'], df_chrom_sizes), axis=1)
df_raw = df_raw[~df_raw['outbounds']].copy().drop('outbounds', axis=1).reset_index(drop=True)

# Extract IDs of confirmed SVs
df_tgenvar = df_tgenvar[(df_tgenvar['StatusSimple'] == 'Confirmed') | (df_tgenvar['StatusSimple'] == 'ConfirmedPublic')].copy().reset_index(drop=True)
confirmed_calls = {'delly' : [], 'manta' : [], 'lumpy' : []}

for i in range(len(df_tgenvar)):
    sub_graph = df_tgenvar.loc[i, 'sub_graph'][1:-1].split(', ')
    for entry in sub_graph:
        for method in METHODS:
            if entry[1:-1].startswith(method):
                confirmed_calls[method].append(entry[1:-1].split('_')[-1])

# Add confirmation status
df_raw['confirmed'] = 0
for i in range(len(df_raw)):
    if df_raw.loc[i, 'id'] in confirmed_calls[df_raw.loc[i, 'method']]:
        df_raw.loc[i, 'confirmed'] = 1 

# Save
df_raw.to_csv(OUTFILE, index=False, sep='\t', na_rep='NA')