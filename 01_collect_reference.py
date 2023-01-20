from lib.utils import read_parameters
from lib.feature_collection_reference import ReferenceAnnotator
import sys

SAMPLE = sys.argv[1]
REF = sys.argv[2]
INFILE = sys.argv[3]
OUTFILE = sys.argv[4]
PARAMS_FILE = sys.argv[5]

PARAMS = read_parameters(PARAMS_FILE, SAMPLE, REF)

RA = ReferenceAnnotator(PARAMS, INFILE, OUTFILE, SAMPLE, REF)

print('Repeats')
RA.annotate_repeats()
print('VNTRs')
RA.annotate_vntrs()
print('STRs')
RA.annotate_strs()
print('CGIs')
RA.annotate_cpg_islands()
print('Centromeres')
RA.annotate_centromeres()
print('Assembly Gaps')
RA.annotate_asmb_gaps()
print('Alternative Haplotypes')
RA.annotate_alt_haps()
print('GC Content')
RA.annotate_gc_content()

RA.to_csv()