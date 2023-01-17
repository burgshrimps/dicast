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
RA.annotate_repeats()
RA.annotate_vntrs()
RA.annotate_strs()
RA.annotate_cpg_islands()
RA.annotate_centromeres()
RA.annotate_asmb_gaps()
RA.annotate_alt_haps()
RA.annotate_gc_content()
RA.to_csv()