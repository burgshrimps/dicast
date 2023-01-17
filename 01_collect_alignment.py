import sys

from lib.utils import read_parameters
from lib.feature_collection_illumina import AlignmentAnnotatorIllumina

SAMPLE = sys.argv[1]
REF = sys.argv[2]
INFILE = sys.argv[3]
OUTFILE = sys.argv[4]
PARAMS_FILE = sys.argv[5]
CHROM = sys.argv[6]

PARAMS = read_parameters(PARAMS_FILE, SAMPLE, REF)

AAI = AlignmentAnnotatorIllumina(PARAMS, INFILE, OUTFILE, SAMPLE, REF, chrom=CHROM)
AAI.calculate_coverage_baseline()
AAI.calculate_insertsize_baseline()
AAI.calculate_mapping_quality_baseline()
AAI.annotate_coverage()
AAI.annotate_read_based_features()
AAI.to_csv()