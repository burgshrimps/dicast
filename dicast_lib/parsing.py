import argparse
import sys

def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='dicast')

    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    parser_call = subparsers.add_parser('call', help='Feature extraction and variant prediction')
    parser_call.add_argument('--cohort', help='Cohort name')
    parser_call.add_argument('--sample', help='Sample name')
    parser_call.add_argument('--ref', help='Reference genome name', default='hg38')
    parser_call.add_argument('--technology', help='Sequencing technology name', default='ill')
    parser_call.add_argument('--workdir', help='Working and output directory')
    parser_call.add_argument('--fai', help='FAI file of the reference genome')
    parser_call.add_argument('--repeats', help='TSV file with repeats annotated by repeatmasker')
    parser_call.add_argument('--cgis', help='TSV file with CpG island annotations')
    parser_call.add_argument('--centromeres', help='TSV file with centromere annotations')
    parser_call.add_argument('--gaps', help='TSV file with assembly gap annotations')
    parser_call.add_argument('--althaps', help='TSV file with alternative haplotype annotations')
    parser_call.add_argument('--vntrs', help='BED file with VNTR regions from Chaisson')
    parser_call.add_argument('--strs', help='BED file with STR regions from Chaisson')
    parser_call.add_argument('--gc', help='BIGWIG file with GC content')
    parser_call.add_argument('--bam', help='BAM file')
    parser_call.add_argument('--vcfs', nargs='*', type=lambda kv: kv.split('='), help='List of VCF files. Needs to be in the format method=vcf_file')
    parser_call.add_argument('--models', help='Directory with trained models')
    parser_call.add_argument('--threads', help='Number of threads to use', default=1, type=int)

    return parser.parse_args(arguments)