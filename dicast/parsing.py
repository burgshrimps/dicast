import argparse
import os
import sys

from dicast import annotations

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Trained models ship inside the package so pip installs find them too.
DEFAULT_MODELS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'models')

# Canonical filenames of the hg38 annotation files shipped under --annot-dir.
ANNOT_CANONICAL_NAMES = {
    'repeats': 'hg38_repeatmasker.tsv',
    'cgis': 'hg38_cpg_islands.tsv',
    'centromeres': 'hg38_centromeres.tsv',
    'gaps': 'hg38_asmb_gaps.tsv',
    'althaps': 'hg38_alt_haps.tsv',
    'vntrs': 'hg38_vntrs_chaisson.bed',
    'strs': 'hg38_strs_chaisson.bed',
    'gc': 'hg38_gc_content.bw',
}

POP_CATALOG_NAME = 'pav_catalog_hg38.vcf.gz'


def _parse_multi_vcf_token(token: str) -> list:
    """ Parses a 'sample:caller=vcf_file' --vcfs token used by the 'multi' subcommand. """

    sample_and_caller, _, vcf_file = token.partition('=')
    sample, _, caller = sample_and_caller.partition(':')
    return [sample, caller, vcf_file]


def resolve_annotation_paths(args: argparse.Namespace) -> argparse.Namespace:
    """ Fills in unset annotation file paths from --annot-dir using canonical hg38 filenames.

    An explicitly passed flag always wins; --annot-dir is only used to fill in gaps.

    Args:
        args (argparse.Namespace): Parsed command line arguments

    Returns:
        argparse.Namespace: Arguments with annotation paths resolved
    """

    if args.command not in ('call', 'multi'):
        return args

    defaulted_flags = []
    if args.annot_dir is None:
        args.annot_dir = str(annotations.annot_dir(ask=True))
        defaulted_flags = [flag for flag in ANNOT_CANONICAL_NAMES if getattr(args, flag) is None]

    for flag, filename in ANNOT_CANONICAL_NAMES.items():
        if getattr(args, flag) is None:
            setattr(args, flag, os.path.join(args.annot_dir, filename))

    # Files the store manages are downloaded on first use (release assets, not
    # shipped in git); explicitly passed paths are the user's responsibility.
    if defaulted_flags:
        annotations.ensure_annotations(
            [ANNOT_CANONICAL_NAMES[flag] for flag in defaulted_flags], args.annot_dir)

    if args.pop_catalog is None:
        args.pop_catalog = os.path.join(args.annot_dir, POP_CATALOG_NAME)

    if args.pop:
        if args.command == 'call':
            args.vcfs = (args.vcfs or []) + [['pav', args.pop_catalog]]
        else:
            samples = sorted({sample for sample, _ in (args.bams or [])})
            args.vcfs = (args.vcfs or []) + [[sample, 'pav', args.pop_catalog] for sample in samples]

    return args


def _check_annotation_files(arguments: argparse.Namespace, problems: list):
    """ Appends a problem for every annotation flag that points at a missing file. """

    for flag, filename in ANNOT_CANONICAL_NAMES.items():
        path = getattr(arguments, flag)
        if path is None or os.path.isfile(path):
            continue
        problems.append(f'--{flag} file not found: {path} '
                        '(run: dicast-fetch-annotations, or point --annot-dir at your copies)')


def _check_model_files(arguments: argparse.Namespace, problems: list):
    """ Appends a problem for every SV type whose model file is missing. """

    sv_types_to_check = arguments.sv_types if arguments.sv_types else ['DEL', 'DUP', 'INS']
    for sv_type in sv_types_to_check:
        if arguments.pop and sv_type in ('DEL', 'INS') and os.path.isfile(os.path.join(arguments.models, f'dicast_{sv_type}_pop.json')):
            continue
        model_filename = os.path.join(arguments.models, f'dicast_{sv_type}.json')
        if not os.path.isfile(model_filename):
            problems.append(f'Model file not found for SV type {sv_type}: {model_filename}')


def make_workdir_tree(root: str):
    """ Creates the fixed input/features/output subtree under one sample's
    workdir root (== --workdir for `call`, == --workdir/<sample> for `multi`).

    Args:
        root (str): Root working directory for one sample.
    """

    os.makedirs(os.path.join(root, 'input'), exist_ok=True)
    os.makedirs(os.path.join(root, 'features', 'ref'), exist_ok=True)
    os.makedirs(os.path.join(root, 'features', 'aln'), exist_ok=True)
    os.makedirs(os.path.join(root, 'output'), exist_ok=True)


def _check_bam_index(bam_file: str, problems: list):
    """ Appends a problem if no index can be found next to a BAM file. """

    bai_candidates = [bam_file + '.bai', bam_file + '.csi']
    if bam_file.endswith('.bam'):
        bai_candidates.append(bam_file[:-len('.bam')] + '.bai')
    if not any(os.path.isfile(candidate) for candidate in bai_candidates):
        problems.append(f'BAM index not found for {bam_file} (run: samtools index {bam_file})')


def validate_inputs(arguments: argparse.Namespace):
    """ Validates inputs, printing every problem found and exiting on failure.

    Args:
        arguments (argparse.Namespace): Parsed command line arguments
    """

    if arguments.command == 'call':
        _validate_call_inputs(arguments)
    elif arguments.command == 'multi':
        _validate_multi_inputs(arguments)


def _validate_call_inputs(arguments: argparse.Namespace):
    """ Validates inputs of the 'call' subcommand. """

    problems = []

    def check_file(path, label):
        if path is None:
            problems.append(f'{label} is not set.')
        elif not os.path.isfile(path) or not os.access(path, os.R_OK):
            problems.append(f'{label} not found or not readable: {path}')

    check_file(arguments.fai, '--fai file')
    check_file(arguments.bam, '--bam file')

    _check_annotation_files(arguments, problems)

    if not arguments.vcfs:
        problems.append('--vcfs is required (format: caller=vcf_file).')
    else:
        for token in arguments.vcfs:
            if len(token) != 2:
                problems.append(f"Malformed --vcfs entry (expected caller=vcf_file): {'='.join(token)}")
            elif token[0] == 'pav' and arguments.pop:
                continue  # covered by the --pop-catalog check below
            elif not os.path.isfile(token[1]):
                problems.append(f'VCF file for caller {token[0]} not found: {token[1]}')

    if arguments.pop and not os.path.isfile(arguments.pop_catalog):
        problems.append(f'--pop-catalog file not found: {arguments.pop_catalog} '
                        '(the PAV population catalog will be published with the release assets; '
                        'until then point --pop-catalog at your copy)')

    if arguments.bam and os.path.isfile(arguments.bam):
        _check_bam_index(arguments.bam, problems)

    _check_model_files(arguments, problems)

    if problems:
        for problem in problems:
            print(f'ERROR: {problem}', file=sys.stderr)
        sys.exit(1)

    make_workdir_tree(arguments.workdir)


def _validate_multi_inputs(arguments: argparse.Namespace):
    """ Validates inputs of the 'multi' subcommand. """

    problems = []

    def check_file(path, label):
        if path is None:
            problems.append(f'{label} is not set.')
        elif not os.path.isfile(path) or not os.access(path, os.R_OK):
            problems.append(f'{label} not found or not readable: {path}')

    check_file(arguments.fai, '--fai file')

    if not arguments.bams:
        problems.append('--bams is required (format: sample=bam_file).')

    bam_samples = set()
    for token in (arguments.bams or []):
        if len(token) != 2:
            problems.append(f"Malformed --bams entry (expected sample=bam_file): {'='.join(token)}")
            continue
        sample, bam_file = token
        bam_samples.add(sample)
        check_file(bam_file, f'--bams BAM file for sample {sample}')
        if os.path.isfile(bam_file):
            _check_bam_index(bam_file, problems)

    if len(bam_samples) < 2:
        problems.append('--bams must specify at least two samples for multi-sample rescue mode.')

    if not arguments.vcfs:
        problems.append('--vcfs is required (format: sample:caller=vcf_file).')

    vcf_samples = set()
    for sample, caller, vcf_file in (arguments.vcfs or []):
        if not sample or not caller or not vcf_file:
            problems.append(f"Malformed --vcfs entry (expected sample:caller=vcf_file): {sample}:{caller}={vcf_file}")
            continue
        vcf_samples.add(sample)
        if caller == 'pav' and arguments.pop:
            continue  # covered by the --pop-catalog check below
        if not os.path.isfile(vcf_file):
            problems.append(f'VCF file for sample {sample}, caller {caller} not found: {vcf_file}')

    missing_vcfs = bam_samples - vcf_samples
    if missing_vcfs:
        problems.append(f'No --vcfs entries for sample(s) with a --bams entry: {", ".join(sorted(missing_vcfs))}')
    unknown_bams = vcf_samples - bam_samples
    if unknown_bams:
        problems.append(f'--vcfs references sample(s) without a --bams entry: {", ".join(sorted(unknown_bams))}')

    if arguments.pop and not os.path.isfile(arguments.pop_catalog):
        problems.append(f'--pop-catalog file not found: {arguments.pop_catalog} '
                        '(the PAV population catalog will be published with the release assets; '
                        'until then point --pop-catalog at your copy)')

    _check_annotation_files(arguments, problems)
    _check_model_files(arguments, problems)

    if problems:
        for problem in problems:
            print(f'ERROR: {problem}', file=sys.stderr)
        sys.exit(1)

    os.makedirs(arguments.workdir, exist_ok=True)
    for sample in bam_samples:
        make_workdir_tree(os.path.join(arguments.workdir, sample))


def parse_arguments(arguments = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='dicast')

    subparsers = parser.add_subparsers(dest='command', help='Subcommands')

    parser_call = subparsers.add_parser('call', help='Feature extraction and variant prediction')
    parser_call.add_argument('--cohort', help='Cohort name', default='none')
    parser_call.add_argument('--sample', help='Sample name', required=True)
    parser_call.add_argument('--chrom', help='Chromosomes', nargs='+', default='all')
    parser_call.add_argument('--ref', help='Reference genome name', default='hg38')
    parser_call.add_argument('--technology', help='Sequencing technology name', default='ill')
    parser_call.add_argument('--workdir', help='Working and output directory', required=True)
    parser_call.add_argument('--fai', help='FAI file of the reference genome', required=True)
    parser_call.add_argument('--annot-dir', help='Directory with the hg38 annotation files (default: the dicast data directory, downloaded on first use)', default=None)
    parser_call.add_argument('--repeats', help='TSV file with repeats annotated by repeatmasker', default=None)
    parser_call.add_argument('--cgis', help='TSV file with CpG island annotations', default=None)
    parser_call.add_argument('--centromeres', help='TSV file with centromere annotations', default=None)
    parser_call.add_argument('--gaps', help='TSV file with assembly gap annotations', default=None)
    parser_call.add_argument('--althaps', help='TSV file with alternative haplotype annotations', default=None)
    parser_call.add_argument('--vntrs', help='BED file with VNTR regions from Chaisson', default=None)
    parser_call.add_argument('--strs', help='BED file with STR regions from Chaisson', default=None)
    parser_call.add_argument('--gc', help='BIGWIG file with GC content', default=None)
    parser_call.add_argument('--bam', help='BAM file', required=True)
    parser_call.add_argument('--vcfs', nargs='*', type=lambda kv: kv.split('='), help='List of VCF files. Needs to be in the format method=vcf_file', required=True)
    parser_call.add_argument('--models', help='Directory with trained models', default=DEFAULT_MODELS_DIR)
    parser_call.add_argument('--threads', help='Number of threads to use', default=1, type=int)
    parser_call.add_argument('--pop', help='Add the PAV population catalog as an additional caller and prefer population-aware models', action='store_true')
    parser_call.add_argument('--pop-catalog', help='VCF file with the PAV population catalog', default=None)
    parser_call.add_argument('--benchmark', help='Path to write a per-stage TSV with wall-time, CPU-time and peak RSS (feature_collection, prediction, total). If unset, no benchmark is written.', default=None)
    parser_call.add_argument('--sv_types', help='Restrict feature extraction and prediction to these SV types. Default: DEL DUP INS.', nargs='+', choices=['DEL', 'DUP', 'INS'], default=None)

    parser_multi = subparsers.add_parser('multi', help='Multi-sample scoring with cross-sample rescue (e.g. a trio)')
    parser_multi.add_argument('--cohort', help='Cohort name', default='none')
    parser_multi.add_argument('--bams', nargs='*', type=lambda kv: kv.split('='), help='Per-sample BAM files. Format: sample=bam_file', required=True)
    parser_multi.add_argument('--vcfs', nargs='*', type=_parse_multi_vcf_token, help='Per-sample, per-caller VCF files. Format: sample:caller=vcf_file', required=True)
    parser_multi.add_argument('--chrom', help='Chromosomes', nargs='+', default='all')
    parser_multi.add_argument('--ref', help='Reference genome name', default='hg38')
    parser_multi.add_argument('--technology', help='Sequencing technology name', default='ill')
    parser_multi.add_argument('--workdir', help='Working and output directory', required=True)
    parser_multi.add_argument('--fai', help='FAI file of the reference genome', required=True)
    parser_multi.add_argument('--annot-dir', help='Directory with the hg38 annotation files (default: the dicast data directory, downloaded on first use)', default=None)
    parser_multi.add_argument('--repeats', help='TSV file with repeats annotated by repeatmasker', default=None)
    parser_multi.add_argument('--cgis', help='TSV file with CpG island annotations', default=None)
    parser_multi.add_argument('--centromeres', help='TSV file with centromere annotations', default=None)
    parser_multi.add_argument('--gaps', help='TSV file with assembly gap annotations', default=None)
    parser_multi.add_argument('--althaps', help='TSV file with alternative haplotype annotations', default=None)
    parser_multi.add_argument('--vntrs', help='BED file with VNTR regions from Chaisson', default=None)
    parser_multi.add_argument('--strs', help='BED file with STR regions from Chaisson', default=None)
    parser_multi.add_argument('--gc', help='BIGWIG file with GC content', default=None)
    parser_multi.add_argument('--models', help='Directory with trained models', default=DEFAULT_MODELS_DIR)
    parser_multi.add_argument('--threads', help='Number of threads to use', default=1, type=int)
    parser_multi.add_argument('--pop', help='Add the PAV population catalog as an additional caller (for every sample) and prefer population-aware models', action='store_true')
    parser_multi.add_argument('--pop-catalog', help='VCF file with the PAV population catalog', default=None)
    parser_multi.add_argument('--benchmark', help='Path to write a per-stage TSV with wall-time, CPU-time and peak RSS (feature_collection, prediction, total). If unset, no benchmark is written.', default=None)
    parser_multi.add_argument('--sv_types', help='Restrict feature extraction and prediction to these SV types. Default: DEL DUP INS.', nargs='+', choices=['DEL', 'DUP', 'INS'], default=None)

    args = parser.parse_args(arguments)
    args = resolve_annotation_paths(args)

    return args
