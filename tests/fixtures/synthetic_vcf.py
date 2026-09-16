"""Parametrisable VCF/FAI text writers for ``tests/unit/test_vcf_input.py``.

``dicast/vcf_input.py`` reads real caller VCFs with wildly different
dialects (declared vs. undeclared contigs/filters, symbolic vs.
sequence-resolved ALTs, breakend bracket notation, 0/1/many sample columns,
tuple-valued SVLEN...), so the tests need fine control over the exact text of
each record rather than a single fixed header + record set. :func:`write_vcf`
takes a list of plain record dicts and writes them out as a minimal, valid
VCF; :func:`write_fai` writes a matching ``.fai``.

Two VCF-coordinate facts every EXPECTED value in the tests is derived from
(never from running the parser first):

* a VCF POS of ``P`` is 1-based; ``rec.pos == P`` and ``rec.start == P - 1``.
* for a *sequence-resolved* (non-symbolic) ALT, htslib sets
  ``rec.stop == rec.start + len(REF)`` even with no INFO/END -- this is what
  lets a sequence-resolved DEL fall through vcf_input's normal
  ``end = rec.stop`` DEL/DUP path with no special-casing.
"""
from __future__ import annotations

import pathlib
from typing import Dict, List, Optional, Sequence, Tuple

import pysam

# INFO/FORMAT fields declared by default so a test only needs `extra_header`
# for something genuinely unusual (a custom FILTER, a pre-existing DQ tag, an
# explicit ##ALT line, ...). SVLEN is declared Number=. on purpose: real
# caller headers do this, and pysam then hands SVLEN back as a tuple, which
# is exactly the case parse_svlen has to unwrap.
_BASE_HEADER_LINES = [
    '##fileformat=VCFv4.2',
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
    '##INFO=<ID=SIMPLE_TYPE,Number=1,Type=String,Description="gridss simple-event-annotation.R type">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">',
    '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


def _record_fields(rec: dict, sample_names: Sequence[str], has_format_col: bool) -> List[str]:
    """Renders one record dict as the tab-separated VCF data-line fields."""
    fields = [
        rec['chrom'],
        str(rec['pos']),
        rec.get('id', '.'),
        rec.get('ref', 'N'),
        rec['alt'],
        str(rec.get('qual', '.')),
        rec.get('filter', 'PASS'),
        rec.get('info', '.'),
    ]
    if has_format_col:
        fields.append(rec.get('format', 'GT'))
        rec_samples = rec.get('samples', {})
        for name in sample_names:
            fields.append(rec_samples.get(name, './.'))
    return fields


def write_vcf(tmp_path, name: str, records: List[dict], samples: Optional[Sequence[str]] = None,
              contigs: Optional[Sequence[Tuple[str, Optional[int]]]] = None,
              extra_header: Optional[Sequence[str]] = None, bgzip: bool = False) -> str:
    """Writes a minimal VCF from plain record dicts and returns its path.

    Args:
        tmp_path: a directory (e.g. pytest's ``tmp_path``) to write into.
        name (str): file basename, e.g. ``'caller.vcf'``.
        records (list): one dict per record. Recognized keys: ``chrom``,
            ``pos`` (required); ``id`` (default ``'.'``), ``ref`` (default
            ``'N'``), ``alt`` (required), ``qual`` (default ``'.'``),
            ``filter`` (default ``'PASS'``), ``info`` (raw INFO text,
            default ``'.'``), ``format`` (default ``'GT'``, only written
            when a FORMAT column exists at all -- see `samples` below), and
            ``samples`` (dict of sample name -> that record's raw FORMAT
            value string, e.g. ``{'S1': '0/1'}``; a sample missing from this
            dict gets ``'./.'``).
        samples (list, optional): sample column names for the ``#CHROM``
            line. A FORMAT column (and per-sample columns) is written only
            when `samples` is non-empty; passing none produces a sites-only
            VCF, exercising the "0 samples" contract case.
        contigs (list, optional): ``(name, length)`` pairs for ``##contig``
            lines; `length` may be None to omit the length attribute. Tests
            that want an *undeclared* contig simply omit it here while still
            using that chrom name in a record -- vcf_input discovers contigs
            from records, not from the header.
        extra_header (list, optional): additional raw ``##...`` lines
            (custom FILTER/INFO/ALT declarations, a pre-existing DQ tag...).
        bgzip (bool): if True, bgzip + tabix-index the file and return the
            ``.gz`` path (tabix requires `records` to already be given in
            ascending (contig-group, POS) order); if False (the default),
            write plain text, which pysam reads with no ordering
            requirement at all -- the right choice whenever a test wants
            record order that is not also chrom/pos order.

    Returns:
        str: the written file's path (``.gz`` when ``bgzip=True``).
    """
    samples = list(samples or [])
    contigs = list(contigs or [])
    extra_header = list(extra_header or [])

    lines = list(_BASE_HEADER_LINES)
    for contig_name, contig_length in contigs:
        if contig_length is None:
            lines.append(f'##contig=<ID={contig_name}>')
        else:
            lines.append(f'##contig=<ID={contig_name},length={contig_length}>')
    lines.extend(extra_header)

    has_format_col = bool(samples)
    header_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']
    if has_format_col:
        header_cols.append('FORMAT')
        header_cols.extend(samples)
    lines.append('\t'.join(header_cols))

    for rec in records:
        lines.append('\t'.join(_record_fields(rec, samples, has_format_col)))

    path = pathlib.Path(tmp_path) / name
    path.write_text('\n'.join(lines) + '\n')

    if not bgzip:
        return str(path)
    return pysam.tabix_index(str(path), preset='vcf', force=True)


def write_fai(tmp_path, name: str, contigs: Sequence[Tuple[str, int]]) -> str:
    """Writes a minimal ``.fai`` (name, length, dummy offset/linebases/
    linewidth columns) and returns its path.

    Args:
        tmp_path: a directory to write into.
        name (str): file basename, e.g. ``'ref.fa.fai'``.
        contigs (list): ``(name, length)`` pairs, in file order.

    Returns:
        str: the written file's path.
    """
    path = pathlib.Path(tmp_path) / name
    lines = [f'{contig_name}\t{length}\t0\t60\t61' for contig_name, length in contigs]
    path.write_text('\n'.join(lines) + '\n')
    return str(path)


# ---------------------------------------------------------------------------
# MIXED_RECORDS: one record set exercising most of read_caller_vcf's
# branches at once, with expected outcomes hand-derived record by record
# (never from running the parser). Used by the small number of integration-
# style tests that check the whole per-file FileReport tally; individual
# branches each get their own minimal, single- or two-record VCF elsewhere
# in test_vcf_input.py.
#
# Every chrom is written *without* a 'chr' prefix (numeric-style, like
# delly/lumpy GRCh37 output) against MIXED_CANONICAL = ['chr1', 'chr2'], so
# every kept record also exercises the '+chr' contig-mapping path.
# ---------------------------------------------------------------------------

MIXED_SAMPLES = ['S1', 'S2']
MIXED_CANONICAL = ['chr1', 'chr2']
# Deliberately no ##contig lines: MIXED_RECORDS also exercises discovering
# contigs from records alone (cnvnator2VCF / raw lumpy have none).
MIXED_CONTIGS: List[Tuple[str, Optional[int]]] = []

MIXED_RECORDS = [
    # DEL1: symbolic ALT, INFO END+SVLEN, S1 has the alt allele.
    {
        'chrom': '1', 'pos': 1000, 'id': 'DEL1', 'alt': '<DEL>',
        'info': 'SVTYPE=DEL;END=1500;SVLEN=-500',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # INS1: symbolic ALT, SVLEN only, but ABSENT in S1 (0/0) -- dropped by
    # the multi-sample "keep only if GT has an alt allele" policy.
    {
        'chrom': '1', 'pos': 2000, 'id': 'INS1', 'alt': '<INS>',
        'info': 'SVTYPE=INS;SVLEN=250',
        'samples': {'S1': '0/0', 'S2': '0/1'},
    },
    # INV1: out of contract -- dropped as INV regardless of S1 genotype.
    {
        'chrom': '1', 'pos': 3000, 'id': 'INV1', 'alt': '<INV>',
        'info': 'SVTYPE=INV;END=3800',
        'samples': {'S1': '1/1', 'S2': '0/0'},
    },
    # CNV1: out of contract -- dropped as CNV.
    {
        'chrom': '1', 'pos': 4000, 'id': 'CNV1', 'alt': '<CN0>',
        'info': 'SVTYPE=CNV',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # SEQDEL: no SVTYPE at all; sequence-resolved REF(61)/ALT(1) -> DEL via
    # the length-diff fallback, end = rec.stop = pos + len(REF) - 1.
    {
        'chrom': '1', 'pos': 5000, 'id': 'SEQDEL', 'ref': 'A' * 61, 'alt': 'A',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # SEQINS: no SVTYPE; sequence-resolved REF(1)/ALT(61) -> INS, sv_len =
    # len(ALT) - len(REF) = 60.
    {
        'chrom': '1', 'pos': 6000, 'id': 'SEQINS', 'ref': 'A', 'alt': 'A' + 'C' * 60,
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # SMALL: sequence-resolved, |diff| = 9 < 50 -> small_variant, dropped.
    {
        'chrom': '1', 'pos': 7000, 'id': 'SMALL', 'ref': 'ACGTACGTAC', 'alt': 'A',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # BNDPLAIN: SVTYPE=BND, no SIMPLE_TYPE -- out of contract, dropped BND.
    {
        'chrom': '1', 'pos': 8000, 'id': 'BNDPLAIN', 'alt': 'N[2:9000[',
        'info': 'SVTYPE=BND',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # GRIDSSLO / GRIDSSHI: a gridss breakend pair both annotated
    # SIMPLE_TYPE=DEL. GRIDSSLO is the lower breakend (pos 9000 < mate
    # 9500, same contig) -> kept as DEL, end = mate pos = 9500. GRIDSSHI is
    # the same event's other end (pos 9500 >= its mate 9000) -> dropped
    # bnd_mate so the event is not double-counted.
    {
        'chrom': '1', 'pos': 9000, 'id': 'GRIDSSLO', 'alt': 'N[1:9500[',
        'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    {
        'chrom': '1', 'pos': 9500, 'id': 'GRIDSSHI', 'alt': ']1:9000]N',
        'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # NOEND: symbolic DEL with neither END nor SVLEN -> missing_end.
    {
        'chrom': '1', 'pos': 10000, 'id': 'NOEND', 'alt': '<DEL>',
        'info': 'SVTYPE=DEL',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
    # CONTIGX: same as DEL1 but on the second (also undeclared, +chr-mapped)
    # contig, to prove both contigs get discovered and mapped.
    {
        'chrom': '2', 'pos': 100, 'id': 'CONTIGX', 'alt': '<DEL>',
        'info': 'SVTYPE=DEL;END=600;SVLEN=-500',
        'samples': {'S1': '0/1', 'S2': '0/0'},
    },
]

# Hand-derived per read_caller_vcf(..., sample='S1', canonical=MIXED_CANONICAL).
# Records read == len(MIXED_RECORDS) == 12.
MIXED_EXPECTED_KEPT = {'DEL': 4, 'DUP': 0, 'INS': 1}
MIXED_EXPECTED_DROPPED = {
    'contig': 0, 'no_svtype': 0, 'INV': 1, 'BND': 1, 'CNV': 1,
    'other_svtype': 0, 'small_variant': 1, 'missing_end': 1, 'bnd_mate': 1,
    'absent_in_sample': 1,
}
# ids expected to survive into the output dataframe, in MIXED_RECORDS order.
MIXED_KEPT_IDS = ['DEL1', 'SEQDEL', 'SEQINS', 'GRIDSSLO', 'CONTIGX']
