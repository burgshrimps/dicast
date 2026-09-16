"""Single entry point for turning an arbitrary SV-caller VCF into the
dataframe rows dicast's model needs, plus a per-file report of what was kept
and what was dropped.

dicast was written against one caller's VCF dialect (`INFO/SVTYPE` present,
`chr`-prefixed hg38 contigs, unique non-empty record IDs, `END` present for
symbolic ALTs, a sample column named exactly like `--sample`). Real caller
output varies on every one of those points, so this module defines the
minimal contract dicast needs from a record --
``chrom, pos, id, sv_type in {DEL, DUP, INS}, end, sv_len, qual, filter,
genotype`` -- and normalizes every input file into it in one place. INV,
BND, CNV, TRA and small (<50 bp) variants are out of scope and are counted,
not silently dropped.

Every file is iterated with plain ``for rec in vcf`` (never ``vcf.fetch()``,
which needs an index and reorders records by contig -- the internal id
``f"{caller}:{ordinal}"`` must follow plain file order to stay stable).

Some pysam 0.24 behaviours this module relies on and works around:

* ``'END' in rec.info`` is always ``False`` -- INFO/END is only visible via
  ``rec.stop``, and htslib >= 1.20 derives ``rec.stop`` from SVLEN for
  symbolic non-INS ALTs even when END is absent.
* ``rec.info['SVLEN']`` is a tuple when the header declares ``Number=.``
  (common in caller VCFs); :func:`parse_svlen` unwraps it.
* ``rec.samples[s].get('GT')`` must be used instead of ``rec.samples[s]['GT']``
  -- raw lumpy records have no GT FORMAT field at all and ``[...]`` raises.
* ``sample in rec.samples`` must be checked before ``rec.samples[sample]`` at
  all -- manta (and presumably other callers) can emit records with a
  completely empty FORMAT/sample column (not merely a missing GT value);
  ``rec.samples`` is then a zero-length container for that one record, and
  indexing it by name raises ``IndexError: Invalid sample index`` rather
  than behaving like a dict that simply lacks the key.
* ``rec.alts`` is ``None`` for sites with no ALT (rare, but seen in gVCF-ish
  output); every ALT access here guards for that.
* htslib does not round-trip records faithfully (INFO/END is re-derived
  from SVLEN for symbolic ALTs, GT phasing is lost when records are rebuilt),
  so :func:`write_dq_tagged_vcf` works on the raw text and never re-writes a
  parsed record.
"""
from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd
import pysam

# Chromosomes dicast processes by default, and the SV types its model
# supports. Moved here from cli.py so the CLI, `dicast check` and this
# module share one definition.
DEFAULT_CHROMS = [
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
    'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
    'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX',
]
SUPPORTED_SV_TYPES = ['DEL', 'DUP', 'INS']

# Column set (and order) of the dataframe read_caller_vcf returns: the
# columns the rest of the pipeline (prepare.py, collect_*.py, model.py,
# merge.py) has always consumed, plus `vcf_id` right after `id`.
RAW_COLUMNS = [
    'id', 'vcf_id', 'sv_type', 'chrom', 'start', 'chrom_2', 'end', 'sv_len',
    'filter', 'qual', 'genotype', 'cohort', 'sample', 'technology', 'caller',
    'reference',
]

# Reasons a record can be dropped, and the per-file FileReport.dropped keys.
DROP_REASONS = (
    'contig', 'no_svtype', 'INV', 'BND', 'CNV', 'other_svtype',
    'small_variant', 'missing_end', 'bnd_mate', 'absent_in_sample',
)
# Informational per-file counters: these do not drop the record.
NOTE_KEYS = ('ins_len_unknown', 'no_gt')

# Sequence-resolved REF/ALT pairs shorter than this are not structural
# variants by dicast's contract and are dropped as 'small_variant'.
MIN_SV_LEN = 50

# Column dtypes to force whenever a dicast TSV (raw / annot / scores) is read
# back: `vcf_id` holds callers' original record IDs, which are purely numeric
# for lumpy/svtyper/smoove/dysgu and would otherwise come back as int64, or
# as float64 ('249.0') as soon as one row has none (a '.' ID, a rescue row).
TSV_DTYPES = {'sample': str, 'cohort_samples': str, 'vcf_id': str}

# chr1 length in the hg38/GRCh38 primary assembly; used only to warn when an
# FAI does not look like hg38.
HG38_CHR1_LENGTH = 248956422

_BND_BRACKET_RE = re.compile(r'[\[\]]([^:\[\]]+):(\d+)[\[\]]')
_INFO_END_RE = re.compile(r'(?:^|;)END=(\d+)(?:;|$)')


class VcfInputError(Exception):
    """Raised for problems that make an input VCF (or a whole run) unusable:
    an FAI/`--chrom` combination with no overlap, or (by callers of this
    module, not this module itself) zero usable records across all files."""


@dataclass
class NormalizedRecord:
    """One VCF record normalized into dicast's contract. Field names and
    semantics are the RAW_COLUMNS downstream code (prepare.py, model.py,
    merge.py) consumes."""

    id: str
    vcf_id: Optional[str]
    sv_type: str
    chrom: str
    start: int
    chrom_2: float
    end: int
    sv_len: float
    filter: str
    qual: object
    genotype: tuple


@dataclass
class FileReport:
    """Per-file summary of what :func:`read_caller_vcf` understood: how many
    records were read, kept (by SV type) and dropped (by reason), plus
    contig-mapping and sample-column bookkeeping. `dicast check` renders
    this via :func:`report_table` / :func:`reports_to_dataframe`."""

    path: str
    caller: str
    records_read: int = 0
    kept: Dict[str, int] = field(
        default_factory=lambda: {t: 0 for t in SUPPORTED_SV_TYPES})
    dropped: Dict[str, int] = field(
        default_factory=lambda: {r: 0 for r in DROP_REASONS})
    notes: Dict[str, int] = field(
        default_factory=lambda: {n: 0 for n in NOTE_KEYS})
    contig_action: str = 'as-is'
    contig_mapped: int = 0
    contig_total: int = 0
    sample_used: Optional[str] = None
    warnings: List[str] = field(default_factory=list)
    error: Optional[str] = None

    @property
    def usable(self) -> int:
        """Total records that survived normalization (sum of `kept`)."""
        return sum(self.kept.values())

    @property
    def status(self) -> str:
        """One of 'ERROR' (file could not be read at all), 'ZERO USABLE'
        (read fine but nothing survived) or 'OK'."""
        if self.error:
            return 'ERROR'
        if self.usable == 0:
            return 'ZERO USABLE'
        return 'OK'


def read_fai_contigs(fai_path: str) -> list:
    """ Reads an .fai index into a list of (contig_name, length) in file order. """

    contigs = []
    with open(fai_path) as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            contigs.append((fields[0], int(fields[1])))
    return contigs


def canonical_chroms(requested: List[str], fai_contigs: List[Tuple[str, int]]) -> Tuple[List[str], List[str]]:
    """Restricts the requested chromosome list to those actually present in
    the FAI, keeping the requested order.

    Args:
        requested (list): chromosome names to run on (e.g. `--chrom` or
            DEFAULT_CHROMS).
        fai_contigs (list): (name, length) pairs from read_fai_contigs.

    Returns:
        tuple: (kept, warnings). `kept` is `requested` intersected with the
        FAI's contig names, in `requested`'s order.

    Raises:
        VcfInputError: if none of `requested` is present in the FAI (a run
        with no usable chromosomes cannot proceed).
    """
    fai_names = {name for name, _ in fai_contigs}
    kept = [c for c in requested if c in fai_names]
    missing = [c for c in requested if c not in fai_names]

    warnings = []
    if missing:
        warnings.append(
            f"Requested chromosome(s) not found in the FAI, skipped: {', '.join(missing)}"
        )
    if not kept:
        raise VcfInputError(
            f"None of the requested chromosomes ({', '.join(requested)}) are present in "
            f"the FAI (its first contigs are: {', '.join(name for name, _ in fai_contigs[:5])})"
        )
    return kept, warnings


def fai_warnings(fai_contigs: List[Tuple[str, int]]) -> List[str]:
    """Sanity-checks an FAI against dicast's hg38, chr-named assumptions.

    Args:
        fai_contigs (list): (name, length) pairs from read_fai_contigs.

    Returns:
        list: human-readable warnings (empty if nothing looks off). This
        never raises -- non-hg38 input is allowed, just flagged loudly, per
        the "warn and continue" decision for non-hg38 FAIs.
    """
    warnings = []
    fai_lengths = dict(fai_contigs)

    missing_default = [c for c in DEFAULT_CHROMS if c not in fai_lengths]
    if missing_default:
        warnings.append(
            "FAI is missing these chr-named chromosomes: "
            f"{', '.join(missing_default)}. dicast's BAM/FAI handling expects "
            "chr1..chr22 and chrX (collect_illumina.py parses int(chrom[3:]))."
        )

    chr1_length = fai_lengths.get('chr1')
    if chr1_length is not None and chr1_length != HG38_CHR1_LENGTH:
        warnings.append(
            f"FAI chr1 length is {chr1_length}, not {HG38_CHR1_LENGTH} -- "
            "this does not look like hg38/GRCh38."
        )
    return warnings


def _map_one_contig(raw_chrom: str, canonical: set) -> Optional[str]:
    """Tries identity, then a 'chr' prefix, then stripping one, against the
    canonical chromosome set. Returns None if none of the three match."""
    if raw_chrom in canonical:
        return raw_chrom
    prefixed = 'chr' + raw_chrom
    if prefixed in canonical:
        return prefixed
    if raw_chrom.startswith('chr') and raw_chrom[3:] in canonical:
        return raw_chrom[3:]
    return None


def _classify_mapping(raw_chrom: str, mapped_chrom: str) -> str:
    """Labels how `raw_chrom` was mapped, for FileReport.contig_action."""
    if raw_chrom == mapped_chrom:
        return 'as-is'
    if mapped_chrom == 'chr' + raw_chrom:
        return '+chr'
    return '-chr'


def _first(value):
    """Unwraps a pysam INFO value that may come back as a tuple (headers
    declaring Number=. do this even for effectively-scalar fields)."""
    if isinstance(value, tuple):
        return value[0] if value else None
    return value


def _info_get(rec, key: str):
    """Safe form of ``rec.info.get(key)`` for this pysam build.

    Unlike a plain dict, ``VariantRecordInfo.get()`` only tolerates a
    *declared* key that happens to be absent on this record -- looking up a
    key the header never declares at all raises ``ValueError: Invalid
    header`` instead of returning None. Real caller VCFs routinely omit the
    declaration for keys this module probes speculatively (SVTYPE,
    SIMPLE_TYPE, SVLEN are each declared by some callers and not others), so
    every INFO lookup in this module goes through this guard instead of
    ``rec.info.get`` directly.
    """
    if key not in rec.header.info:
        return None
    return rec.info.get(key)


def bracket_mate(alt: Optional[str]) -> Optional[Tuple[str, int]]:
    """Extracts the mate (chrom, pos) from a breakend ALT's bracket
    notation, e.g. ``'N[chr5:9000['`` -> ``('chr5', 9000)``. Covers all four
    breakend orientations (`t[p[`, `t]p]`, `[p[t`, `]p]t`).

    Args:
        alt (str | None): one ALT allele string.

    Returns:
        tuple | None: (mate_chrom, mate_pos), or None if `alt` carries no
        bracket-mate notation (i.e. is not a breakend ALT).
    """
    if not alt:
        return None
    match = _BND_BRACKET_RE.search(alt)
    if not match:
        return None
    return match.group(1), int(match.group(2))


def derive_sv_type(rec) -> Optional[str]:
    """Determines a record's raw SV type string, following dicast's
    precedence chain, with the subtype stripped at the first colon
    (`'DUP:TANDEM'` -> `'DUP'`).

    Precedence: INFO/SVTYPE, unless it is 'BND' -- in that case INFO/
    SIMPLE_TYPE (gridss's `simple-event-annotation.R`) is tried first, so a
    gridss-annotated file is not left looking breakend-only just because
    every record is technically still a BND line -- then INFO/SIMPLE_TYPE,
    then 'BND' as SVTYPE's literal value, then a symbolic ALT
    (`<DUP:TANDEM>`, `<DEL:ME:ALU>`, `<CN0>`, ...), then bracket-mate ALT
    notation (-> 'BND').

    Args:
        rec (pysam.VariantRecord): the record to classify.

    Returns:
        str | None: the raw type, or None when none of SVTYPE/SIMPLE_TYPE/
        symbolic-ALT/bracket-ALT yields anything -- the caller then tries
        the sequence-resolved REF/ALT length-difference fallback, which
        also needs the small-variant size threshold and so is not folded in
        here.
    """
    svtype = _first(_info_get(rec, 'SVTYPE'))
    simple_type = _first(_info_get(rec, 'SIMPLE_TYPE'))

    alt = rec.alts[0] if rec.alts else None
    symbolic = alt[1:-1].split(':')[0] if alt and alt.startswith('<') and alt.endswith('>') else None

    if svtype == 'CNV' and symbolic in ('DEL', 'DUP'):
        # DRAGEN-style CNV output: SVTYPE=CNV but the ALT says which.
        return symbolic
    if svtype and svtype != 'BND':
        return str(svtype).split(':')[0]
    if simple_type:
        return str(simple_type).split(':')[0]
    if svtype == 'BND':
        return 'BND'

    if symbolic is not None:
        return symbolic
    if alt and bracket_mate(alt) is not None:
        return 'BND'
    return None


def _seq_resolved_sv_type(rec) -> Tuple[Optional[str], Optional[str]]:
    """Classifies a record with no SVTYPE/SIMPLE_TYPE/symbolic/bracket ALT
    by REF/ALT sequence length difference (>=50 bp: DEL if REF is longer,
    INS otherwise; smaller is dropped as 'small_variant').

    Returns:
        tuple: (sv_type, drop_reason), exactly one not None. 'no_svtype' is
        used when there is no usable ALT sequence to compare at all (missing
        ALT, an unrecognized symbolic ALT, or the `*` overlapping-deletion
        allele).
    """
    alt = rec.alts[0] if rec.alts else None
    if not alt or alt.startswith('<') or alt == '*':
        return None, 'no_svtype'
    ref = rec.ref or ''
    diff = len(alt) - len(ref)
    if abs(diff) < MIN_SV_LEN:
        return None, 'small_variant'
    return ('DEL' if diff < 0 else 'INS'), None


def info_end(rec) -> Optional[int]:
    """Returns INFO/END as written in the file, or None if absent.

    pysam hides END behind ``rec.stop``, and htslib >= 1.20 sets ``rec.stop``
    to max(END, POS + |SVLEN|) for symbolic ALTs -- so a caller whose SVLEN
    is one larger than END - POS (cnvnator2VCF writes POS=1, END=10000,
    SVLEN=-10000) would get an end one base past its own END. END is the
    authoritative field per the VCF spec, so it is read from the record's
    text instead.

    Args:
        rec (pysam.VariantRecord): the record to read END from.

    Returns:
        int | None: the INFO/END value, or None when the record has none.
    """
    match = _INFO_END_RE.search(str(rec).split('\t', 8)[7])
    return int(match.group(1)) if match else None


def parse_svlen(rec) -> Optional[int]:
    """Returns INFO/SVLEN as a plain int, or None if absent.

    Args:
        rec (pysam.VariantRecord): the record to read SVLEN from.

    Returns:
        int | None: pysam represents SVLEN as a tuple when the header
        declares Number=. (common in caller VCFs); this takes the first
        element.
    """
    svlen = _first(_info_get(rec, 'SVLEN'))
    if svlen is None:
        return None
    return int(svlen)


def normalize_record(rec, ordinal: int, caller: str, chrom: str, sample: Optional[str]
                      ) -> Tuple[Optional[NormalizedRecord], Optional[str], List[str]]:
    """Normalizes one already contig-mapped VCF record into dicast's
    contract (id, sv_type, chrom, start, end, sv_len, filter, qual,
    genotype). Contig mapping and the "keep only if GT has an alt allele"
    multi-sample policy are file-level decisions and are applied by
    :func:`read_caller_vcf` around this call, not here.

    Args:
        rec (pysam.VariantRecord): the record, with its original (unmapped)
            rec.chrom.
        ordinal (int): 0-based index of this record in `for rec in vcf`
            order; becomes the internal id `f"{caller}:{ordinal}"`.
        caller (str): caller label, used only to build the internal id.
        chrom (str): the canonical chrom name rec.chrom was mapped to.
        sample (str | None): the resolved sample column to read GT from, or
            None when the file has no samples (e.g. a population catalog).

    Returns:
        tuple: (NormalizedRecord | None, drop_reason | None, notes).
        Exactly one of the first two is not None. `notes` are informational
        counters ('ins_len_unknown', 'no_gt') that apply alongside a kept
        record, never a drop.
    """
    notes: List[str] = []
    raw_type = derive_sv_type(rec)
    svtype_is_bnd = _first(_info_get(rec, 'SVTYPE')) == 'BND'
    # True only when SVTYPE said BND but SIMPLE_TYPE overrode it -- these
    # are still physically BND lines (bracket-mate ALT), so they need the
    # lower-breakend/mate handling below instead of the plain END/SVLEN path.
    uses_simple_type = svtype_is_bnd and raw_type not in (None, 'BND')

    if raw_type is None:
        raw_type, drop_reason = _seq_resolved_sv_type(rec)
        if raw_type is None:
            return None, drop_reason, notes

    if raw_type == 'INV':
        return None, 'INV', notes
    if raw_type == 'CNV':
        return None, 'CNV', notes
    if raw_type == 'BND' and not uses_simple_type:
        return None, 'BND', notes
    if raw_type not in SUPPORTED_SV_TYPES:
        return None, 'other_svtype', notes

    pos = rec.pos

    svlen = parse_svlen(rec)

    if uses_simple_type:
        alt = rec.alts[0] if rec.alts else None
        mate = bracket_mate(alt) if alt else None
        # Keep only the lower breakend (this record's pos before its
        # mate's, on the same contig); the other end of the pair is the
        # same event seen from the other side and would double-count it.
        if mate is None or mate[0] != rec.chrom or pos >= mate[1]:
            return None, 'bnd_mate', notes
        if raw_type == 'INS':
            end = pos + 1
            # simple-event-annotation.R's SVLEN for an insertion is the
            # inserted length; a zero-length "insertion" is a breakpoint
            # with nothing inserted and nothing deleted, i.e. not an SV.
            sv_len = abs(svlen) if svlen is not None else float('nan')
            if svlen is None:
                notes.append('ins_len_unknown')
            elif sv_len == 0:
                return None, 'small_variant', notes
        elif svlen is not None:
            # simple-event-annotation.R writes SVLEN as the event length;
            # for a deletion breakend p[q[ the deleted bases are p+1..q-1,
            # so this equals the END a symbolic record would carry.
            end = pos + abs(svlen)
            sv_len = end - pos
        else:
            end = mate[1] - 1
            sv_len = end - pos
    elif raw_type in ('DEL', 'DUP'):
        # END is authoritative; rec.stop also covers a long REF allele
        # (sequence-resolved deletion) and htslib's SVLEN-derived stop.
        end = info_end(rec) if svlen is not None else None
        if end is None and rec.stop is not None and rec.stop > pos:
            end = rec.stop
        if end is None:
            if svlen is None:
                return None, 'missing_end', notes
            end = pos + abs(svlen)
        sv_len = end - pos
    else:  # INS
        end = pos + 1
        if svlen is not None:
            sv_len = abs(svlen)
        elif rec.alts and rec.ref and not rec.alts[0].startswith('<'):
            sv_len = len(rec.alts[0]) - len(rec.ref)
        else:
            sv_len = float('nan')
            notes.append('ins_len_unknown')

    vcf_id = rec.id if rec.id not in (None, '.') else None

    genotype = (None, None)
    if sample is not None:
        # A record can have a completely empty FORMAT/sample column (seen in
        # real manta output): rec.samples[sample] then raises IndexError
        # rather than behaving like a dict missing the key, so membership is
        # checked first instead of relying on .get() to catch it.
        gt = rec.samples[sample].get('GT') if sample in rec.samples else None
        if gt is None:
            notes.append('no_gt')
        else:
            genotype = tuple(gt)

    normalized = NormalizedRecord(
        id=f'{caller}:{ordinal}',
        vcf_id=vcf_id,
        sv_type=raw_type,
        chrom=chrom,
        start=pos,
        chrom_2=float('nan'),
        end=end,
        sv_len=sv_len,
        filter=', '.join(rec.filter.keys()),
        qual=rec.qual,
        genotype=genotype,
    )
    return normalized, None, notes


def _rows_to_dataframe(rows: List[NormalizedRecord], cohort: str, sample: str,
                        technology: str, caller: str, reference: str) -> pd.DataFrame:
    """Assembles the kept NormalizedRecords plus the constant, args-derived
    columns into a dataframe with exactly RAW_COLUMNS."""
    if not rows:
        return pd.DataFrame(columns=RAW_COLUMNS)

    df = pd.DataFrame({
        'id': [r.id for r in rows],
        'vcf_id': [r.vcf_id for r in rows],
        'sv_type': [r.sv_type for r in rows],
        'chrom': [r.chrom for r in rows],
        'start': [r.start for r in rows],
        'chrom_2': [r.chrom_2 for r in rows],
        'end': [r.end for r in rows],
        'sv_len': [r.sv_len for r in rows],
        'filter': [r.filter for r in rows],
        'qual': [r.qual for r in rows],
        'genotype': [r.genotype for r in rows],
    })
    df['cohort'] = cohort
    df['sample'] = sample
    df['technology'] = technology
    df['caller'] = caller
    df['reference'] = reference
    return df[RAW_COLUMNS]


def _resolve_sample(vcf_samples: List[str], requested_sample: str, report: FileReport
                     ) -> Tuple[Optional[str], bool]:
    """Implements dicast's sample-column policy: 0 samples (a catalog VCF
    such as PAV) -> no genotype column at all; 1 sample -> use it whatever
    its name, warning if it differs from --sample; >1 samples -> --sample
    must name one of them, else this is a per-file error.

    Returns:
        tuple: (sample_used, filter_by_gt). `filter_by_gt` is True only in
        the >1-sample case, where records absent in the chosen sample
        (`GT` with no allele > 0) are dropped by the caller.
    """
    if len(vcf_samples) == 0:
        return None, False
    if len(vcf_samples) == 1:
        sample_used = vcf_samples[0]
        if requested_sample is not None and sample_used != requested_sample:
            report.warnings.append(
                f"using this file's only sample column '{sample_used}' "
                f"(differs from --sample '{requested_sample}')"
            )
        return sample_used, False
    if requested_sample in vcf_samples:
        return requested_sample, True
    if requested_sample is None:
        report.error = (
            f"this file has {len(vcf_samples)} sample columns "
            f"({', '.join(vcf_samples)}); pass --sample to pick one"
        )
    else:
        report.error = (
            f"--sample '{requested_sample}' is not among this file's sample "
            f"columns: {', '.join(vcf_samples)}"
        )
    return None, False


def read_caller_vcf(path: str, caller: str, sample: str, canonical: List[str],
                     cohort: str, reference: str, technology: str
                     ) -> Tuple[pd.DataFrame, FileReport]:
    """Reads one caller's VCF and normalizes every record into dicast's
    contract in a single pass. This is the only place in dicast that reads a
    VCF for variant calls (see write_dq_tagged_vcf for the separate DQ-
    tagging output path).

    Args:
        path (str): VCF path (plain, bgzipped, or bgzipped+tabixed).
        caller (str): caller label; used to build internal ids
            (f"{caller}:{ordinal}") and stored in the 'caller' column.
        sample (str): the --sample the run is for; only used to pick a
            column when the file has more than one sample.
        canonical (list): canonical chromosome names contig names in this
            file must map onto (membership only, order irrelevant here).
        cohort, reference, technology (str): passed straight through to the
            like-named output columns.

    Returns:
        tuple: (DataFrame with exactly RAW_COLUMNS, FileReport). The
        DataFrame is empty (but has RAW_COLUMNS) when the file could not be
        opened, no sample column matched --sample, or no contig in the file
        could be mapped to `canonical` -- in every such case `report.error`
        explains why, rather than this function raising.
    """
    report = FileReport(path=path, caller=caller)
    canonical_set = set(canonical)
    contig_map: Dict[str, Optional[str]] = {}
    rows: List[NormalizedRecord] = []

    save = pysam.set_verbosity(0)
    try:
        vcf = pysam.VariantFile(path)
    except (OSError, ValueError) as exc:
        pysam.set_verbosity(save)
        report.error = f'could not open VCF: {exc}'
        return pd.DataFrame(columns=RAW_COLUMNS), report

    sample_used, filter_by_gt = _resolve_sample(list(vcf.header.samples), sample, report)
    if report.error:
        pysam.set_verbosity(save)
        return pd.DataFrame(columns=RAW_COLUMNS), report
    report.sample_used = sample_used

    distinct_raw_contigs = set()
    mapped_raw_contigs = set()
    mapping_kind_counts: Dict[str, int] = {}

    for ordinal, rec in enumerate(vcf):
        report.records_read += 1
        raw_chrom = rec.chrom
        distinct_raw_contigs.add(raw_chrom)

        if raw_chrom in contig_map:
            chrom = contig_map[raw_chrom]
        else:
            chrom = _map_one_contig(raw_chrom, canonical_set)
            contig_map[raw_chrom] = chrom
            if chrom is not None:
                mapped_raw_contigs.add(raw_chrom)
                kind = _classify_mapping(raw_chrom, chrom)
                mapping_kind_counts[kind] = mapping_kind_counts.get(kind, 0) + 1

        if chrom is None:
            report.dropped['contig'] += 1
            continue

        normalized, drop_reason, notes = normalize_record(rec, ordinal, caller, chrom, sample_used)
        for note in notes:
            report.notes[note] += 1

        if normalized is None:
            report.dropped[drop_reason] += 1
            continue

        if filter_by_gt and not any(a is not None and a > 0 for a in normalized.genotype):
            report.dropped['absent_in_sample'] += 1
            continue

        report.kept[normalized.sv_type] += 1
        rows.append(normalized)

    vcf.close()
    pysam.set_verbosity(save)

    report.contig_total = len(distinct_raw_contigs)
    report.contig_mapped = len(mapped_raw_contigs)
    if mapping_kind_counts:
        report.contig_action = max(mapping_kind_counts, key=mapping_kind_counts.get)
    if distinct_raw_contigs and not mapped_raw_contigs:
        report.error = (
            "none of this file's contigs could be mapped to the canonical "
            f"chromosomes (found: {', '.join(sorted(distinct_raw_contigs))})"
        )
        return pd.DataFrame(columns=RAW_COLUMNS), report

    df = _rows_to_dataframe(rows, cohort, sample, technology, caller, reference)
    return df, report


def check_files(vcfs: List[Tuple[str, str]], canonical: List[str], sample: str,
                 cohort: str = 'na', reference: str = 'na', technology: str = 'na'
                 ) -> List[FileReport]:
    """Runs :func:`read_caller_vcf` over every (caller, path) pair and
    returns just the reports -- the backend for `dicast check`, which never
    needs the dataframes themselves.

    Args:
        vcfs (list): (caller, path) pairs, as parsed from `--vcfs`.
        canonical (list): canonical chromosome names, from canonical_chroms.
        sample (str): the --sample to check files against.
        cohort, reference, technology (str): placeholders passed through to
            read_caller_vcf; `check` does not use the resulting columns.

    Returns:
        list: one FileReport per input file, in input order.
    """
    reports = []
    for caller, path in vcfs:
        _, report = read_caller_vcf(path, caller, sample, canonical, cohort, reference, technology)
        reports.append(report)
    return reports


def report_table(reports: List[FileReport]) -> str:
    """Renders FileReports as a human-readable table for `dicast check` and
    the log line prepare.py emits after reading every input file: one
    summary row per file (caller, status, records read/usable, path), then
    indented detail lines for kept/dropped counts, notes, contig handling,
    the sample column used, and warnings/errors.

    Args:
        reports (list): FileReports, in the order they should be printed.

    Returns:
        str: the rendered table (no trailing newline). A file with zero
        usable records and only BND drops gets an explicit hint, since that
        combination means the file is breakend-only (gridss without
        simple-event-annotation.R, or delly TRA/BND-only output) and dicast
        does not pair breakends itself.
    """
    caller_width = max([len('CALLER')] + [len(r.caller) for r in reports]) + 2
    header = f"{'CALLER':<{caller_width}}{'STATUS':<14}{'READ':>7}{'USABLE':>8}  PATH"
    lines = [header, '-' * len(header)]

    for r in reports:
        lines.append(f"{r.caller:<{caller_width}}{r.status:<14}{r.records_read:>7}{r.usable:>8}  {r.path}")
        if r.error:
            lines.append(f'    error: {r.error}')
            continue

        kept_str = ', '.join(f'{t}={r.kept.get(t, 0)}' for t in SUPPORTED_SV_TYPES)
        lines.append(f'    kept:    {kept_str}')

        dropped_nonzero = {k: v for k, v in r.dropped.items() if v}
        if dropped_nonzero:
            lines.append('    dropped: ' + ', '.join(f'{k}={v}' for k, v in dropped_nonzero.items()))

        notes_nonzero = {k: v for k, v in r.notes.items() if v}
        if notes_nonzero:
            lines.append('    notes:   ' + ', '.join(f'{k}={v}' for k, v in notes_nonzero.items()))

        lines.append(f'    contigs: {r.contig_action} ({r.contig_mapped}/{r.contig_total} mapped)')
        if r.sample_used is not None:
            lines.append(f'    sample:  {r.sample_used}')
        for w in r.warnings:
            lines.append(f'    warning: {w}')

        if r.usable == 0 and r.dropped.get('absent_in_sample', 0) > 0:
            lines.append(
                f"    hint: no record carries an alt allele for sample '{r.sample_used}'. "
                'If the file is not genotyped (GT ./. throughout), extract this '
                "sample's calls into a single-sample VCF instead."
            )
        if r.usable == 0 and r.dropped.get('BND', 0) > 0:
            lines.append(
                '    hint: this file looks breakend-only (every record dropped as BND). '
                'dicast does not pair breakends itself -- for gridss, annotate the raw VCF '
                'with scripts/gridss_simple_event_annotation.R first; for delly, use its '
                'DEL/DUP/INS calls instead of TRA/BND.'
            )
    return '\n'.join(lines)


def reports_to_dataframe(reports: List[FileReport]) -> pd.DataFrame:
    """Flattens FileReports into one fixed-column row per file, for `dicast
    check --out <tsv>`. The column set (every kept/dropped/note counter,
    always present even at 0) is stable across runs so downstream tooling
    can rely on it.

    Args:
        reports (list): FileReports, in the order they should appear.

    Returns:
        pandas.DataFrame: one row per file.
    """
    rows = []
    for r in reports:
        row = {
            'path': r.path,
            'caller': r.caller,
            'records_read': r.records_read,
            'usable': r.usable,
            'status': r.status,
            'sample_used': r.sample_used,
            'contig_action': r.contig_action,
            'contig_mapped': r.contig_mapped,
            'contig_total': r.contig_total,
            'error': r.error or '',
            'warnings': '; '.join(r.warnings),
        }
        for t in SUPPORTED_SV_TYPES:
            row[f'kept_{t}'] = r.kept.get(t, 0)
        for reason in DROP_REASONS:
            row[f'dropped_{reason}'] = r.dropped.get(reason, 0)
        for note in NOTE_KEYS:
            row[f'notes_{note}'] = r.notes.get(note, 0)
        rows.append(row)
    return pd.DataFrame(rows)


def _open_vcf_text(path: str):
    """Opens a plain or (b)gzipped VCF as a text stream (the DQ tagger works
    on raw lines, see write_dq_tagged_vcf)."""
    with open(path, 'rb') as probe:
        gzipped = probe.read(2) == b'\x1f\x8b'
    if gzipped:
        return gzip.open(path, 'rt')
    return open(path, 'r')


def write_dq_tagged_vcf(vcf_path: str, out_path: str, caller: str, scores: dict) -> int:
    """Writes a copy of `vcf_path` with an INFO/DQ tag added to every
    record, scored by the record's internal id (f"{caller}:{ordinal}",
    matching read_caller_vcf's id scheme over plain file order).

    This works on the raw text rather than through pysam on purpose: the
    output must be the user's VCF plus one tag, and htslib rewrites what it
    parses (INFO/END is dropped or re-derived from SVLEN for symbolic ALTs,
    GT phasing is lost when records are rebuilt, undeclared FORMAT fields
    with missing values fail to round-trip). Copying lines verbatim and
    appending `DQ=` to the INFO column keeps every caller dialect intact and
    never fails on a header that htslib would only tolerate leniently.

    Args:
        vcf_path (str): input VCF (plain or bgzipped).
        out_path (str): output VCF path (plain text).
        caller (str): caller label used to rebuild each record's internal id
            to look its score up in `scores`.
        scores (dict): internal id -> dicast_qual score. A record with no
            entry (e.g. dropped upstream by read_caller_vcf) gets -1, as the
            previous vcfpy-based tagging did.

    Returns:
        int: number of records written.
    """
    dq_header_line = '##INFO=<ID=DQ,Number=1,Type=String,Description="Dicast Quality Score">'
    ordinal = 0
    dq_declared = False
    with _open_vcf_text(vcf_path) as vcf_in, open(out_path, 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith('##'):
                if line.startswith('##INFO=<ID=DQ,'):
                    dq_declared = True
                vcf_out.write(line)
                continue
            if line.startswith('#'):
                if not dq_declared:
                    vcf_out.write(dq_header_line + '\n')
                vcf_out.write(line)
                continue
            if not line.strip():
                continue

            fields = line.rstrip('\r\n').split('\t')
            score = scores.get(f'{caller}:{ordinal}')
            dq = f'DQ={score if score is not None else -1}'
            # Replace a DQ left by an earlier dicast run rather than stacking.
            info = [f for f in fields[7].split(';') if f and f != '.' and not f.startswith('DQ=')]
            fields[7] = ';'.join(info + [dq])
            vcf_out.write('\t'.join(fields) + '\n')
            ordinal += 1
    return ordinal
