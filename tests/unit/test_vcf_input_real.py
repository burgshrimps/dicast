"""Tests for :mod:`dicast.vcf_input` against real SV-caller output.

``tests/unit/test_vcf_input.py`` proves the normalizer's rules against
synthetic, single-purpose fixtures. This module proves the same rules
against the real caller files in ``tests/data/callers/`` (see that
directory's ``README.md`` for full provenance) -- messy real headers, real
ID schemes, real sample-column dialects, all at once per file.

How CALLER_FIXTURES' expected numbers were derived
----------------------------------------------------
Every ``read``/``kept``/``dropped``/``notes``/``sample_used`` value below was
computed by an independent, from-scratch re-implementation of the plan's
normalization rules (contig mapping, the SVTYPE/SIMPLE_TYPE/symbolic-ALT/
bracket-ALT/sequence-resolved precedence chain, the END/SVLEN-per-type
rules, the 0/1/many sample-column policy) driven directly by pysam over each
file's raw records -- written and run *before* this module existed, without
importing ``dicast.vcf_input`` at all, so the numbers are not "run the
parser and copy its output". That standalone script is not part of the test
suite (it lives in scratch space); reproducing its logic here would just
recreate ``dicast/vcf_input.py`` a second time, so instead each entry below
documents in a comment what the file's own dialect implies about the count,
citing the SVTYPE/GT/END facts already recorded in
``tests/data/callers/README.md`` and the caller-inventory notes it summarizes.

Deriving these numbers independently is exactly what surfaced a real bug,
fixed as part of this task: ``manta_hg38_mixed.vcf.gz`` has 5 of 11 records
with a completely empty FORMAT/sample column (no FORMAT field at all, not
merely a missing GT value); ``rec.samples[sample_name]`` raises
``IndexError: Invalid sample index`` for such a record instead of behaving
like a dict that simply lacks the key, which crashed
``read_caller_vcf``. Fixed in ``dicast/vcf_input.py::normalize_record`` by
checking ``sample in rec.samples`` before indexing (see the fix's comment
there, and item 11 in the callers README's "Problems encountered").

Canonical chromosomes for every fixture are DEFAULT_CHROMS (chr1..chr22,
chrX) intersected with ``tests/data/hg38.fa.fai`` -- the same chr-prefixed,
hg38 FAI the shipped demo dataset uses. This is deliberate for the GRCh37/
b37-numeric fixtures (``delly_grch37_hg002_del``, ``manta_grch37_hg002``,
``cnvnator_grch37_hg002``, ``lumpy_raw_b37``, ``lumpy_svtyper_b37``,
``gridss_raw_grch37``): none of their contigs are declared with a ``chr``
prefix, so every kept record is only reachable through the ``'chr'+name``
mapping branch, exercising it against real data instead of only the
synthetic fixtures in ``test_vcf_input.py``.
"""
from __future__ import annotations

import re
import subprocess
import sys

import pandas as pd
import pytest

from tests.conftest import REPO_DIR
from dicast import vcf_input as vi

CALLERS_DIR = REPO_DIR / 'tests/data/callers'
FAI_PATH = REPO_DIR / 'tests/data/hg38.fa.fai'


def _canonical():
    fai_contigs = vi.read_fai_contigs(str(FAI_PATH))
    canonical, _warnings = vi.canonical_chroms(vi.DEFAULT_CHROMS, fai_contigs)
    return canonical


# ---------------------------------------------------------------------------
# CALLER_FIXTURES: one entry per file in tests/data/callers/.
#
# `sample` is the --sample this file is checked/run against (per the task:
# "single-sample files: any name works" -- picked as a deliberately
# unrelated 'test' string for every single-sample file below, to prove the
# "use the file's one sample column whatever its name is, just warn"
# policy; trio files instead pick one real sample name so --sample actually
# has to match one of several columns; sites-only files' sample value is
# unused. `kept`/`dropped`/`notes` list only non-zero entries -- every other
# SUPPORTED_SV_TYPES/DROP_REASONS/NOTE_KEYS key is implicitly 0.
# ---------------------------------------------------------------------------

CALLER_FIXTURES = [
    dict(
        file='delly_hg38_trio.vcf.gz', caller='delly', sample='NA12878_S1',
        # 9 records: SVTYPE DEL=2, INS=1, INV=3, DUP=2, TRA=1 (README). TRA is
        # not INV/CNV/BND/DEL/DUP/INS -> 'other_svtype' (1). INV (3) drops as
        # 'INV'. Of the 2 DUP records, trio genotyping means one is homref/
        # missing-alt in NA12878_S1 specifically (verified by hand: rec 6's
        # NA12878_S1 GT has no ALT allele) -> 'absent_in_sample' (1), leaving
        # DUP=1 kept. DEL=2/INS=1 keep whole (no per-sample filtering removes
        # them). kept total 2+1+1=4; dropped 3+1+1=5; read=9=4+5.
        sample_used='NA12878_S1', status='OK', read=9,
        kept={'DEL': 2, 'DUP': 1, 'INS': 1},
        dropped={'INV': 3, 'other_svtype': 1, 'absent_in_sample': 1},
        notes={},
    ),
    dict(
        file='delly_grch37_hg002_del.vcf.gz', caller='delly', sample='test',
        # 100 records, SVTYPE=DEL for all 100 (README); single sample -> used
        # regardless of --sample='test'. delly has no INFO/SVLEN in this
        # header at all, but every record's END/rec.stop is present
        # (end_from_stop_present=100/100 per the inventory), so every DEL
        # keeps via the rec.stop path -- 0 dropped.
        sample_used='HG002-NA24385-50x.70_percent.markdup.realigned',
        status='OK', read=100,
        kept={'DEL': 100, 'DUP': 0, 'INS': 0}, dropped={}, notes={},
    ),
    dict(
        file='manta_hg38_mixed.vcf.gz', caller='manta', sample='test',
        # 11 records: SVTYPE BND=2, DEL=3, INS=3, INV=1, DUP=2 (one
        # <DUP:TANDEM>, subtype stripped). DEL/INS/DUP (8) all keep; BND (2)
        # and INV (1) drop by name. Single sample -> used regardless of
        # --sample. 5 of 11 records have a completely empty FORMAT/sample
        # column (README item 11): all 5 are among the 8 kept DEL/INS/DUP
        # records (verified by hand against the record dump), so they
        # contribute 'no_gt' notes (genotype stays (None,None), not a drop)
        # rather than changing the kept/dropped counts. 8/11 SVLEN values are
        # pysam tuples; one kept INS record has neither SVLEN nor a
        # sequence-resolved ALT length to fall back on -> 'ins_len_unknown'.
        sample_used='/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/sv_benchmark/scripts/../data.chm/chm1.1.fq',
        status='OK', read=11,
        kept={'DEL': 3, 'DUP': 2, 'INS': 3}, dropped={'BND': 2, 'INV': 1},
        notes={'ins_len_unknown': 1, 'no_gt': 5},
    ),
    dict(
        file='manta_grch37_hg002.vcf.gz', caller='manta', sample='test',
        # 100 records: SVTYPE INS=37, DUP=11, DEL=45, BND=4, INV=3 (README).
        # DEL/DUP/INS (93) keep; BND (4) and INV (3) drop. Single sample ->
        # used regardless of --sample. 88/100 SVLEN values are pysam tuples;
        # 8 kept INS records have neither SVLEN nor a usable sequence-
        # resolved length -> 'ins_len_unknown' (verified against the
        # inventory's "88/100 SVLEN as tuple" note: 12 records have no SVLEN
        # at all, of which the DEL/DUP ones fall back to rec.stop fine and 8
        # are the affected INS).
        sample_used='HG002-NA24385-50x.70_percent', status='OK', read=100,
        kept={'DEL': 45, 'DUP': 11, 'INS': 37}, dropped={'BND': 4, 'INV': 3},
        notes={'ins_len_unknown': 8},
    ),
    dict(
        file='lumpy_raw_b37.vcf.gz', caller='lumpy', sample='test',
        # 100 records: SVTYPE DEL=96, INV=3, BND=1 (README). Sites-only (0
        # samples) -> sample_used=None, genotype always (None,None), no
        # per-sample filtering possible, and no 'no_gt' note either (that
        # note only fires when a sample column *was* being read from).
        # DEL (96) keeps; INV (3) and BND (1) drop by name.
        sample_used=None, status='OK', read=100,
        kept={'DEL': 96, 'DUP': 0, 'INS': 0}, dropped={'INV': 3, 'BND': 1},
        notes={},
    ),
    dict(
        file='lumpy_svtyper_b37.vcf.gz', caller='lumpy', sample='test',
        # Same 100-record window as lumpy_raw_b37 except the file's one BND
        # record falls just outside it (README) -> SVTYPE DEL=97, INV=3
        # only. Now genotyped (1 sample, NA12878) -> used regardless of
        # --sample. DEL (97) keeps; INV (3) drops.
        sample_used='NA12878', status='OK', read=100,
        kept={'DEL': 97, 'DUP': 0, 'INS': 0}, dropped={'INV': 3}, notes={},
    ),
    dict(
        file='lumpy_hg38_trio.vcf.gz', caller='lumpy', sample='NA12878_S1',
        # 12 records: SVTYPE DEL=2, BND=8, DUP=2 (README, "BND-heavy"). BND
        # (8) drops by name. Of the 2 DUP records, one is absent-in-sample
        # for NA12878_S1 specifically (verified by hand, same trio-
        # genotyping situation as delly_hg38_trio) -> 'absent_in_sample' (1),
        # DUP=1 kept. DEL=2 keeps whole. kept 2+1+0=3; dropped 8+1=9; read=12.
        sample_used='NA12878_S1', status='OK', read=12,
        kept={'DEL': 2, 'DUP': 1, 'INS': 0},
        dropped={'BND': 8, 'absent_in_sample': 1}, notes={},
    ),
    dict(
        file='smoove_hg38_header.vcf.gz', caller='smoove', sample='test',
        # 1 record, DEL, sites-only (0 samples) -- this fixture's whole point
        # is the 3366-##contig-line header, not the record count.
        sample_used=None, status='OK', read=1,
        kept={'DEL': 1, 'DUP': 0, 'INS': 0}, dropped={}, notes={},
    ),
    dict(
        file='gridss_raw_hg19chr.vcf.gz', caller='gridss', sample='test',
        # 8 records, SVTYPE=BND for all 8, no SIMPLE_TYPE declared at all
        # (raw gridss, pre-simple-event-annotation.R) -> every record drops
        # as 'BND', usable=0 -> ZERO USABLE with the breakend-only hint. This
        # is the concrete "gridss raw must yield ZERO USABLE" case the task
        # calls for. Sites-only (0 samples).
        sample_used=None, status='ZERO USABLE', read=8,
        kept={'DEL': 0, 'DUP': 0, 'INS': 0}, dropped={'BND': 8}, notes={},
    ),
    dict(
        file='gridss_raw_grch37.vcf.gz', caller='gridss', sample='test',
        # Same situation as gridss_raw_hg19chr, GRCh37 dialect: 60 records,
        # all SVTYPE=BND, no SIMPLE_TYPE -> ZERO USABLE. Sites-only.
        sample_used=None, status='ZERO USABLE', read=60,
        kept={'DEL': 0, 'DUP': 0, 'INS': 0}, dropped={'BND': 60}, notes={},
    ),
    dict(
        file='cnvnator_grch37_hg002.vcf.gz', caller='cnvnator', sample='test',
        # 100 records, SVTYPE DEL=48/DUP=52 (README), no INS. No ##contig
        # lines and no visible INFO/END via `'END' in rec.info` (htslib
        # quirk, see README item 10), but rec.stop is correct for all 100
        # (end_from_stop_present=100/100) so every record keeps.
        sample_used='cnvnator', status='OK', read=100,
        kept={'DEL': 48, 'DUP': 52, 'INS': 0}, dropped={}, notes={},
    ),
    # --- gameboycolor hg38 chr21 fixtures ---------------------------------
    dict(
        file='delly_hg38_chr21_hg002.vcf.gz', caller='delly', sample='test',
        # 27 records: SVTYPE DEL=15, INS=11, DUP=1 (README), no TRA/INV on
        # this input. INFO/END, INFO/SVLEN, FORMAT/GT all present on every
        # record (100% coverage per the inventory) -> nothing drops.
        sample_used='caller', status='OK', read=27,
        kept={'DEL': 15, 'DUP': 1, 'INS': 11}, dropped={}, notes={},
    ),
    dict(
        file='manta_hg38_chr21_hg002.vcf.gz', caller='manta', sample='test',
        # 28 records: SVTYPE DEL=17, INS=11 (README), no BND/DUP/INV on this
        # chr21 slice -> nothing drops. All 28 have INFO/SVLEN; 18 have
        # END/stop!=pos, meaning 10 do not (imprecise symbolic records) --
        # among the 11 INS, 1 has neither a usable SVLEN nor a sequence-
        # resolved length (verified by hand) -> 'ins_len_unknown'.
        sample_used='SAMPLE1', status='OK', read=28,
        kept={'DEL': 17, 'DUP': 0, 'INS': 11}, dropped={},
        notes={'ins_len_unknown': 1},
    ),
    dict(
        file='manta_candidatesv_hg38_chr21_hg002.vcf.gz', caller='manta',
        sample='test',
        # 100 trimmed records (raw candidateSV.vcf had 302: DEL=168, INS=134,
        # no BND/DUP/INV -- README). Sites-only (0 samples, pre-genotyping
        # candidate-stage VCF, no FORMAT column at all) -> sample_used=None.
        # All 100 kept, split DEL=58/INS=42 (verified against the actual
        # trimmed file: it's a positionally-selected subset of the raw
        # DEL/INS-only 302, so the ratio need not match the raw file's); 2
        # kept INS records lack a usable length -> 'ins_len_unknown'.
        sample_used=None, status='OK', read=100,
        kept={'DEL': 58, 'DUP': 0, 'INS': 42}, dropped={},
        notes={'ins_len_unknown': 2},
    ),
    dict(
        file='smoove_hg38_chr21_hg002.vcf.gz', caller='smoove', sample='test',
        # 4 records, symbolic <DEL> only (README), single sample HG002 with
        # non-missing GT on all 4 -> all keep, nothing drops.
        sample_used='HG002', status='OK', read=4,
        kept={'DEL': 4, 'DUP': 0, 'INS': 0}, dropped={}, notes={},
    ),
    dict(
        file='lumpy_raw_hg38_chr21_hg002.vcf.gz', caller='lumpy',
        sample='test',
        # 3 records, symbolic <DEL> only, single sample HG002 but GT is
        # always './.' (raw lumpyexpress, pre-svtyper genotyping -- README).
        # pysam represents './.' as the tuple (None, None), not a missing GT
        # field -- rec.samples[s].get('GT') returns that tuple rather than
        # None, so the 'no_gt' note (which only fires when GT is genuinely
        # absent, as in manta's empty-FORMAT records) does NOT fire here; all
        # 3 keep with genotype (None, None) and no note.
        sample_used='HG002', status='OK', read=3,
        kept={'DEL': 3, 'DUP': 0, 'INS': 0}, dropped={}, notes={},
    ),
    dict(
        file='gridss_hg38_chr21_hg002.vcf.gz', caller='gridss', sample='test',
        # 100 trimmed records, SVTYPE=BND for all of them (this fixture is
        # deliberately the raw, pre-simple-event-annotation.R gridss run --
        # README), no SIMPLE_TYPE declared -> every record drops as 'BND',
        # ZERO USABLE with the breakend-only hint. This is the second
        # required "gridss raw -> ZERO USABLE" case, and the reason this
        # file is deliberately left OUT of the multi-caller E2E `dicast
        # call` run in test_dicast.py (its companion
        # gridss_simple_hg38_chr21_hg002.vcf.gz is used there instead).
        sample_used='caller.bam', status='ZERO USABLE', read=100,
        kept={'DEL': 0, 'DUP': 0, 'INS': 0}, dropped={'BND': 100}, notes={},
    ),
    dict(
        file='gridss_simple_hg38_chr21_hg002.vcf.gz', caller='gridss',
        sample='test',
        # Same 912-record base as the raw fixture, trimmed to the same 100
        # records, now carrying INFO/SIMPLE_TYPE (README full-file counts:
        # unset=266, DEL=236, DUP=26, INS=370, INV=14, CTX=0). Every
        # SIMPLE_TYPE-bearing record is still a bracket-ALT BND line, so only
        # the lower breakend of each resolvable pair is kept
        # (normalize_record's uses_simple_type path: pos before the mate's
        # pos, same contig); its mate drops as 'bnd_mate'. Hand-counting this
        # 100-record window by SIMPLE_TYPE/pair-orientation gives 33 lower
        # breakends (DEL=13, DUP=2, INS=18), 32 of their upper mates
        # also inside the trim window (dropping as 'bnd_mate' -- the 33rd
        # pair's other end fell outside the trim window entirely, so it
        # contributes no drop here), 1 record with SIMPLE_TYPE=INV (drops as
        # 'INV'), and the remaining 34 records with SIMPLE_TYPE='.'
        # (unresolved single/unpaired breakends), which fall back to
        # SVTYPE='BND' (drops as 'BND'). Two of the 18 INS lower breakends
        # carry SVLEN=0 (gridss0fb_2o, gridss1fb_38o: nothing inserted,
        # nothing deleted) and drop as 'small_variant', leaving INS=16.
        sample_used='caller.bam', status='OK', read=100,
        kept={'DEL': 13, 'DUP': 2, 'INS': 16},
        dropped={'BND': 34, 'bnd_mate': 32, 'INV': 1, 'small_variant': 2}, notes={},
    ),
    dict(
        file='cnvnator_hg38_chr21_hg002.vcf.gz', caller='cnvnator',
        sample='test',
        # 38 records: SVTYPE DEL=21, DUP=17 (README), no INS. Same no-##contig/
        # no-visible-INFO/END-but-correct-rec.stop dialect as the GRCh37
        # cnvnator fixture -> nothing drops.
        sample_used='cnvnator', status='OK', read=38,
        kept={'DEL': 21, 'DUP': 17, 'INS': 0}, dropped={}, notes={},
    ),
    dict(
        file='dysgu_hg38_chr21_hg002.vcf.gz', caller='dysgu', sample='test',
        # 68 records: SVTYPE DEL=35, INS=33 (README), no DUP/INV/BND. INFO/
        # END (via rec.stop) present for 61/68 (DEL + precise INS); INFO/
        # SVLEN present for all 68, so every INS has a usable length and
        # nothing drops.
        sample_used='caller.dysgu_reads', status='OK', read=68,
        kept={'DEL': 35, 'DUP': 0, 'INS': 33}, dropped={}, notes={},
    ),
]

_FIXTURE_IDS = [f['file'] for f in CALLER_FIXTURES]


# ---------------------------------------------------------------------------
# Exact-match test through read_caller_vcf.
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.parametrize('fx', CALLER_FIXTURES, ids=_FIXTURE_IDS)
def test_read_caller_vcf_matches_hand_derived_counts(fx):
    canonical = _canonical()
    path = str(CALLERS_DIR / fx['file'])

    df, report = vi.read_caller_vcf(
        path, fx['caller'], fx['sample'], canonical, 'na', 'hg38', 'na')

    assert report.records_read == fx['read']
    assert report.sample_used == fx['sample_used']
    assert report.status == fx['status']
    assert report.error is None

    assert report.kept == fx['kept']
    expected_dropped = {r: 0 for r in vi.DROP_REASONS}
    expected_dropped.update(fx['dropped'])
    assert report.dropped == expected_dropped

    expected_notes = {n: 0 for n in vi.NOTE_KEYS}
    expected_notes.update(fx['notes'])
    assert report.notes == expected_notes

    # The dataframe agrees with the report: row counts per sv_type, and
    # every RAW_COLUMNS column present with nothing extra.
    assert list(df.columns) == vi.RAW_COLUMNS
    assert len(df) == report.usable
    for sv_type, count in fx['kept'].items():
        assert (df['sv_type'] == sv_type).sum() == count

    # Internal ids follow f"{caller}:{ordinal}" over plain file-order
    # iteration and are unique; vcf_id is either None or a real VCF ID.
    if report.usable:
        assert df['id'].is_unique
        assert (df['id'].str.fullmatch(rf'{re.escape(fx["caller"])}:\d+')).all()
        assert df['vcf_id'].apply(lambda v: v is None or isinstance(v, str)).all()


# ---------------------------------------------------------------------------
# Invariants that must hold for every real fixture regardless of dialect.
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.parametrize('fx', CALLER_FIXTURES, ids=_FIXTURE_IDS)
def test_read_caller_vcf_invariants_on_real_fixtures(fx):
    canonical = _canonical()
    path = str(CALLERS_DIR / fx['file'])

    df, report = vi.read_caller_vcf(
        path, fx['caller'], fx['sample'], canonical, 'na', 'hg38', 'na')

    # read == usable + every dropped reason, always.
    assert report.records_read == report.usable + sum(report.dropped.values())

    if df.empty:
        return

    # ids are unique and follow the internal ordinal scheme.
    assert df['id'].is_unique
    id_pattern = re.compile(rf'^{re.escape(fx["caller"])}:\d+$')
    assert df['id'].apply(lambda v: bool(id_pattern.match(v))).all()

    del_dup = df[df['sv_type'].isin(['DEL', 'DUP'])]
    assert (del_dup['sv_len'] == del_dup['end'] - del_dup['start']).all()

    # Every INS (plain or SIMPLE_TYPE-derived) is a point event with a
    # positive (or unknown) inserted length.
    ins = df[df['sv_type'] == 'INS']
    assert (ins['end'] == ins['start'] + 1).all()
    assert ((ins['sv_len'] > 0) | ins['sv_len'].isna()).all()


# ---------------------------------------------------------------------------
# Same expectations, through the `dicast check --out` subprocess.
# ---------------------------------------------------------------------------

@pytest.mark.integration
@pytest.mark.parametrize('fx', CALLER_FIXTURES, ids=_FIXTURE_IDS)
def test_check_subcommand_matches_hand_derived_counts(fx, tmp_path):
    path = str(CALLERS_DIR / fx['file'])
    out_tsv = tmp_path / 'report.tsv'

    cmd = [
        sys.executable, '-m', 'dicast', 'check',
        '--vcfs', f'{fx["caller"]}={path}',
        '--fai', str(FAI_PATH),
        '--sample', fx['sample'],
        '--out', str(out_tsv),
    ]
    result = subprocess.run(
        cmd, cwd=str(REPO_DIR), capture_output=True, text=True, timeout=60)

    expected_rc = 1 if fx['status'] == 'ZERO USABLE' else 0
    assert result.returncode == expected_rc, (
        f'dicast check exited {result.returncode}\n'
        f'stdout:\n{result.stdout}\nstderr:\n{result.stderr}')

    if fx['status'] == 'ZERO USABLE' and fx['dropped'].get('BND'):
        assert 'breakend-only' in result.stdout

    report = pd.read_csv(out_tsv, sep='\t', keep_default_na=False)
    assert len(report) == 1
    row = report.iloc[0]

    assert row['caller'] == fx['caller']
    assert int(row['records_read']) == fx['read']
    assert row['status'] == fx['status']
    sample_used = row['sample_used']
    expected_sample_used = fx['sample_used'] if fx['sample_used'] is not None else ''
    assert (sample_used if not pd.isna(sample_used) else '') == expected_sample_used

    for sv_type, count in fx['kept'].items():
        assert int(row[f'kept_{sv_type}']) == count
    for reason in vi.DROP_REASONS:
        assert int(row[f'dropped_{reason}']) == fx['dropped'].get(reason, 0)
