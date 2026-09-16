"""Unit tests for :mod:`dicast.vcf_input`.

This module is the only place dicast reads a VCF for variant calls, so the
tests exercise the whole precedence chain that turns an arbitrary record
into (or out of) dicast's contract: SV-type derivation (SVTYPE/SIMPLE_TYPE/
symbolic-ALT/bracket-ALT/sequence-resolved), END/sv_len rules per type,
contig mapping, the sample-column policy, id assignment, and the FileReport/
`dicast check` reporting surface, plus :func:`write_dq_tagged_vcf`.

Fixtures come from :mod:`tests.fixtures.synthetic_vcf`, which writes a
minimal VCF from plain record dicts (see its docstring for the schema).
EXPECTED values throughout are derived from what each test writes, never
from running the parser first.
"""
from __future__ import annotations

import pandas as pd
import pysam
import pytest

from dicast import vcf_input as vi
from tests.fixtures import synthetic_vcf as sv


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DEFAULT_CONTIGS = [('chr1', 100_000), ('chr2', 100_000)]
DEFAULT_CANONICAL = ['chr1', 'chr2']


def _read(tmp_path, records, *, samples=('S1',), sample='S1', canonical=None,
          contigs=None, extra_header=None, caller='testcaller', name='x.vcf'):
    """Writes `records` to a VCF and runs read_caller_vcf over it with
    sensible defaults (one sample 'S1', two declared canonical contigs)."""
    path = sv.write_vcf(
        tmp_path, name, records, samples=list(samples),
        contigs=DEFAULT_CONTIGS if contigs is None else contigs,
        extra_header=extra_header,
    )
    return vi.read_caller_vcf(
        path, caller, sample,
        DEFAULT_CANONICAL if canonical is None else canonical,
        'mycohort', 'hg38', 'illumina',
    )


def _open_records(tmp_path, records, *, samples=(), contigs=None, extra_header=None):
    """Writes `records` and returns the parsed pysam.VariantRecord list, for
    tests that exercise a single helper function directly rather than the
    full read_caller_vcf pipeline."""
    path = sv.write_vcf(
        tmp_path, 'x.vcf', records, samples=list(samples),
        contigs=DEFAULT_CONTIGS if contigs is None else contigs,
        extra_header=extra_header,
    )
    save = pysam.set_verbosity(0)
    vcf = pysam.VariantFile(path)
    recs = list(vcf)
    vcf.close()
    pysam.set_verbosity(save)
    return recs


# ---------------------------------------------------------------------------
# read_fai_contigs / canonical_chroms / fai_warnings
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_read_fai_contigs_parses_name_length_pairs_in_file_order(tmp_path):
    fai = sv.write_fai(tmp_path, 'ref.fa.fai', [('chr2', 200), ('chr1', 100)])
    assert vi.read_fai_contigs(fai) == [('chr2', 200), ('chr1', 100)]


@pytest.mark.unit
def test_canonical_chroms_keeps_requested_order_intersected_with_fai():
    fai_contigs = [('chr2', 100), ('chr1', 100), ('chr3', 100)]
    kept, warnings = vi.canonical_chroms(['chr1', 'chr2'], fai_contigs)
    assert kept == ['chr1', 'chr2']  # requested order, not FAI order
    assert warnings == []


@pytest.mark.unit
def test_canonical_chroms_warns_about_requested_chroms_missing_from_fai():
    fai_contigs = [('chr1', 100)]
    kept, warnings = vi.canonical_chroms(['chr1', 'chrX'], fai_contigs)
    assert kept == ['chr1']
    assert len(warnings) == 1
    assert 'chrX' in warnings[0]


@pytest.mark.unit
def test_canonical_chroms_raises_when_intersection_is_empty():
    fai_contigs = [('chr5', 100)]
    with pytest.raises(vi.VcfInputError):
        vi.canonical_chroms(['chr1', 'chr2'], fai_contigs)


@pytest.mark.unit
def test_fai_warnings_flags_missing_default_chromosomes():
    fai_contigs = [(c, 1_000_000) for c in vi.DEFAULT_CHROMS if c != 'chrX']
    warnings = vi.fai_warnings(fai_contigs)
    assert any('chrX' in w for w in warnings)


@pytest.mark.unit
def test_fai_warnings_flags_non_hg38_chr1_length():
    fai_contigs = [('chr1', 249_250_621)]  # GRCh37 chr1 length, not hg38's
    warnings = vi.fai_warnings(fai_contigs)
    assert any('hg38' in w for w in warnings)


@pytest.mark.unit
def test_fai_warnings_empty_for_clean_hg38_fai():
    fai_contigs = [(c, 1) for c in vi.DEFAULT_CHROMS]
    fai_contigs = [(name, vi.HG38_CHR1_LENGTH if name == 'chr1' else 1) for name, _ in fai_contigs]
    assert vi.fai_warnings(fai_contigs) == []


# ---------------------------------------------------------------------------
# bracket_mate
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize('alt, expected', [
    ('N[chr5:9000[', ('chr5', 9000)),
    ('N]chr5:9000]', ('chr5', 9000)),
    ('[chr5:9000[N', ('chr5', 9000)),
    (']chr5:9000]N', ('chr5', 9000)),
])
def test_bracket_mate_parses_all_four_breakend_orientations(alt, expected):
    assert vi.bracket_mate(alt) == expected


@pytest.mark.unit
def test_bracket_mate_none_for_symbolic_alt():
    assert vi.bracket_mate('<DEL>') is None


@pytest.mark.unit
def test_bracket_mate_none_for_missing_alt():
    assert vi.bracket_mate(None) is None


# ---------------------------------------------------------------------------
# parse_svlen
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_parse_svlen_unwraps_tuple_valued_svlen(tmp_path):
    # The base fixture header declares SVLEN Number=., so pysam hands it
    # back as a tuple even for one value.
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<DEL>', 'info': 'SVTYPE=DEL;SVLEN=-500'},
    ])
    assert isinstance(recs[0].info['SVLEN'], tuple)
    assert vi.parse_svlen(recs[0]) == -500


@pytest.mark.unit
def test_parse_svlen_none_when_absent(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<DEL>', 'info': 'SVTYPE=DEL'},
    ])
    assert vi.parse_svlen(recs[0]) is None


# ---------------------------------------------------------------------------
# derive_sv_type
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_derive_sv_type_reads_svtype_directly(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<DEL>', 'info': 'SVTYPE=DEL'},
    ])
    assert vi.derive_sv_type(recs[0]) == 'DEL'


@pytest.mark.unit
def test_derive_sv_type_simple_type_wins_when_svtype_is_bnd(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': 'N[chr1:200[', 'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL'},
    ])
    assert vi.derive_sv_type(recs[0]) == 'DEL'


@pytest.mark.unit
def test_derive_sv_type_bnd_literal_when_no_simple_type(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': 'N[chr1:200[', 'info': 'SVTYPE=BND'},
    ])
    assert vi.derive_sv_type(recs[0]) == 'BND'


@pytest.mark.unit
def test_derive_sv_type_symbolic_alt_subtype_stripped_at_first_colon(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<DUP:TANDEM>'},
        {'chrom': 'chr1', 'pos': 200, 'alt': '<DEL:ME:ALU>'},
    ])
    assert vi.derive_sv_type(recs[0]) == 'DUP'
    assert vi.derive_sv_type(recs[1]) == 'DEL'


@pytest.mark.unit
def test_derive_sv_type_bracket_alt_without_svtype_is_bnd(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': 'N[chr1:200['},
    ])
    assert vi.derive_sv_type(recs[0]) == 'BND'


@pytest.mark.unit
def test_derive_sv_type_none_for_plain_sequence_alt(tmp_path):
    recs = _open_records(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'ref': 'ACGTACGTAC', 'alt': 'A'},
    ])
    assert vi.derive_sv_type(recs[0]) is None


# ---------------------------------------------------------------------------
# SV-type sources, end-to-end through read_caller_vcf
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_svtype_info_field_is_used_directly(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ])
    assert report.kept == {'DEL': 1, 'DUP': 0, 'INS': 0}
    assert df.loc[0, 'sv_type'] == 'DEL'


@pytest.mark.unit
def test_sequence_resolved_del_end_comes_from_rec_stop(tmp_path):
    ref = 'A' * 61  # len 61 -> rec.stop = start + 61 = (pos-1)+61 = pos+60
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'SD', 'ref': ref, 'alt': 'A'},
    ])
    assert report.kept['DEL'] == 1
    row = df.set_index('vcf_id').loc['SD']
    assert row['sv_type'] == 'DEL'
    assert row['end'] == 1060
    assert row['sv_len'] == 60


@pytest.mark.unit
def test_sequence_resolved_ins_len_is_alt_minus_ref(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 2000, 'id': 'SI', 'ref': 'A', 'alt': 'A' + 'C' * 60},
    ])
    assert report.kept['INS'] == 1
    row = df.set_index('vcf_id').loc['SI']
    assert row['end'] == 2001
    assert row['sv_len'] == 60


@pytest.mark.unit
def test_sequence_resolved_below_threshold_dropped_as_small_variant(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 3000, 'id': 'SM', 'ref': 'ACGTACGTAC', 'alt': 'A'},  # |diff|=9
    ])
    assert report.usable == 0
    assert report.dropped['small_variant'] == 1


@pytest.mark.unit
def test_no_usable_alt_at_all_is_no_svtype(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 4000, 'id': 'STAR', 'alt': '*'},
    ])
    assert report.dropped['no_svtype'] == 1


@pytest.mark.unit
def test_inv_is_dropped(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<INV>', 'info': 'SVTYPE=INV;END=500'},
    ])
    assert report.dropped['INV'] == 1
    assert report.usable == 0


@pytest.mark.unit
def test_cnv_is_dropped(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<CN0>', 'info': 'SVTYPE=CNV'},
    ])
    assert report.dropped['CNV'] == 1


@pytest.mark.unit
def test_unrecognized_svtype_is_other_svtype(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<TRA>', 'info': 'SVTYPE=TRA'},
    ])
    assert report.dropped['other_svtype'] == 1


@pytest.mark.unit
def test_cnv_with_symbolic_del_dup_alt_takes_alt_type(tmp_path):
    # DRAGEN CNV output: SVTYPE=CNV but the ALT is an unambiguous <DEL>/<DUP>.
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'D', 'alt': '<DEL>', 'info': 'SVTYPE=CNV;END=600'},
        {'chrom': 'chr1', 'pos': 1000, 'id': 'U', 'alt': '<DUP>', 'info': 'SVTYPE=CNV;END=1600'},
        {'chrom': 'chr1', 'pos': 2000, 'id': 'C', 'alt': '<CNV>', 'info': 'SVTYPE=CNV;END=2600'},
    ])
    assert report.kept == {'DEL': 1, 'DUP': 1, 'INS': 0}
    assert report.dropped['CNV'] == 1
    assert df['sv_type'].tolist() == ['DEL', 'DUP']


@pytest.mark.unit
def test_end_is_taken_from_info_end_not_svlen(tmp_path):
    # cnvnator2VCF writes POS=1, END=10000, SVLEN=-10000; htslib's rec.stop
    # would be 10001 (max of END and POS+|SVLEN|). END is authoritative.
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1, 'id': 'CNVnator_del_1', 'alt': '<DEL>',
         'info': 'END=10000;SVTYPE=DEL;SVLEN=-10000'},
        {'chrom': 'chr1', 'pos': 20000, 'id': 'X', 'alt': '<DEL>',
         'info': 'SVTYPE=DEL;END=21000;SVLEN=-1500'},
    ])
    assert df['end'].tolist() == [10000, 21000]
    assert df['sv_len'].tolist() == [9999, 1000]


# ---------------------------------------------------------------------------
# SIMPLE_TYPE bracket pairs
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_simple_type_bracket_pair_keeps_lower_breakend_drops_mate(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'LO', 'alt': 'N[chr1:1500[',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL'},
        {'chrom': 'chr1', 'pos': 1500, 'id': 'HI', 'alt': ']chr1:1000]N',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL'},
    ])
    assert report.kept == {'DEL': 1, 'DUP': 0, 'INS': 0}
    assert report.dropped['bnd_mate'] == 1
    assert df['vcf_id'].tolist() == ['LO']
    row = df.iloc[0]
    # Deleted bases are 1001..1499, so END is the mate position minus one.
    assert row['end'] == 1499
    assert row['sv_len'] == 499


@pytest.mark.unit
def test_simple_type_lengths_come_from_svlen(tmp_path):
    # simple-event-annotation.R writes SVLEN; for DEL/DUP it defines end, for
    # INS it is the inserted length (the breakend distance is always 1 there).
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'D', 'alt': 'N[chr1:1500[',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL;SVLEN=-499'},
        {'chrom': 'chr1', 'pos': 2000, 'id': 'I', 'alt': 'NACGTACGT[chr1:2001[',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=INS;SVLEN=34'},
        {'chrom': 'chr1', 'pos': 3000, 'id': 'U', 'alt': 'N[chr1:3001[',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=INS'},
        {'chrom': 'chr1', 'pos': 4000, 'id': 'Z', 'alt': 'N[chr1:4001[',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=INS;SVLEN=0'},
    ])
    assert report.kept == {'DEL': 1, 'DUP': 0, 'INS': 2}
    assert report.dropped['small_variant'] == 1
    assert df.set_index('vcf_id').loc['D', ['end', 'sv_len']].tolist() == [1499, 499]
    assert df.set_index('vcf_id').loc['I', ['end', 'sv_len']].tolist() == [2001, 34]
    assert df.set_index('vcf_id').loc['U', 'end'] == 3001
    assert pd.isna(df.set_index('vcf_id').loc['U', 'sv_len'])
    assert report.notes['ins_len_unknown'] == 1


@pytest.mark.unit
def test_simple_type_bracket_pair_cross_contig_mate_is_dropped(tmp_path):
    # SIMPLE_TYPE=DEL implies an intra-chromosomal event; if the bracket
    # mate is on a different contig, dicast cannot trust the annotation and
    # drops it defensively rather than emitting a bogus end.
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'LO', 'alt': 'N[chr2:1500[',
         'info': 'SVTYPE=BND;SIMPLE_TYPE=DEL'},
    ])
    assert report.dropped['bnd_mate'] == 1
    assert report.usable == 0


# ---------------------------------------------------------------------------
# BND-only files
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_bnd_only_file_is_zero_usable_status(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': 'N[chr2:200[', 'info': 'SVTYPE=BND'},
        {'chrom': 'chr1', 'pos': 300, 'alt': 'N[chr2:400[', 'info': 'SVTYPE=BND'},
    ])
    assert report.status == 'ZERO USABLE'
    assert report.usable == 0
    assert report.dropped['BND'] == 2


@pytest.mark.unit
def test_report_table_hints_breakend_only_file(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': 'N[chr2:200[', 'info': 'SVTYPE=BND'},
    ])
    table = vi.report_table([report])
    assert 'breakend-only' in table
    assert 'gridss_simple_event_annotation.R' in table


@pytest.mark.unit
def test_report_table_no_hint_when_file_has_usable_records(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ])
    assert 'breakend-only' not in vi.report_table([report])


# ---------------------------------------------------------------------------
# END / SVLEN rules
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_del_end_from_info_end(tmp_path):
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=1500'},
    ])
    row = df.iloc[0]
    assert row['end'] == 1500
    assert row['sv_len'] == 500


@pytest.mark.unit
def test_del_end_from_svlen_when_end_absent(tmp_path):
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;SVLEN=-500'},
    ])
    row = df.iloc[0]
    assert row['end'] == 1500
    assert row['sv_len'] == 500


@pytest.mark.unit
def test_del_end_equal_to_pos_is_treated_as_missing(tmp_path):
    # A caller writing END=POS (no real span) must not silently produce
    # sv_len == 0; rec.stop == pos fails the 'rec.stop > pos' check and,
    # with no SVLEN to fall back to, the record is dropped.
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=1000'},
    ])
    assert report.dropped['missing_end'] == 1
    assert report.usable == 0


@pytest.mark.unit
def test_del_with_neither_end_nor_svlen_is_missing_end(tmp_path):
    _, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 1000, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL'},
    ])
    assert report.dropped['missing_end'] == 1


@pytest.mark.unit
def test_ins_end_is_always_pos_plus_one(tmp_path):
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 5000, 'id': 'A', 'alt': '<INS>', 'info': 'SVTYPE=INS;SVLEN=250'},
    ])
    assert df.iloc[0]['end'] == 5001
    assert df.iloc[0]['sv_len'] == 250


@pytest.mark.unit
def test_ins_len_unknown_noted_when_no_svlen_and_symbolic_alt(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 5000, 'id': 'A', 'alt': '<INS>', 'info': 'SVTYPE=INS'},
    ])
    assert report.notes['ins_len_unknown'] == 1
    assert pd.isna(df.iloc[0]['sv_len'])


# ---------------------------------------------------------------------------
# Contig mapping
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_contig_plus_chr_mapping(tmp_path):
    df, report = _read(
        tmp_path,
        [{'chrom': '1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'}],
        contigs=[],
    )
    assert report.contig_action == '+chr'
    assert report.contig_mapped == 1
    assert df.iloc[0]['chrom'] == 'chr1'


@pytest.mark.unit
def test_contig_minus_chr_mapping():
    # canonical names are themselves un-prefixed here, so a 'chr'-prefixed
    # record contig must be stripped to match.
    canonical = ['1', '2']
    assert vi._map_one_contig('chr1', set(canonical)) == '1'


@pytest.mark.unit
def test_contig_unmapped_is_dropped_and_counted(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chrUn_extra', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
        {'chrom': 'chr1', 'pos': 100, 'id': 'B', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ])
    assert report.dropped['contig'] == 1
    assert report.usable == 1
    assert df.iloc[0]['chrom'] == 'chr1'


@pytest.mark.unit
def test_contig_discovered_without_any_header_contig_lines(tmp_path):
    # cnvnator2VCF / raw lumpy output has no ##contig lines at all.
    df, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'}],
        contigs=[],
    )
    assert report.usable == 1
    assert report.contig_action == 'as-is'


@pytest.mark.unit
def test_all_contigs_unmapped_sets_file_error(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chrZZ', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ])
    assert report.error is not None
    assert 'chrZZ' in report.error
    assert report.status == 'ERROR'
    assert len(df) == 0


# ---------------------------------------------------------------------------
# ids
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_ordinal_ids_follow_plain_file_order_and_vcf_id_is_preserved(tmp_path):
    # File order deliberately not sorted by position, to prove ordinals
    # come from `for rec in vcf` order, not a re-sort.
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 5000, 'id': 'LATER', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=5500'},
        {'chrom': 'chr1', 'pos': 1000, 'id': 'EARLIER', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=1500'},
    ])
    assert df['id'].tolist() == ['testcaller:0', 'testcaller:1']
    assert df['vcf_id'].tolist() == ['LATER', 'EARLIER']


@pytest.mark.unit
def test_ordinal_counts_dropped_records_too(tmp_path):
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'DROPPED', 'alt': '<INV>', 'info': 'SVTYPE=INV;END=200'},
        {'chrom': 'chr1', 'pos': 1000, 'id': 'KEPT', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=1500'},
    ])
    assert df.iloc[0]['id'] == 'testcaller:1'  # index 1, not 0 -- ordinal 0 was the dropped INV


@pytest.mark.unit
def test_dot_id_becomes_none_vcf_id(tmp_path):
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': '.', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ])
    assert df.iloc[0]['vcf_id'] is None


@pytest.mark.unit
def test_ids_are_unique_within_a_file(tmp_path):
    records = [
        {'chrom': 'chr1', 'pos': 1000 + i * 1000, 'id': f'R{i}', 'alt': '<DEL>',
         'info': 'SVTYPE=DEL;END=' + str(1500 + i * 1000)}
        for i in range(5)
    ]
    df, _ = _read(tmp_path, records)
    assert df['id'].is_unique


# ---------------------------------------------------------------------------
# Sample-column policy
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_zero_samples_gets_missing_genotype(tmp_path):
    df, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'}],
        samples=(),
    )
    assert report.sample_used is None
    assert tuple(df.iloc[0]['genotype']) == (None, None)
    assert report.usable == 1  # not filtered out despite the missing genotype


@pytest.mark.unit
def test_single_sample_used_even_if_name_differs_from_requested(tmp_path):
    df, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
          'samples': {'CATALOG_SAMPLE': '0/1'}}],
        samples=('CATALOG_SAMPLE',), sample='requested_name',
    )
    assert report.sample_used == 'CATALOG_SAMPLE'
    assert any('CATALOG_SAMPLE' in w for w in report.warnings)
    assert tuple(df.iloc[0]['genotype']) == (0, 1)


@pytest.mark.unit
def test_multi_sample_without_match_is_a_file_error(tmp_path):
    df, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
          'samples': {'S1': '0/1', 'S2': '0/0'}}],
        samples=('S1', 'S2'), sample='not_present',
    )
    assert report.status == 'ERROR'
    assert 'S1' in report.error and 'S2' in report.error
    assert len(df) == 0


@pytest.mark.unit
def test_multi_sample_gt_filtering_keeps_partial_genotype_with_an_alt(tmp_path):
    df, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
          'samples': {'S1': './1', 'S2': '0/0'}}],
        samples=('S1', 'S2'), sample='S1',
    )
    assert report.usable == 1
    assert tuple(df.iloc[0]['genotype']) == (None, 1)


@pytest.mark.unit
def test_multi_sample_gt_filtering_drops_homref_and_counts_absent_in_sample(tmp_path):
    _, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
          'samples': {'S1': '0/0', 'S2': '0/1'}}],
        samples=('S1', 'S2'), sample='S1',
    )
    assert report.dropped['absent_in_sample'] == 1
    assert report.usable == 0


@pytest.mark.unit
def test_multi_sample_no_gt_at_all_is_noted_and_dropped(tmp_path):
    # A record whose FORMAT doesn't even include GT (raw lumpy-style):
    # rec.samples[s].get('GT') is None, not a KeyError.
    df, report = _read(
        tmp_path,
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
          'format': 'DP', 'samples': {'S1': '10', 'S2': '12'}}],
        samples=('S1', 'S2'), sample='S1',
        extra_header=['##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">'],
    )
    assert report.notes['no_gt'] == 1
    assert report.dropped['absent_in_sample'] == 1


@pytest.mark.unit
def test_sites_only_file_has_no_samples_and_keeps_records(tmp_path):
    path = sv.write_vcf(
        tmp_path, 'sites.vcf',
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'}],
        samples=[], contigs=DEFAULT_CONTIGS,
    )
    df, report = vi.read_caller_vcf(path, 'testcaller', 'anything', DEFAULT_CANONICAL, 'coh', 'hg38', 'illumina')
    assert report.sample_used is None
    assert report.usable == 1
    assert tuple(df.iloc[0]['genotype']) == (None, None)


@pytest.mark.unit
def test_sample_omitted_single_sample_no_warning_multi_sample_error(tmp_path):
    # `dicast check` makes --sample optional: a single-sample file is used
    # silently, a multi-sample file asks for the flag.
    _, single = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
         'samples': {'HG002': '0/1'}},
    ], samples=('HG002',), sample=None)
    assert single.usable == 1 and single.warnings == [] and single.sample_used == 'HG002'

    _, trio = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
         'samples': {'kid': '0/1', 'mom': '0/0', 'dad': '0/0'}},
    ], samples=('kid', 'mom', 'dad'), sample=None, name='trio.vcf')
    assert trio.status == 'ERROR'
    assert 'pass --sample' in trio.error and 'kid, mom, dad' in trio.error


# ---------------------------------------------------------------------------
# Whole-file report surface: kept/dropped tally, columns, TSV
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_mixed_records_end_to_end_matches_hand_derived_counts(tmp_path):
    path = sv.write_vcf(tmp_path, 'mixed.vcf', sv.MIXED_RECORDS, samples=sv.MIXED_SAMPLES, contigs=sv.MIXED_CONTIGS)
    df, report = vi.read_caller_vcf(path, 'testcaller', 'S1', sv.MIXED_CANONICAL, 'coh', 'hg38', 'illumina')
    assert report.records_read == len(sv.MIXED_RECORDS)
    assert report.kept == sv.MIXED_EXPECTED_KEPT
    assert report.dropped == sv.MIXED_EXPECTED_DROPPED
    assert df['vcf_id'].tolist() == sv.MIXED_KEPT_IDS
    assert list(df.columns) == vi.RAW_COLUMNS
    assert df['id'].is_unique


@pytest.mark.unit
def test_read_caller_vcf_output_has_exactly_raw_columns(tmp_path):
    df, _ = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ])
    assert list(df.columns) == vi.RAW_COLUMNS
    assert df['chrom_2'].isna().all()
    assert (df['cohort'] == 'mycohort').all()
    assert (df['sample'] == 'S1').all()
    assert (df['technology'] == 'illumina').all()
    assert (df['caller'] == 'testcaller').all()
    assert (df['reference'] == 'hg38').all()


@pytest.mark.unit
def test_read_caller_vcf_empty_result_still_has_raw_columns(tmp_path):
    df, report = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<INV>', 'info': 'SVTYPE=INV;END=200'},
    ])
    assert len(df) == 0
    assert list(df.columns) == vi.RAW_COLUMNS


@pytest.mark.unit
def test_read_caller_vcf_missing_file_sets_report_error(tmp_path):
    df, report = vi.read_caller_vcf(
        str(tmp_path / 'does_not_exist.vcf'), 'testcaller', 'S1', DEFAULT_CANONICAL, 'coh', 'hg38', 'illumina',
    )
    assert report.status == 'ERROR'
    assert report.error is not None
    assert len(df) == 0


@pytest.mark.unit
def test_reports_to_dataframe_has_stable_fixed_columns(tmp_path):
    _, report_a = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ], name='a.vcf')
    _, report_b = _read(tmp_path, [
        {'chrom': 'chr1', 'pos': 100, 'alt': '<INV>', 'info': 'SVTYPE=INV;END=200'},
    ], name='b.vcf')
    table = vi.reports_to_dataframe([report_a, report_b])
    assert list(table['status']) == ['OK', 'ZERO USABLE']
    for reason in vi.DROP_REASONS:
        assert f'dropped_{reason}' in table.columns
    for sv_type in vi.SUPPORTED_SV_TYPES:
        assert f'kept_{sv_type}' in table.columns


@pytest.mark.unit
def test_check_files_reads_every_file_in_order(tmp_path):
    path_a = sv.write_vcf(tmp_path, 'a.vcf', [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600'},
    ], samples=['S1'], contigs=DEFAULT_CONTIGS)
    path_b = sv.write_vcf(tmp_path, 'b.vcf', [
        {'chrom': 'chr1', 'pos': 100, 'id': 'B', 'alt': '<INS>', 'info': 'SVTYPE=INS;SVLEN=100'},
    ], samples=['S1'], contigs=DEFAULT_CONTIGS)
    reports = vi.check_files([('caller_a', path_a), ('caller_b', path_b)], DEFAULT_CANONICAL, 'S1')
    assert [r.caller for r in reports] == ['caller_a', 'caller_b']
    assert [r.usable for r in reports] == [1, 1]


# ---------------------------------------------------------------------------
# write_dq_tagged_vcf
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_write_dq_tagged_vcf_tags_records_by_internal_ordinal_id(tmp_path):
    path = sv.write_vcf(tmp_path, 'in.vcf', [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
         'samples': {'S1': '0/1'}},
        {'chrom': 'chr1', 'pos': 700, 'id': 'B', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=900',
         'samples': {'S1': '0/1'}},
    ], samples=['S1'], contigs=DEFAULT_CONTIGS)
    out_path = str(tmp_path / 'out.vcf')

    vi.write_dq_tagged_vcf(path, out_path, 'testcaller', {'testcaller:0': 0.87})

    out = pysam.VariantFile(out_path)
    records = list(out)
    assert records[0].info['DQ'] == '0.87'
    # No score entry for ordinal 1 -> -1, matching today's add_info_tag_to_vcf.
    assert records[1].info['DQ'] == '-1'


@pytest.mark.unit
def test_write_dq_tagged_vcf_handles_undeclared_contig_and_filter(tmp_path):
    # Written with no ##contig line and a FILTER value the header never
    # declares -- exactly the cnvnator2VCF / raw-caller dialect this needs
    # to survive.
    path = sv.write_vcf(tmp_path, 'in.vcf', [
        {'chrom': 'chrUndeclared', 'pos': 100, 'id': 'A', 'alt': '<DEL>',
         'info': 'SVTYPE=DEL;END=600', 'filter': 'lowQ', 'samples': {'S1': '0/1'}},
    ], samples=['S1'], contigs=[])
    out_path = str(tmp_path / 'out.vcf')

    assert vi.write_dq_tagged_vcf(path, out_path, 'testcaller', {'testcaller:0': 0.5}) == 1

    rec = next(iter(pysam.VariantFile(out_path)))
    assert rec.chrom == 'chrUndeclared'
    assert rec.info['DQ'] == '0.5'
    assert list(rec.filter.keys()) == ['lowQ']


@pytest.mark.unit
def test_write_dq_tagged_vcf_is_a_verbatim_copy_plus_dq(tmp_path):
    # The output must be the user's file plus one tag: htslib would rewrite
    # END for symbolic ALTs and lose GT phasing, so the tagger works on text.
    path = sv.write_vcf(tmp_path, 'in.vcf', [
        {'chrom': 'chr1', 'pos': 1, 'id': 'CNVnator_del_1', 'alt': '<DEL>',
         'info': 'END=10000;SVTYPE=DEL;SVLEN=-10000;IMPRECISE', 'samples': {'S1': '0|1'}},
        {'chrom': 'chr1', 'pos': 500, 'id': 'NOINFO', 'alt': '<DEL>', 'info': '.',
         'samples': {'S1': '0/1'}},
    ], samples=['S1'], contigs=DEFAULT_CONTIGS)
    out_path = str(tmp_path / 'out.vcf')

    vi.write_dq_tagged_vcf(path, out_path, 'c', {'c:0': 0.25})

    in_lines = open(path).read().splitlines()
    out_lines = open(out_path).read().splitlines()
    in_header = [l for l in in_lines if l.startswith('##')]
    out_header = [l for l in out_lines if l.startswith('##')]
    assert out_header == in_header + ['##INFO=<ID=DQ,Number=1,Type=String,Description="Dicast Quality Score">']
    in_records = [l for l in in_lines if not l.startswith('#')]
    out_records = [l for l in out_lines if not l.startswith('#')]
    assert out_records[0] == in_records[0].replace('SVLEN=-10000;IMPRECISE', 'SVLEN=-10000;IMPRECISE;DQ=0.25')
    assert out_records[1].split('\t')[7] == 'DQ=-1'
    assert out_records[1].split('\t')[9] == '0/1'


@pytest.mark.unit
def test_write_dq_tagged_vcf_does_not_duplicate_pre_existing_dq_field(tmp_path):
    path = sv.write_vcf(
        tmp_path, 'in.vcf',
        [{'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
          'samples': {'S1': '0/1'}}],
        samples=['S1'], contigs=DEFAULT_CONTIGS,
        extra_header=['##INFO=<ID=DQ,Number=1,Type=String,Description="pre-existing">'],
    )
    out_path = str(tmp_path / 'out.vcf')

    vi.write_dq_tagged_vcf(path, out_path, 'testcaller', {'testcaller:0': 0.1})

    header_text = str(pysam.VariantFile(out_path).header)
    assert header_text.count('##INFO=<ID=DQ,') == 1
    rec = next(iter(pysam.VariantFile(out_path)))
    assert rec.info['DQ'] == '0.1'


@pytest.mark.unit
def test_write_dq_tagged_vcf_preserves_genotypes(tmp_path):
    path = sv.write_vcf(tmp_path, 'in.vcf', [
        {'chrom': 'chr1', 'pos': 100, 'id': 'A', 'alt': '<DEL>', 'info': 'SVTYPE=DEL;END=600',
         'samples': {'S1': '0/1'}},
    ], samples=['S1'], contigs=DEFAULT_CONTIGS)
    out_path = str(tmp_path / 'out.vcf')

    vi.write_dq_tagged_vcf(path, out_path, 'testcaller', {})

    rec = next(iter(pysam.VariantFile(out_path)))
    assert tuple(rec.samples['S1']['GT']) == (0, 1)
