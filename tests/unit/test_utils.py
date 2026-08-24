"""Unit tests for :mod:`dicast_lib.utils`.

These cover the pure helpers (``compute_overlap``, ``cigartuples_to_array``,
``pad_alignment_matrices``, ``mad``, ``replace_filename``, ``read_parameters``)
and the VCF-parser function (``caller_vcf_to_dataframe``), which together
account for the bulk of the module's executable lines.

The BAM-/file-heavy functions ``compute_aln_matrix``, ``compute_cov_df`` and
``compute_rep_df`` are intentionally out of scope: the first needs a real
indexed BAM fetch, and the latter two call ``replace_filename`` with the wrong
arity (a single ``params`` dict instead of ``filename, sample, ref``), which
raises ``TypeError`` before doing anything testable.

NOTE on divergence from the older lucid/dicast dev line this was ported from:
this repo removed cohort VCF/CSV mode entirely (replaced by the ``multi``
subcommand's cross-sample rescue in ``dicast_lib/multi.py``), so
``dicast_lib/utils.py`` no longer defines ``sample_vcf_to_dataframe`` at all.
All tests for that function (and its supporting ``make_sample_vcf`` fixture)
are dropped -- see the module docstring in ``tests/fixtures/synthetic_vcf.py``
for the matching fixture-side removal.

``dicast_lib`` imports cleanly, so ``utils`` is imported directly. EXPECTED
values for the VCF tests are derived from how the fixture is constructed
(``tests/fixtures/synthetic_vcf.py``), never from running the parser first.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from dicast_lib import utils
from tests.fixtures import synthetic_vcf as sv


# ---------------------------------------------------------------------------
# read_parameters(file)  ->  json.load of the file
#
# This function is still defined in dicast_lib/utils.py, but a repo-wide grep
# turned up no callers left in dicast_lib/ or dicast.py -- it looks like dead
# code left over from a removed parameter-file workflow. It is still simple,
# pure, and part of the module's public surface, so it is covered here as new
# coverage rather than dropped.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_read_parameters_roundtrips_json(tmp_path):
    params = {
        "sample": "S1",
        "reference": "hg38",
        "nested": {"a": 1, "b": [1, 2, 3]},
    }
    path = tmp_path / "params.json"
    path.write_text(json.dumps(params))
    out = utils.read_parameters(str(path))
    assert out == params


@pytest.mark.unit
def test_read_parameters_returns_dict(tmp_path):
    path = tmp_path / "p.json"
    path.write_text('{"k": "v"}')
    out = utils.read_parameters(str(path))
    assert isinstance(out, dict)
    assert out["k"] == "v"


# ---------------------------------------------------------------------------
# replace_filename(filename, sample, ref)
# Replaces the literal token 'SAMPLE' with `sample` and 'REF' with `ref`.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_replace_filename_substitutes_both_tokens():
    out = utils.replace_filename("SAMPLE.REF.bam", "NA12878", "hg38")
    assert out == "NA12878.hg38.bam"


@pytest.mark.unit
def test_replace_filename_replaces_all_occurrences():
    out = utils.replace_filename("dir/SAMPLE/SAMPLE_REF.vcf", "s1", "grch38")
    assert out == "dir/s1/s1_grch38.vcf"


@pytest.mark.unit
def test_replace_filename_no_tokens_is_unchanged():
    assert utils.replace_filename("plain_name.txt", "s1", "hg38") == "plain_name.txt"


# ---------------------------------------------------------------------------
# compute_overlap(s1, s2, e1, e2)  ->  max(0, min(e1, e2) - max(s1, s2))
# Note the argument order: both starts first, then both ends.
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize(
    "s1, s2, e1, e2, expected",
    [
        # Overlapping: [0,10) and [5,15) -> overlap of 5
        (0, 5, 10, 15, 5),
        # Disjoint (seg1 entirely left of seg2): [0,5) and [10,20) -> 0
        (0, 10, 5, 20, 0),
        # Touching (end of seg1 == start of seg2): [0,5) and [5,10) -> 0
        (0, 5, 5, 10, 0),
        # Nested: [0,20) fully contains [5,10) -> overlap is inner width 5
        (0, 5, 20, 10, 5),
        # Identical segments: [3,9) and [3,9) -> full width 6
        (3, 3, 9, 9, 6),
    ],
)
def test_compute_overlap(s1, s2, e1, e2, expected):
    assert utils.compute_overlap(s1, s2, e1, e2) == expected


@pytest.mark.unit
def test_compute_overlap_never_negative():
    # Fully disjoint with a gap must clamp to 0, not a negative number.
    assert utils.compute_overlap(0, 100, 10, 110) == 0


# ---------------------------------------------------------------------------
# cigartuples_to_array(cigartuples)
# Each (op, length) tuple expands to `op` repeated `length` times, flattened.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_cigartuples_to_array_typical():
    # 3M = op0 x3, 2I = op1 x2, 1D = op2 x1
    out = utils.cigartuples_to_array([(0, 3), (1, 2), (2, 1)])
    assert isinstance(out, np.ndarray)
    np.testing.assert_array_equal(out, np.array([0, 0, 0, 1, 1, 2]))


@pytest.mark.unit
def test_cigartuples_to_array_single_op():
    out = utils.cigartuples_to_array([(0, 4)])
    np.testing.assert_array_equal(out, np.array([0, 0, 0, 0]))


@pytest.mark.unit
def test_cigartuples_to_array_empty():
    out = utils.cigartuples_to_array([])
    assert isinstance(out, np.ndarray)
    assert out.shape == (0,)


@pytest.mark.unit
def test_cigartuples_to_array_length_is_sum_of_lengths():
    tuples = [(0, 10), (4, 5), (0, 3)]
    out = utils.cigartuples_to_array(tuples)
    assert len(out) == sum(length for _, length in tuples)


# ---------------------------------------------------------------------------
# pad_alignment_matrices(left, right)
# Pads the shorter matrix (by row count) with the -1 sentinel so both share the
# same number of rows. Column count is left unchanged.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_pad_alignment_matrices_right_shorter():
    left = np.zeros((5, 4))
    right = np.zeros((2, 4))
    out_left, out_right = utils.pad_alignment_matrices(left, right)
    assert out_left.shape == (5, 4)
    assert out_right.shape == (5, 4)
    # Original rows preserved, padded rows are the -1 sentinel.
    np.testing.assert_array_equal(out_right[:2], np.zeros((2, 4)))
    np.testing.assert_array_equal(out_right[2:], -1 * np.ones((3, 4)))


@pytest.mark.unit
def test_pad_alignment_matrices_left_shorter():
    left = np.zeros((1, 3))
    right = np.zeros((4, 3))
    out_left, out_right = utils.pad_alignment_matrices(left, right)
    assert out_left.shape == (4, 3)
    assert out_right.shape == (4, 3)
    np.testing.assert_array_equal(out_left[1:], -1 * np.ones((3, 3)))


@pytest.mark.unit
def test_pad_alignment_matrices_equal_heights_unchanged():
    left = np.ones((3, 2))
    right = np.zeros((3, 2))
    out_left, out_right = utils.pad_alignment_matrices(left, right)
    assert out_left.shape == (3, 2)
    assert out_right.shape == (3, 2)
    np.testing.assert_array_equal(out_left, np.ones((3, 2)))
    np.testing.assert_array_equal(out_right, np.zeros((3, 2)))


# ---------------------------------------------------------------------------
# mad(arr)  ->  median(|arr - median(arr)|)
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_mad_simple_array():
    # arr = [1, 2, 3, 4, 5]; median = 3; |dev| = [2, 1, 0, 1, 2]; median = 1.
    arr = np.array([1, 2, 3, 4, 5])
    assert utils.mad(arr) == pytest.approx(1.0)


@pytest.mark.unit
def test_mad_with_outlier():
    # arr = [1, 1, 2, 2, 100]; median = 2; |dev| = [1, 1, 0, 0, 98];
    # sorted |dev| = [0, 0, 1, 1, 98]; median = 1. The outlier barely moves it.
    arr = np.array([1, 1, 2, 2, 100])
    assert utils.mad(arr) == pytest.approx(1.0)


@pytest.mark.unit
def test_mad_constant_array_is_zero():
    # No deviation from the median anywhere -> MAD is 0.
    assert utils.mad(np.array([7, 7, 7, 7])) == pytest.approx(0.0)


@pytest.mark.unit
def test_mad_even_length_uses_interpolated_median():
    # arr = [1, 2, 3, 4]; median = 2.5; |dev| = [1.5, 0.5, 0.5, 1.5];
    # sorted = [0.5, 0.5, 1.5, 1.5]; median = 1.0.
    arr = np.array([1, 2, 3, 4])
    assert utils.mad(arr) == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# caller_vcf_to_dataframe(vcf, cohort, sample, reference, technology, caller,
#                         canonical_chroms, check_genotype=False)
#
# Builds a DataFrame from a pysam.VariantFile. See synthetic_vcf.py for the
# record set: 5 records survive both the canonical-chrom and supported-svtype
# filters (DEL1, INS1, INV1, DUP1, BND1); a non-canonical DEL and an
# unsupported CNV are dropped.
#
# Unchanged from the lucid dev line: caller_vcf_to_dataframe's signature and
# per-svtype branches (DEL/INS/INV/DUP/BND) are identical in this repo, so
# these tests port over test-by-test with no behavioral adaptation needed.
# The BND branch is still exercised here at the utils-module unit level even
# though the CLI's default sv_types list no longer reaches BND end-to-end.
# ---------------------------------------------------------------------------

def _caller_df(tmp_path, *, check_genotype=False, name="caller.vcf"):
    vcf = sv.make_caller_vcf(tmp_path, name=name)
    return utils.caller_vcf_to_dataframe(
        vcf,
        cohort="mycohort",
        sample=sv.SAMPLE,
        reference="hg38",
        technology="illumina",
        caller="manta",
        canonical_chroms=sv.CANONICAL_CHROMS,
        check_genotype=check_genotype,
    )


@pytest.mark.unit
def test_caller_vcf_keeps_only_canonical_supported_records(tmp_path):
    df = _caller_df(tmp_path)
    # Non-canonical chrUn DEL and unsupported CNV are filtered out.
    assert len(df) == len(sv.KEPT_IDS)
    assert set(df["id"]) == set(sv.KEPT_IDS)


@pytest.mark.unit
def test_caller_vcf_has_expected_columns(tmp_path):
    df = _caller_df(tmp_path)
    expected_cols = {
        "id", "sv_type", "chrom", "start", "chrom_2", "end", "sv_len",
        "filter", "qual", "genotype", "cohort", "sample", "technology",
        "caller", "reference",
    }
    assert expected_cols.issubset(set(df.columns))


@pytest.mark.unit
def test_caller_vcf_constant_columns_filled_from_args(tmp_path):
    df = _caller_df(tmp_path)
    assert (df["cohort"] == "mycohort").all()
    assert (df["sample"] == sv.SAMPLE).all()
    assert (df["technology"] == "illumina").all()
    assert (df["caller"] == "manta").all()
    assert (df["reference"] == "hg38").all()


@pytest.mark.unit
def test_caller_vcf_deletion_coordinates(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    src = sv.record_by_id("DEL1")
    row = df.loc["DEL1"]
    assert row["sv_type"] == "DEL"
    assert row["chrom"] == src["chrom"]
    # start is the 1-based VCF POS (rec.start + 1).
    assert row["start"] == src["pos"]
    # DEL end is the INFO END (rec.stop); caller sv_len = end - start.
    assert row["end"] == src["end"]
    assert row["sv_len"] == src["end"] - src["pos"]
    assert pd.isna(row["chrom_2"])


@pytest.mark.unit
def test_caller_vcf_insertion_uses_svlen_and_start_plus_two(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    src = sv.record_by_id("INS1")
    row = df.loc["INS1"]
    assert row["sv_type"] == "INS"
    assert row["start"] == src["pos"]
    # Insertions: end = rec.start + 2 = (pos - 1) + 2 = pos + 1.
    assert row["end"] == src["pos"] + 1
    # sv_len comes straight from the INFO SVLEN value.
    assert row["sv_len"] == src["svlen"]
    assert pd.isna(row["chrom_2"])


@pytest.mark.unit
def test_caller_vcf_inversion_coordinates(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    src = sv.record_by_id("INV1")
    row = df.loc["INV1"]
    assert row["sv_type"] == "INV"
    assert row["start"] == src["pos"]
    assert row["end"] == src["end"]
    assert row["sv_len"] == src["end"] - src["pos"]


@pytest.mark.unit
def test_caller_vcf_duplication_coordinates(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    src = sv.record_by_id("DUP1")
    row = df.loc["DUP1"]
    assert row["sv_type"] == "DUP"
    assert row["start"] == src["pos"]
    assert row["end"] == src["end"]
    assert row["sv_len"] == src["end"] - src["pos"]


@pytest.mark.unit
def test_caller_vcf_breakend_parses_mate_from_alt(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    row = df.loc["BND1"]
    assert row["sv_type"] == "BND"
    # ALT was 'N[chr5:9000[' -> chrom_2 'chr5', end '9000' (kept as a string).
    assert row["chrom_2"] == "chr5"
    assert str(row["end"]) == "9000"
    assert pd.isna(row["sv_len"])


@pytest.mark.unit
def test_caller_vcf_genotype_is_record_gt(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    # INS1 was written 0/0 and DEL1 0/1; pysam returns GT as an allele tuple.
    assert tuple(df.loc["DEL1"]["genotype"]) == (0, 1)
    assert tuple(df.loc["INS1"]["genotype"]) == (0, 0)
    assert tuple(df.loc["INV1"]["genotype"]) == (1, 1)


@pytest.mark.unit
def test_caller_vcf_check_genotype_drops_absent_samples(tmp_path):
    df = _caller_df(tmp_path, check_genotype=True)
    # Only records whose GT contains allele 1 survive; INS1 (0/0) drops out.
    assert set(df["id"]) == set(sv.KEPT_PRESENT_IDS)
    assert "INS1" not in set(df["id"])


@pytest.mark.unit
def test_caller_vcf_filter_and_qual(tmp_path):
    df = _caller_df(tmp_path).set_index("id")
    # All records were written with FILTER=PASS and QUAL=50.
    assert (df["filter"] == "PASS").all()
    assert df.loc["DEL1"]["qual"] == pytest.approx(50.0)


# ---------------------------------------------------------------------------
# sample_vcf_to_dataframe -- DROPPED.
#
# This repo removed cohort VCF/CSV mode entirely (replaced by the `multi`
# subcommand's cross-sample rescue in dicast_lib/multi.py); dicast_lib/utils.py
# no longer defines sample_vcf_to_dataframe at all (confirmed by reading the
# module and by a repo-wide grep turning up zero definitions). All tests that
# exercised it in the lucid dev line -- keeps-only-canonical-supported-records,
# has-sample-specific-columns, deletion/insertion/breakend coordinate tests,
# caller lowercased/comma-joined, cohort_ac int coercion, single-supp-sample
# roundtrip, check_genotype drop, and constant-columns-filled -- are dropped
# here rather than adapted, since there is no function left to call.
# ---------------------------------------------------------------------------
