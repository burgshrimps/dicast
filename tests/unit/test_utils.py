"""Unit tests for :mod:`dicast.utils`.

These cover the pure helpers (``compute_overlap``, ``cigartuples_to_array``,
``pad_alignment_matrices``, ``mad``, ``replace_filename``, ``read_parameters``),
which account for the bulk of the module's executable lines that are testable
without a real BAM.

The BAM-/file-heavy functions ``compute_aln_matrix``, ``compute_cov_df`` and
``compute_rep_df`` are intentionally out of scope: the first needs a real
indexed BAM fetch, and the latter two call ``replace_filename`` with the wrong
arity (a single ``params`` dict instead of ``filename, sample, ref``), which
raises ``TypeError`` before doing anything testable.

NOTE on divergence from the older lucid/dicast dev line this was ported from:
this repo removed cohort VCF/CSV mode entirely (replaced by the ``multi``
subcommand's cross-sample rescue in ``dicast/multi.py``), so
``dicast/utils.py`` no longer defines ``sample_vcf_to_dataframe`` at all.
All tests for that function (and its supporting ``make_sample_vcf`` fixture)
are dropped.

``caller_vcf_to_dataframe`` -- DROPPED along with the function itself. VCF
reading now lives entirely in ``dicast/vcf_input.py::read_caller_vcf`` (see
``tests/unit/test_vcf_input.py``); ``dicast/prepare.py`` was the only caller
of ``caller_vcf_to_dataframe`` and now calls ``read_caller_vcf`` instead.

``dicast`` imports cleanly, so ``utils`` is imported directly.
"""
from __future__ import annotations

import json

import numpy as np
import pytest

from dicast import utils


# ---------------------------------------------------------------------------
# read_parameters(file)  ->  json.load of the file
#
# This function is still defined in dicast/utils.py, but a repo-wide grep
# turned up no callers left in dicast/ or dicast.py -- it looks like dead
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

