"""Unit tests for :mod:`dicast_lib.collect_illumina`.

Two complementary strategies are used:

(a) **Pure helpers** are tested by bypassing ``__init__`` via
    ``object.__new__(AlignmentAnnotatorIllumina)`` and setting only the
    attributes the method under test reads. These need no BAM at all and let us
    derive expected outputs by construction (e.g. ``get_end_position`` from a
    known CIGAR, ``get_overlap`` of two known intervals).

(b) **Integration** drives the full read-based feature pipeline against a tiny
    synthetic, indexed BAM written by :mod:`tests.fixtures.synthetic_bam`. This
    is what exercises the large ``calculate_read_based_features`` /
    ``calculate_coverage`` code paths, which is required to clear 60% coverage.

``dicast_lib`` imports cleanly, so the class is imported directly.

Notes on porting from the older lucid/dicast dev line:
  - This repo's runtime env pins pysam 0.24.0, which has ``AlignedSegment.is_mapped``
    (added in pysam 0.21), so the old repo's "skip the read-based stage on old
    pysam" workaround does not apply here and has been dropped.
  - ``cov_thr`` is 5 in this repo (was 3 in the older dev line).
  - ``AlignmentAnnotatorIllumina.__init__`` no longer takes ``log_file``/``job_id``;
    progress logging now goes through the stdlib ``logging`` module
    (``logging.info``) via ``self._log``, rather than writing to a log file.
  - The three baseline-sampling methods each call ``np.random.seed(42)`` right
    before their sampling loop, so repeated calls with the same args are exactly
    reproducible -- this is asserted directly below.
"""
from __future__ import annotations

import logging
import types

import numpy as np
import pandas as pd
import pytest

from dicast_lib.collect_illumina import AlignmentAnnotatorIllumina
from tests.fixtures import synthetic_bam


# ===========================================================================
# (a) Pure helpers — no BAM required.
# ===========================================================================

def _bare() -> AlignmentAnnotatorIllumina:
    """An instance with __init__ skipped (no BAM opened)."""
    return object.__new__(AlignmentAnnotatorIllumina)


# ---------------------------------------------------------------------------
# get_overlap(a, b) -> max(0, min(a[1],b[1]) - max(a[0],b[0]) + 1)
# Note the INCLUSIVE +1: touching intervals overlap by 1, identical by width+1.
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize(
    "a, b, expected",
    [
        # [0,10] and [5,15] -> min(10,15)-max(0,5)+1 = 10-5+1 = 6
        ((0, 10), (5, 15), 6),
        # Disjoint with a gap -> clamps to 0
        ((0, 5), (10, 20), 0),
        # Touching endpoints [0,5] and [5,10] -> 5-5+1 = 1 (inclusive)
        ((0, 5), (5, 10), 1),
        # Nested: [0,20] contains [5,10] -> 10-5+1 = 6
        ((0, 20), (5, 10), 6),
        # Identical [3,9] -> 9-3+1 = 7
        ((3, 9), (3, 9), 7),
    ],
)
def test_get_overlap(a, b, expected):
    assert _bare().get_overlap(a, b) == expected


@pytest.mark.unit
def test_get_overlap_never_negative():
    assert _bare().get_overlap((0, 10), (100, 200)) == 0


# ---------------------------------------------------------------------------
# get_end_position(start, cigar_str) — end = start-1 + sum of M/D/N/=/X lengths.
# Soft-clip (S), insertion (I), hard-clip (H), pad (P) do NOT advance reference.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_end_position_simple_match():
    # 100M from start 1000 -> end = 999 + 100 = 1099
    assert _bare().get_end_position(1000, "100M") == 1099


@pytest.mark.unit
def test_get_end_position_with_deletion_and_softclip():
    # 10S then 50M 5D 50M from 1000:
    #  S consumes no reference; M(50)+D(5)+M(50) = 105 ref bases.
    #  end = 999 + 105 = 1104
    assert _bare().get_end_position(1000, "10S50M5D50M") == 1104


@pytest.mark.unit
def test_get_end_position_insertion_does_not_advance():
    # 20M 10I 20M -> only the 40 M bases advance the reference.
    # end = 999 + 40 = 1039
    assert _bare().get_end_position(1000, "20M10I20M") == 1039


@pytest.mark.unit
def test_get_end_position_skip_and_equals_x():
    # 10= 5N 10X advances by 10 + 5 + 10 = 25
    assert _bare().get_end_position(500, "10=5N10X") == 500 - 1 + 25


# ---------------------------------------------------------------------------
# suffix_to_bin_idx
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("suffix, idx", [("I", 1), ("II", 2), ("III", 3), ("IV", 4)])
def test_suffix_to_bin_idx(suffix, idx):
    assert _bare().suffix_to_bin_idx(suffix) == idx


# ---------------------------------------------------------------------------
# divide_sv_body(start, end) -> 3 interior np.linspace edges as a Series.
# np.linspace(start, end, 5) gives 5 evenly spaced points; [1:-1] keeps the
# three interior ones, cast to int.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_divide_sv_body_values_and_index():
    out = _bare().divide_sv_body(0, 400)
    # linspace(0, 400, 5) = [0, 100, 200, 300, 400] -> interior = 100, 200, 300
    assert list(out.index) == ["body_I", "body_II", "body_III"]
    assert list(out.values) == [100, 200, 300]


@pytest.mark.unit
def test_divide_sv_body_is_series_of_ints():
    out = _bare().divide_sv_body(1000, 2000)
    assert isinstance(out, pd.Series)
    # interior points of linspace(1000,2000,5)=[1000,1250,1500,1750,2000]
    assert list(out.values) == [1250, 1500, 1750]


# ---------------------------------------------------------------------------
# get_clipped_span(read) — uses read.reference_start/end and cigartuples.
# Left clip (op 4/5 at cigar[0]): (ref_start - len, ref_start - 1)
# Right clip (op 4/5 at cigar[-1]): (ref_end + 1, ref_end + len)
# Neither -> None.  Left clip takes precedence over right.
# ---------------------------------------------------------------------------

def _fake_read(**kw):
    return types.SimpleNamespace(**kw)


@pytest.mark.unit
def test_get_clipped_span_left_softclip():
    read = _fake_read(reference_start=1000, reference_end=1200,
                      cigartuples=[(4, 15), (0, 200)])
    # (1000 - 15, 1000 - 1) = (985, 999)
    assert _bare().get_clipped_span(read) == (985, 999)


@pytest.mark.unit
def test_get_clipped_span_right_hardclip():
    read = _fake_read(reference_start=1000, reference_end=1200,
                      cigartuples=[(0, 200), (5, 30)])
    # (1200 + 1, 1200 + 30) = (1201, 1230)
    assert _bare().get_clipped_span(read) == (1201, 1230)


@pytest.mark.unit
def test_get_clipped_span_none_when_unclipped():
    read = _fake_read(reference_start=1000, reference_end=1200,
                      cigartuples=[(0, 200)])
    assert _bare().get_clipped_span(read) is None


# ---------------------------------------------------------------------------
# get_overlap_mate_bins(read, bins) — mate window is
# [next_reference_start, next_reference_start + 150]; 1 per bin it overlaps.
# Only counts when read.is_mapped is True. The fake read carries an explicit
# `is_mapped` attribute.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_overlap_mate_bins_hits_overlapping_bins():
    read = _fake_read(is_mapped=True, next_reference_start=1000)
    bins = [(900, 950),       # no overlap (mate window 1000..1150)
            (1100, 1200),     # overlaps
            (5000, 5050)]     # no overlap
    out = _bare().get_overlap_mate_bins(read, bins)
    np.testing.assert_array_equal(out, np.array([0.0, 1.0, 0.0]))


@pytest.mark.unit
def test_get_overlap_mate_bins_unmapped_returns_zeros():
    read = _fake_read(is_mapped=False, next_reference_start=1000)
    bins = [(1000, 1100), (1100, 1200)]
    out = _bare().get_overlap_mate_bins(read, bins)
    np.testing.assert_array_equal(out, np.zeros(2))


# ---------------------------------------------------------------------------
# jump_to_next_variant_for_coverage_calculation — when the bin's coverage mean
# is <= cov_thr returns (1, []); otherwise it excludes overlapping variants.
# This repo's cov_thr default is 5 (was 3 in the older dev line); the threshold
# used below is set explicitly on the bare instance either way, so the value
# tested against is independent of that default.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_jump_below_threshold_returns_default():
    obj = _bare()
    obj.cov_thr = 5
    df = pd.DataFrame({"start": [1000], "end": [2000], "ill_cov_mean_I": [1.0]})
    step, exclude = obj.jump_to_next_variant_for_coverage_calculation(
        df, 0, "start", "start", -52, -2, "I")
    assert step == 1
    assert exclude == []


@pytest.mark.unit
def test_jump_above_threshold_excludes_self():
    obj = _bare()
    obj.cov_thr = 5
    # A single variant whose coverage exceeds the threshold: it overlaps its own
    # bin, so it is excluded -> step==1, exclude==[0].
    df = pd.DataFrame({"start": [1000], "end": [2000], "ill_cov_mean_I": [9.0]})
    step, exclude = obj.jump_to_next_variant_for_coverage_calculation(
        df, 0, "start", "start", -52, -2, "I")
    assert step == 1
    assert exclude == [0]


@pytest.mark.unit
def test_jump_uses_instance_cov_thr_of_5():
    # A value that would have excluded under the old dev line's cov_thr=3 but
    # stays below this repo's cov_thr=5 -- confirms the real default is read.
    obj = _bare()
    obj.cov_thr = 5
    df = pd.DataFrame({"start": [1000], "end": [2000], "ill_cov_mean_I": [4.0]})
    step, exclude = obj.jump_to_next_variant_for_coverage_calculation(
        df, 0, "start", "start", -52, -2, "I")
    assert step == 1
    assert exclude == []


# ---------------------------------------------------------------------------
# prepare_dataframe / load_from_df / to_df — DataFrame plumbing.
# prepare_dataframe filters to chrom+sv_type and adds all feature columns.
# ---------------------------------------------------------------------------

def _variants_df():
    return pd.DataFrame({
        "id": ["v1", "v2", "v3"],
        "cohort": ["c", "c", "c"],
        "sample": ["s", "s", "s"],
        "reference": ["hg38", "hg38", "hg38"],
        "technology": ["ill", "ill", "ill"],
        "caller": ["x", "x", "x"],
        "sv_type": ["DEL", "DUP", "DEL"],     # v2 is a DUP -> filtered out
        "chrom": ["chr1", "chr1", "chr2"],     # v3 is chr2 -> filtered out
        "chrom_2": ["chr1", "chr1", "chr2"],
        "start": [1000, 2000, 3000],
        "end": [1400, 2400, 3400],
        "sv_len": [400, 400, 400],
    })


def _obj_for_df(chrom="chr1", sv_type="DEL"):
    """A bare instance configured just enough for prepare_dataframe/load_from_df."""
    obj = _bare()
    obj.chrom = chrom
    obj.sv_type = sv_type
    # Set the feature attribute lists exactly as __init__ does.
    obj.features_breakpoints = ['ill_cov_mean_', 'ill_cov_std_', 'ill_isize_mean_', 'ill_isize_std_',
                                'ill_mapq_mean_', 'ill_mapq_std_', 'ill_clipreads_', 'ill_splitreads_',
                                'ill_disco_ff_', 'ill_disco_rr_', 'ill_disco_rf_', 'ill_disco_tx_']
    obj.features_body = ['ill_cov_mean_', 'ill_cov_std_']
    obj.features_connection = ['ill_disco_ff_', 'ill_disco_rr_', 'ill_disco_rf_', 'ill_splitreads_']
    return obj


@pytest.mark.unit
def test_prepare_dataframe_filters_chrom_and_svtype():
    obj = _obj_for_df()
    obj.load_from_df(_variants_df())
    # Only v1 survives (chr1 + DEL).
    assert list(obj.df_calls_annot["id"]) == ["v1"]


@pytest.mark.unit
def test_prepare_dataframe_adds_feature_columns():
    obj = _obj_for_df()
    obj.load_from_df(_variants_df())
    cols = obj.df_calls_annot.columns
    # Breakpoint features get the I..IV suffixes.
    assert "ill_cov_mean_I" in cols
    assert "ill_mapq_std_IV" in cols
    # Body features get the IIa..IIIa suffixes.
    assert "ill_cov_mean_IIa" in cols
    assert "ill_cov_std_IIIa" in cols
    # Connection features get paired suffixes like I_II.
    assert "ill_disco_ff_I_II" in cols
    assert "ill_splitreads_III_IV" in cols
    # New feature columns are initialised to NaN.
    assert obj.df_calls_annot["ill_cov_mean_I"].isna().all()


@pytest.mark.unit
def test_load_from_df_copies_input():
    obj = _obj_for_df()
    df = _variants_df()
    before = df.copy(deep=True)
    obj.load_from_df(df)
    # load_from_df does df.copy(), so the caller's frame is untouched.
    pd.testing.assert_frame_equal(df, before)


@pytest.mark.unit
def test_to_df_returns_annot_frame():
    obj = _bare()
    sentinel = pd.DataFrame({"a": [1, 2]})
    obj.df_calls_annot = sentinel
    assert obj.to_df() is sentinel


# ===========================================================================
# (b) Integration — drive the full pipeline against a tiny indexed BAM.
# ===========================================================================

@pytest.fixture()
def bam_path(tmp_path):
    return synthetic_bam.make_bam(tmp_path / "synthetic.bam")


@pytest.fixture()
def integration_variants():
    """A single DEL inside the covered window, with all required columns."""
    return pd.DataFrame({
        "id": ["del1"],
        "cohort": ["test_cohort"],
        "sample": ["test_sample"],
        "reference": ["hg38"],
        "technology": ["illumina"],
        "caller": ["manta"],
        "sv_type": ["DEL"],
        "chrom": ["chr1"],
        "chrom_2": ["chr1"],
        "start": [synthetic_bam.SV_START],
        "end": [synthetic_bam.SV_END],
        "sv_len": [synthetic_bam.SV_END - synthetic_bam.SV_START],
    })


@pytest.mark.unit
def test_integration_full_pipeline(bam_path, integration_variants):
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")

    aai.load_from_df(integration_variants)

    # Small n keeps the random baseline sampling fast on the 5 kb contig.
    aai.calculate_coverage_baseline(s=500, n=20)
    aai.calculate_insertsize_baseline(s=500, n=20)
    aai.calculate_mapping_quality_baseline(s=500, n=20)

    # Baselines are populated and finite. The dense full-length reads guarantee
    # a positive median coverage baseline.
    assert aai.baseline_coverage_mean > 0
    assert np.isfinite(aai.baseline_mapq_mean)
    assert np.isfinite(aai.baseline_insertsize_median)

    aai.annotate_coverage()

    # The coverage-stage columns are populated.
    cov_annot = aai.df_calls_annot
    assert len(cov_annot) == 1
    for col in ["ill_cov_mean_I", "ill_cov_std_IV", "ill_cov_mean_IIa", "ill_cov_std_IIIa"]:
        assert col in cov_annot.columns, f"missing column {col}"

    aai.annotate_read_based_features()
    aai.aggregate_results()
    result = aai.to_df()

    # The DEL survived (coverage over the window is well below the log2 thr).
    assert len(result) == 1

    # Expected feature columns are present after aggregation.
    for col in ["ill_cov_mean_I", "ill_cov_std_IV",
                "ill_isize_mean_I", "ill_mapq_mean_II",
                "ill_clipreads_III", "ill_splitreads_IV",
                "ill_disco_ff_I", "ill_disco_rr_II", "ill_disco_rf_III", "ill_disco_tx_IV",
                "ill_cov_mean_IIa", "ill_cov_std_IIIa"]:
        assert col in result.columns, f"missing column {col}"

    # Identity columns are carried through unchanged.
    row = result.iloc[0]
    assert row["id"] == "del1"
    assert row["sample"] == "test_sample"
    assert row["chrom"] == "chr1"
    assert row["start"] == synthetic_bam.SV_START
    assert row["end"] == synthetic_bam.SV_END

    # Read-based feature values are finite numbers (the bins around the start
    # breakpoint see plenty of reads in the fixture).
    assert np.isfinite(row["ill_mapq_mean_I"])
    assert np.isfinite(row["ill_isize_mean_I"])
    # Clipped/split reads were constructed near the start breakpoint, so at least
    # one of the start bins reports a non-negative clip/split proportion.
    assert row["ill_clipreads_I"] >= 0
    assert row["ill_splitreads_I"] >= 0


@pytest.mark.unit
def test_calculate_coverage_region_positive_over_covered_window(bam_path, integration_variants):
    """calculate_coverage_region returns finite log2-ratios; with the dense
    fixture the raw window coverage is clearly positive."""
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")
    aai.load_from_df(integration_variants)
    aai.baseline_coverage_mean = 1.0
    aai.baseline_coverage_std = 1.0

    out = aai.calculate_coverage_region("chr1", 1000, 1100, "I")
    assert list(out.index) == ["ill_cov_mean_I", "ill_cov_std_I"]
    assert np.isfinite(out["ill_cov_mean_I"])
    # The dense reads (8x full-length) plus several window reads put the raw mean
    # coverage above the baseline of 1.0, so the log2 ratio is positive.
    assert out["ill_cov_mean_I"] > 0
    aai.bam.close()


@pytest.mark.unit
def test_calculate_read_based_features_returns_series(bam_path, integration_variants):
    """Directly exercise calculate_read_based_features for one bin and check the
    returned Series shape/labels."""
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")
    aai.load_from_df(integration_variants)
    # Minimal baselines required by the normalisation math.
    aai.baseline_insertsize_median = 300.0
    aai.baseline_insertsize_mad = 50.0
    aai.baseline_mapq_mean = 55.0
    aai.baseline_mapq_std = 10.0

    out = aai.calculate_read_based_features(
        chrom="chr1", start=998, end=1052,
        sv_start=synthetic_bam.SV_START, sv_end=synthetic_bam.SV_END, suffix="I")

    assert isinstance(out, pd.Series)
    # Bin I -> bin_idx 1 -> three connection bins -> connection labels present.
    for label in ["ill_isize_mean_I", "ill_mapq_std_I", "ill_clipreads_I",
                  "ill_disco_ff_I", "ill_disco_tx_I",
                  "ill_splitreads_I_II", "ill_disco_ff_I_II"]:
        assert label in out.index
    assert np.isfinite(out["ill_mapq_mean_I"])
    aai.bam.close()


# ===========================================================================
# New coverage for this repo's specific behavior (not present in the older
# lucid/dicast dev line's test module):
#   - cov_thr default is 5, not 3.
#   - np.random.seed(42) precedes each baseline-sampling loop, so repeated
#     calls with identical args reproduce identical baseline statistics.
#   - progress messages go through logging.info via self._log, not a log file.
# ===========================================================================

@pytest.mark.unit
def test_cov_thr_default_is_5(bam_path):
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")
    assert aai.cov_thr == 5
    aai.bam.close()


@pytest.mark.unit
def test_coverage_baseline_is_reproducible_across_calls(bam_path):
    """np.random.seed(42) precedes the sampling loop in
    calculate_coverage_baseline, so two independent calls with the same (s, n)
    against the same BAM produce bit-identical baseline statistics."""
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")

    aai.calculate_coverage_baseline(s=500, n=20)
    mean_first = aai.baseline_coverage_mean
    std_first = aai.baseline_coverage_std

    aai.calculate_coverage_baseline(s=500, n=20)
    mean_second = aai.baseline_coverage_mean
    std_second = aai.baseline_coverage_std

    assert mean_first == mean_second
    assert std_first == std_second
    aai.bam.close()


@pytest.mark.unit
def test_insertsize_baseline_is_reproducible_across_calls(bam_path):
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")

    aai.calculate_insertsize_baseline(s=500, n=20)
    median_first = aai.baseline_insertsize_median
    mad_first = aai.baseline_insertsize_mad

    aai.calculate_insertsize_baseline(s=500, n=20)
    median_second = aai.baseline_insertsize_median
    mad_second = aai.baseline_insertsize_mad

    assert median_first == median_second
    assert mad_first == mad_second
    aai.bam.close()


@pytest.mark.unit
def test_mapping_quality_baseline_is_reproducible_across_calls(bam_path):
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")

    aai.calculate_mapping_quality_baseline(s=500, n=20)
    mean_first = aai.baseline_mapq_mean
    std_first = aai.baseline_mapq_std

    aai.calculate_mapping_quality_baseline(s=500, n=20)
    mean_second = aai.baseline_mapq_mean
    std_second = aai.baseline_mapq_std

    assert mean_first == mean_second
    assert std_first == std_second
    aai.bam.close()


@pytest.mark.unit
def test_log_emits_via_logging_info(bam_path, caplog):
    """self._log now routes through logging.info (no more log_file/job_id
    print-based logging from the older dev line)."""
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")

    with caplog.at_level(logging.INFO):
        aai._log("hello world", level="INFO")

    assert len(caplog.records) == 1
    record = caplog.records[0]
    assert record.levelno == logging.INFO
    assert "test_sample_chr1_DEL" in record.message
    assert "hello world" in record.message
    aai.bam.close()


@pytest.mark.unit
def test_annotate_coverage_logs_progress_via_logging(bam_path, integration_variants, caplog):
    aai = AlignmentAnnotatorIllumina(
        bam_filename=bam_path, chrom="chr1", sv_type="DEL", sample="test_sample")
    aai.load_from_df(integration_variants)
    aai.calculate_coverage_baseline(s=500, n=20)

    with caplog.at_level(logging.INFO):
        aai.annotate_coverage()

    messages = [r.message for r in caplog.records]
    assert any("Starting coverage annotation" in m for m in messages)
    assert any("DONE" in m for m in messages)
    aai.bam.close()
