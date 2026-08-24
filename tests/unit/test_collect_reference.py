"""Unit tests for :mod:`dicast.collect_reference`.

Two complementary strategies are used:

* **Pure helpers** are driven on a bare ``object.__new__(ReferenceAnnotator)``
  instance so the file-reading constructor is skipped. The relevant attributes
  (``df_calls_annot``, ``bw_gc``, the BND split flag, ...) are set by hand and the
  method's output is asserted against values derived purely from construction.
* **Integration** builds a *real* annotator from tiny synthetic reference tracks
  (see :mod:`tests.fixtures.synthetic_reference`) and drives the full
  ``load_from_df -> split_bnd -> annotate_* -> annotate_gc_content ->
  aggregate_results -> to_df`` chain over a single call whose coordinates are
  laid out to overlap several tracks, so the emitted distances are known.

This repo's ``ReferenceAnnotator`` has no ``annotate_genes`` / ``annotate_orphanet``
methods at all (removed along with the gene/orphanet reference-file inputs), so
there is nothing to test or skip for them here -- they simply do not exist.

``collect_reference`` imports ``pyBigWig`` and ``bioframe`` at module top, so it
is only importable in the CI venv that provides them.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dicast.collect_reference import ReferenceAnnotator
from tests.fixtures import synthetic_reference as sref


# ---------------------------------------------------------------------------
# Bare-instance helper: skip the file-reading __init__.
# ---------------------------------------------------------------------------

def _blank_annotator() -> ReferenceAnnotator:
    """A ``ReferenceAnnotator`` with ``__init__`` skipped (no files read)."""
    return object.__new__(ReferenceAnnotator)


# ---------------------------------------------------------------------------
# load_from_df: copies the 11 named columns into df_calls / df_calls_annot.
# ---------------------------------------------------------------------------

def _full_calls_df():
    return pd.DataFrame(
        {
            "id": ["a", "b"],
            "cohort": ["c", "c"],
            "sample": ["s", "s"],
            "reference": ["GRCh38", "GRCh38"],
            "technology": ["ONT", "ONT"],
            "caller": ["sniffles", "sniffles"],
            "sv_type": ["DEL", "INS"],
            "chrom": ["chr1", "chr2"],
            "chrom_2": ["chr1", "chr2"],
            "start": [100, 200],
            "end": [150, 250],
            # an extra column that must NOT be carried into df_calls_annot
            "extra": ["x", "y"],
        }
    )


@pytest.mark.unit
def test_load_from_df_selects_expected_columns():
    obj = _blank_annotator()
    df = _full_calls_df()
    obj.load_from_df(df)

    expected_cols = [
        "id", "cohort", "sample", "reference", "technology", "caller",
        "sv_type", "chrom", "chrom_2", "start", "end",
    ]
    assert list(obj.df_calls_annot.columns) == expected_cols
    # The 'extra' column is dropped from the annotation frame but kept on df_calls.
    assert "extra" not in obj.df_calls_annot.columns
    assert "extra" in obj.df_calls.columns


@pytest.mark.unit
def test_load_from_df_copies_do_not_alias_input():
    obj = _blank_annotator()
    df = _full_calls_df()
    obj.load_from_df(df)
    # Mutating the source frame must not bleed into the stored copy.
    df.loc[0, "start"] = 999_999
    assert obj.df_calls.loc[0, "start"] == 100


# ---------------------------------------------------------------------------
# split_bnd: separates BND rows into bnd1 / bnd2 frames and sets the flag.
# NOTE: the method overwrites its own name with a bool, so each instance can
# only be split once.
# ---------------------------------------------------------------------------

def _annot_frame_with_bnd():
    return pd.DataFrame(
        {
            "id": ["d0", "b0"],
            "cohort": ["c", "c"],
            "sample": ["s", "s"],
            "reference": ["GRCh38", "GRCh38"],
            "technology": ["ONT", "ONT"],
            "caller": ["sniffles", "sniffles"],
            "sv_type": ["DEL", "BND"],
            "chrom": ["chr1", "chr1"],
            "chrom_2": ["chr1", "chr2"],
            "start": [1000, 2000],
            "end": [1100, 3000],
        }
    )


@pytest.mark.unit
def test_split_bnd_with_bnd_present():
    obj = _blank_annotator()
    obj.df_calls_annot = _annot_frame_with_bnd()
    obj.split_bnd()

    # The flag is now True (the method overwrote its own attribute).
    assert obj.split_bnd is True
    # Non-BND rows remain in the main frame.
    assert list(obj.df_calls_annot["sv_type"].unique()) == ["DEL"]

    # bnd1 is a +/-50 window around the BND start (2000) and drops chrom_2.
    bnd1 = obj.df_calls_annot_bnd1
    assert "chrom_2" not in bnd1.columns
    assert bnd1.loc[0, "start"] == 2000 - 50
    assert bnd1.loc[0, "end"] == 2000 + 50
    assert bnd1.loc[0, "chrom"] == "chr1"

    # bnd2 moves to chrom_2 and is a +/-50 window around the BND end (3000).
    bnd2 = obj.df_calls_annot_bnd2
    assert "chrom_2" not in bnd2.columns
    assert bnd2.loc[0, "chrom"] == "chr2"
    assert bnd2.loc[0, "start"] == 3000 - 50
    assert bnd2.loc[0, "end"] == 3000 + 50


@pytest.mark.unit
def test_split_bnd_without_bnd_sets_flag_false():
    obj = _blank_annotator()
    obj.df_calls_annot = pd.DataFrame(
        {
            "sv_type": ["DEL", "INS"],
            "chrom": ["chr1", "chr2"],
            "chrom_2": ["chr1", "chr2"],
            "start": [10, 20],
            "end": [30, 40],
        }
    )
    obj.split_bnd()
    assert obj.split_bnd is False
    # No BND-specific frames were created.
    assert not hasattr(obj, "df_calls_annot_bnd1")


# ---------------------------------------------------------------------------
# calculate_gc_content: reads self.bw_gc.stats(...); except Exception -> NaN.
# ---------------------------------------------------------------------------

class _FakeBigWig:
    """Minimal stand-in: ``stats`` returns a preset list, or raises."""

    def __init__(self, value, raises=False):
        self._value = value
        self._raises = raises
        self.calls = []

    def stats(self, chrom, start, end):
        self.calls.append((chrom, start, end))
        if self._raises:
            raise RuntimeError("boom")
        return [self._value]


@pytest.mark.unit
def test_calculate_gc_content_returns_stat_value():
    obj = _blank_annotator()
    obj.bw_gc = _FakeBigWig(0.55)
    row = pd.Series({"chrom": "chr1", "start": 1000, "end": 2000})
    out = obj.calculate_gc_content(row, "start")
    assert out == 0.55
    # The query window is row[region] +/- 50.
    assert obj.bw_gc.calls == [("chr1", 950, 1050)]


@pytest.mark.unit
def test_calculate_gc_content_uses_requested_region():
    obj = _blank_annotator()
    obj.bw_gc = _FakeBigWig(0.10)
    row = pd.Series({"chrom": "chrX", "start": 100, "end": 5000})
    obj.calculate_gc_content(row, "end")
    assert obj.bw_gc.calls == [("chrX", 4950, 5050)]


@pytest.mark.unit
def test_calculate_gc_content_exception_yields_nan():
    # This repo's calculate_gc_content narrowed the bare `except:` to
    # `except Exception:` -- a RuntimeError from stats() is still caught.
    obj = _blank_annotator()
    obj.bw_gc = _FakeBigWig(None, raises=True)
    row = pd.Series({"chrom": "chr1", "start": 1000, "end": 2000})
    out = obj.calculate_gc_content(row, "start")
    assert np.isnan(out)


# ---------------------------------------------------------------------------
# aggregate_results: no-op when split_bnd is falsy; merges BND frames otherwise.
# ---------------------------------------------------------------------------

# The exact column set aggregate_results' groupby aggregation references.
_REP_COLS = [
    "rep_LINE", "rep_SINE", "rep_LTR", "rep_DNA", "rep_Simple_repeat",
    "rep_Satellite", "rep_Low_complexity", "rep_Retroposon", "rep_snRNA",
    "rep_tRNA", "rep_srpRNA", "rep_rRNA", "rep_RC", "rep_scRNA", "rep_RNA",
    "rep_VNTR", "rep_STR",
]
_DIST_COLS = ["cpg_islands", "centromeres", "asmb_gaps", "alt_haps"]
_META_COLS = [
    "id", "cohort", "sample", "reference", "technology", "caller",
    "sv_type", "chrom", "chrom_2", "start", "end",
]


def _annotated_row(**overrides):
    """A fully-annotated single-row dict with every column aggregate_results reads."""
    row = {
        "id": "x", "cohort": "c", "sample": "s", "reference": "GRCh38",
        "technology": "ONT", "caller": "sniffles", "sv_type": "BND",
        "chrom": "chr1", "chrom_2": "chr2", "start": 1000, "end": 2000,
        "GC_content_left": 0.4, "GC_content_right": 0.6,
    }
    for col in _REP_COLS + _DIST_COLS:
        row[col] = 0
    row.update(overrides)
    return row


@pytest.mark.unit
def test_aggregate_results_noop_when_not_split():
    obj = _blank_annotator()
    obj.split_bnd = False
    before = pd.DataFrame([_annotated_row()])
    obj.df_calls_annot = before.copy()
    obj.aggregate_results()
    # Frame is untouched in the no-op branch.
    pd.testing.assert_frame_equal(obj.df_calls_annot, before)


@pytest.mark.unit
def test_aggregate_results_merges_bnd_frames():
    obj = _blank_annotator()
    obj.split_bnd = True

    # Main (non-BND) frame: one DEL with the same column layout as the BND output.
    main_cols = _META_COLS + _REP_COLS + _DIST_COLS + [
        "GC_content_left", "GC_content_right",
    ]
    main_row = _annotated_row(id="del0", sv_type="DEL")
    obj.df_calls_annot = pd.DataFrame([main_row])[main_cols]

    # Two halves of a single BND (shared id) -> aggregated back into one row.
    # split_bnd drops 'chrom_2' from both BND frames, and aggregate_results
    # recreates it by renaming bnd2's 'chrom' -> 'chrom_2'; so the BND frames
    # must NOT carry a 'chrom_2' column here (otherwise the rename collides).
    bnd1_row = _annotated_row(
        id="bnd0", start=1950, GC_content_left=0.2, GC_content_right=0.4
    )
    bnd2_row = _annotated_row(
        id="bnd0", chrom="chr2", end=3050, GC_content_left=0.6, GC_content_right=0.8
    )
    bnd1_row.pop("chrom_2")
    bnd2_row.pop("chrom_2")
    obj.df_calls_annot_bnd1 = pd.DataFrame([bnd1_row])
    obj.df_calls_annot_bnd2 = pd.DataFrame([bnd2_row])

    obj.aggregate_results()

    out = obj.df_calls_annot
    # Original DEL plus the single re-aggregated BND row.
    assert set(out["id"]) == {"del0", "bnd0"}
    assert len(out) == 2
    # The merged frame keeps the main frame's column layout.
    assert list(out.columns) == main_cols

    bnd = out.loc[out["id"] == "bnd0"].iloc[0]
    # aggregate_results adds 50 back to bnd1.start before the groupby 'first'.
    assert bnd["start"] == 1950 + 50
    # bnd2.end has 50 subtracted; groupby takes 'last' for end.
    assert bnd["end"] == 3050 - 50
    # chrom_2 is recreated from bnd2's chrom ('chr2') and aggregated with 'last'.
    assert bnd["chrom_2"] == "chr2"
    # GC handling: bnd1 sets left = mean(its left, right) and right = NaN; bnd2
    # sets right = mean(its left, right) and left = NaN. groupby 'first' skips
    # NaN, so the merged row takes left from bnd1 and right from bnd2:
    #   left  = (0.2 + 0.4) / 2 = 0.3   (from bnd1)
    #   right = (0.6 + 0.8) / 2 = 0.7   (from bnd2)
    assert bnd["GC_content_left"] == pytest.approx(0.3)
    assert bnd["GC_content_right"] == pytest.approx(0.7)


# ---------------------------------------------------------------------------
# to_df: returns the annotation frame.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_to_df_returns_annotation_frame():
    obj = _blank_annotator()
    frame = pd.DataFrame({"id": ["a"], "start": [1]})
    obj.df_calls_annot = frame
    assert obj.to_df() is frame


# ---------------------------------------------------------------------------
# INTEGRATION: build a real annotator from synthetic tracks and run the chain.
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_full_annotation_pipeline(tmp_path):
    reference_filenames = sref.build_reference_files(tmp_path, gc_value=0.42)
    annotator = ReferenceAnnotator(reference_filenames=reference_filenames)

    annotator.load_from_df(sref.make_calls_df())
    annotator.split_bnd()
    # No BND in the call frame -> the flag is False, exercising the non-split path.
    assert annotator.split_bnd is False

    annotator.annotate_repeats()
    annotator.annotate_vntrs()
    annotator.annotate_strs()
    annotator.annotate_cpg_islands()
    annotator.annotate_centromeres()
    annotator.annotate_asmb_gaps()
    annotator.annotate_alt_haps()
    annotator.annotate_gc_content()
    annotator.aggregate_results()  # no-op (not split), but covers the branch

    out = annotator.to_df()
    assert len(out) == 1
    row = out.iloc[0]

    # Annotation columns are present.
    for col in ["rep_LINE", "rep_VNTR", "rep_STR", "cpg_islands",
                "centromeres", "asmb_gaps", "alt_haps",
                "GC_content_left", "GC_content_right"]:
        assert col in out.columns

    # The call [1000, 1200) overlaps the LINE repeat / VNTR / STR / CpG island /
    # assembly gap, so bf.closest reports distance 0 for each.
    assert row["rep_LINE"] == 0
    assert row["rep_VNTR"] == 0
    assert row["rep_STR"] == 0
    assert row["cpg_islands"] == 0
    assert row["asmb_gaps"] == 0

    # The centromere (40k-45k) and alt haplotype (50k-55k) are far away -> > 0.
    assert row["centromeres"] > 0
    assert row["alt_haps"] > 0

    # GC content comes from the bigWig value we wrote.
    assert row["GC_content_left"] == pytest.approx(0.42)
    assert row["GC_content_right"] == pytest.approx(0.42)


@pytest.mark.unit
def test_full_annotation_pipeline_with_bnd(tmp_path):
    """A BND call drives the BND branches of every annotate_* method.

    With a BND present ``split_bnd`` is True, so each ``annotate_*`` method also
    annotates ``df_calls_annot_bnd1`` / ``df_calls_annot_bnd2`` (the +/-50
    breakpoint windows). ``aggregate_results`` is *not* called here: its groupby
    references a fixed list of ``rep_*`` classes that the tiny repeat fixture
    does not all provide, so the populated merge is covered by the dedicated
    pure test instead.

    Note: this exercises BND handling that still lives inside
    ``collect_reference.py`` at the unit level (bins, ``split_bnd``, chrom_2
    parsing). The repo-wide removal of BND *scoring* only means the CLI's
    default sv_types list and the model layer never invoke it end-to-end --
    the module-level code paths tested here are still live and reachable.
    """
    reference_filenames = sref.build_reference_files(tmp_path, gc_value=0.42)
    annotator = ReferenceAnnotator(reference_filenames=reference_filenames)

    annotator.load_from_df(sref.make_calls_df_with_bnd())
    annotator.split_bnd()
    assert annotator.split_bnd is True
    # The DEL stays in the main frame; the BND is split out into two halves.
    assert list(annotator.df_calls_annot["sv_type"]) == ["DEL"]
    assert len(annotator.df_calls_annot_bnd1) == 1
    assert len(annotator.df_calls_annot_bnd2) == 1

    annotator.annotate_repeats()
    annotator.annotate_vntrs()
    annotator.annotate_strs()
    annotator.annotate_cpg_islands()
    annotator.annotate_centromeres()
    annotator.annotate_asmb_gaps()
    annotator.annotate_alt_haps()
    annotator.annotate_gc_content()

    # The BND breakpoint window [1100, 1200) overlaps the LINE repeat / VNTR /
    # STR / CpG / assembly-gap tracks, so its first-chromosome half is flagged.
    bnd1 = annotator.df_calls_annot_bnd1.iloc[0]
    assert bnd1["rep_LINE"] == 0
    assert bnd1["rep_VNTR"] == 0
    assert bnd1["cpg_islands"] == 0
    # GC content was annotated onto both BND halves.
    assert annotator.df_calls_annot_bnd1.iloc[0]["GC_content_left"] == pytest.approx(0.42)
    assert annotator.df_calls_annot_bnd2.iloc[0]["GC_content_right"] == pytest.approx(0.42)


@pytest.mark.unit
def test_gc_content_out_of_range_is_nan(tmp_path):
    """A call near the contig end queries past the bigWig -> NaN via except Exception."""
    reference_filenames = sref.build_reference_files(tmp_path)
    annotator = ReferenceAnnotator(reference_filenames=reference_filenames)

    df = sref.make_calls_df()
    # Push 'end' just past the contig length so stats() raises -> NaN.
    df.loc[0, "end"] = sref.CONTIG_LENGTH + 10
    annotator.load_from_df(df)
    annotator.split_bnd()
    annotator.annotate_gc_content()

    row = annotator.to_df().iloc[0]
    assert np.isnan(row["GC_content_right"])
