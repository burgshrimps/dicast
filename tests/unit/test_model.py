"""Unit tests for ``dicast_lib/model.py`` (the ``Dicast`` class).

These cover the deterministic, non-xgboost surface of ``Dicast``:

* ``__init__`` feature assembly for every SV type (the dominant ~94-line
  lever). With ``feature_config_file=None`` (the default) no file is read --
  the constructor just builds the hardcoded per-SV-type feature dicts and
  expands bin suffixes into ``self.features``.
* The rule-based scorers ``score_inversions`` / ``score_translocations``
  (no fitted model needed). Note: this repo's CLI never invokes BND scoring
  end-to-end (no BND model ships in models/, and the default --sv-types list
  excludes BND), but the ``Dicast`` class itself is unchanged and still
  supports BND at the unit level (features_aln_bp['BND'], score_translocations),
  so those branches are still exercised here.
* ``impute_missing_values`` (NaN fill behaviour, DEL vs INS branches).
* ``load_from_df`` / ``to_db`` / ``to_df`` (column/shape contracts).
* ``get_feature_importance`` (with a tiny stub model).

The train/predict/save/load methods are xgboost-heavy and intentionally out
of scope here.

``dicast_lib.model`` imports cleanly (given xgboost + pyyaml are installed),
so it is imported directly.
"""
from __future__ import annotations

import types

import numpy as np
import pandas as pd
import pytest

from dicast_lib.model import Dicast


# ---------------------------------------------------------------------------
# Expected feature counts, derived purely from the construction logic in
# __init__ (sv_len + 23 reference features + the expanded alignment features).
#
#   bp features (12) x suffices_bp
#   body features  x suffices_body
#   conn features  x (ordered pairs of bins)
#
# INS uses suffices_bp = [I, II] and has NO body/conn features, so it lands on
# a different suffix branch than the other SV types -- both are covered.
# ---------------------------------------------------------------------------

EXPECTED_TOTAL = {
    "DEL": 104,
    "INS": 48,
    "DUP": 104,
    "INV": 104,
    "BND": 96,
}
EXPECTED_ALN = {
    "DEL": 80,
    "INS": 24,
    "DUP": 80,
    "INV": 80,
    "BND": 72,
}

N_REF = 23  # length of self.features_ref


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("sv_type", ["DEL", "INS", "DUP", "INV", "BND"])
def test_init_sets_basic_attributes(sv_type):
    d = Dicast(sv_type)
    assert d.sv_type == sv_type
    # Log2 coverage threshold default for the WGS model. This repo raised the
    # threshold from 3 (lucid dev line) to 5.
    assert d.cov_thr == 5
    assert d.features_var == ["sv_len"]
    assert len(d.features_ref) == N_REF


@pytest.mark.unit
@pytest.mark.parametrize("sv_type", ["DEL", "INS", "DUP", "INV", "BND"])
def test_init_feature_counts(sv_type):
    d = Dicast(sv_type)
    assert len(d.features_aln) == EXPECTED_ALN[sv_type]
    assert len(d.features) == EXPECTED_TOTAL[sv_type]
    # features = var + ref + aln, in that order.
    assert d.features == d.features_var + d.features_ref + d.features_aln


@pytest.mark.unit
def test_init_default_reads_no_file(tmp_path):
    # feature_config_file defaults to None -> no disk access, no feature_config.
    d = Dicast("DEL")
    assert not hasattr(d, "feature_config")


@pytest.mark.unit
def test_init_del_alignment_suffix_branch():
    # DEL uses the 4-bin breakpoint suffices (I..IV), has body features
    # expanded over the 'a'/'b' body suffices, and connection features
    # expanded over ordered bin pairs.
    d = Dicast("DEL")
    feats = set(d.features_aln)
    # Breakpoint features over all four bins.
    assert "ill_cov_mean_I" in feats
    assert "ill_cov_mean_IV" in feats
    assert "ill_disco_tx_IV" in feats
    # Body features only exist for the cov features over the body suffices.
    assert "ill_cov_mean_IIa" in feats
    assert "ill_cov_std_IIIa" in feats
    # Connection features over ordered bin pairs.
    assert "ill_splitreads_I_II" in feats
    assert "ill_splitreads_III_IV" in feats
    # 'IV' has no outgoing connections.
    assert "ill_splitreads_IV_I" not in feats


@pytest.mark.unit
def test_init_ins_alignment_suffix_branch():
    # INS is the special branch: only two breakpoint bins (I, II), and NO body
    # or connection features at all.
    d = Dicast("INS")
    feats = set(d.features_aln)
    assert "ill_cov_mean_I" in feats
    assert "ill_cov_mean_II" in feats
    # No bins III/IV for INS.
    assert "ill_cov_mean_III" not in feats
    assert "ill_cov_mean_IV" not in feats
    # No body features for INS.
    assert "ill_cov_mean_IIa" not in feats
    # No connection features for INS.
    assert not any("_I_II" in f for f in feats)
    # All alignment features are pure breakpoint features -> 12 bp x 2 bins.
    assert len(d.features_aln) == 24


# ---------------------------------------------------------------------------
# load_from_df
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_load_from_df_filters_by_sv_type_and_resets_index():
    df = pd.DataFrame(
        {
            "sv_type": ["DEL", "INS", "DEL", "DUP"],
            "chrom": ["chr1", "chr1", "chr2", "chr3"],
            "sv_len": [100, 200, 300, 400],
        }
    )
    d = Dicast("DEL")
    d.load_from_df(df)
    # Only the two DEL rows survive.
    assert list(d.variants["sv_len"]) == [100, 300]
    # Index reset to a contiguous range.
    assert list(d.variants.index) == [0, 1]


@pytest.mark.unit
def test_load_from_df_does_not_mutate_input():
    df = pd.DataFrame({"sv_type": ["DEL", "INS"], "sv_len": [10, 20]})
    before = df.copy(deep=True)
    Dicast("DEL").load_from_df(df)
    pd.testing.assert_frame_equal(df, before)


# ---------------------------------------------------------------------------
# impute_missing_values
# ---------------------------------------------------------------------------

def _impute_base_df(sv_type, n=2):
    """A minimal df carrying every column impute_missing_values touches."""
    cols = {
        "sv_type": [sv_type] * n,
        "GC_content_left": [np.nan, 0.4],
        "GC_content_right": [0.6, np.nan],
        "ill_cov_mean_I": [np.nan, 9.0],
        "ill_cov_mean_II": [np.nan, 8.0],
        "ill_cov_std_I": [np.nan, 2.0],
        "ill_cov_std_II": [np.nan, 3.0],
    }
    if sv_type != "INS":
        cols.update(
            {
                "ill_cov_mean_III": [np.nan, 7.0],
                "ill_cov_mean_IV": [np.nan, 6.0],
                "ill_cov_std_III": [np.nan, 4.0],
                "ill_cov_std_IV": [np.nan, 5.0],
            }
        )
    return pd.DataFrame(cols)


@pytest.mark.unit
def test_impute_missing_values_del_branch():
    d = Dicast("DEL")
    d.load_from_df(_impute_base_df("DEL"))
    d.impute_missing_values()
    v = d.variants
    # GC content NaNs filled with the column median.
    # left: [nan, 0.4] -> median 0.4 ; right: [0.6, nan] -> median 0.6
    assert v.loc[0, "GC_content_left"] == pytest.approx(0.4)
    assert v.loc[1, "GC_content_right"] == pytest.approx(0.6)
    # Coverage means NaN -> cov_thr (5 in this repo); coverage std NaN -> 1.
    assert v.loc[0, "ill_cov_mean_I"] == d.cov_thr
    assert v.loc[0, "ill_cov_mean_II"] == d.cov_thr
    assert v.loc[0, "ill_cov_mean_III"] == d.cov_thr
    assert v.loc[0, "ill_cov_mean_IV"] == d.cov_thr
    assert v.loc[0, "ill_cov_std_I"] == 1
    assert v.loc[0, "ill_cov_std_IV"] == 1
    # Non-NaN values are left untouched.
    assert v.loc[1, "ill_cov_mean_III"] == 7.0
    assert v.loc[1, "ill_cov_std_IV"] == 5.0


@pytest.mark.unit
def test_impute_missing_values_ins_branch_skips_iii_iv():
    # The INS branch must NOT touch (or require) bins III/IV.
    d = Dicast("INS")
    d.load_from_df(_impute_base_df("INS"))
    d.impute_missing_values()
    v = d.variants
    assert v.loc[0, "ill_cov_mean_I"] == d.cov_thr
    assert v.loc[0, "ill_cov_std_II"] == 1
    # Bins III/IV columns were never added for INS.
    assert "ill_cov_mean_III" not in v.columns
    assert "ill_cov_mean_IV" not in v.columns


# ---------------------------------------------------------------------------
# score_inversions  (rule-based, no model)
# ---------------------------------------------------------------------------

def _inversion_df():
    """Two rows: row 0 satisfies every inversion mask, row 1 fails them.

    Inversion is scored 1 when ALL of:
      * >=3 of the four clipreads bins  > 0.2
      * >=3 of the eight disco ff/rr cols > 0.1  (relaxed)
      * >=2 of the eight disco ff/rr cols > 0.2  (strict)
      * all four cov_mean bins <= 3.5
      * sv_len < 3_000_000
    """
    bp = ["I", "II", "III", "IV"]
    row_hit = {}
    row_miss = {}
    for b in bp:
        # clipreads: hit has all 4 > 0.2; miss has all 0.
        row_hit[f"ill_clipreads_{b}"] = 0.9
        row_miss[f"ill_clipreads_{b}"] = 0.0
        # disco ff/rr: hit has all > 0.2; miss has all 0.
        row_hit[f"ill_disco_ff_{b}"] = 0.5
        row_miss[f"ill_disco_ff_{b}"] = 0.0
        row_hit[f"ill_disco_rr_{b}"] = 0.5
        row_miss[f"ill_disco_rr_{b}"] = 0.0
        # cov_mean: hit <= 3.5; miss high.
        row_hit[f"ill_cov_mean_{b}"] = 2.0
        row_miss[f"ill_cov_mean_{b}"] = 10.0
    row_hit["sv_len"] = 1000
    row_miss["sv_len"] = 1000
    row_hit["sv_type"] = "INV"
    row_miss["sv_type"] = "INV"
    row_hit["chrom"] = "chr1"
    row_miss["chrom"] = "chr2"
    return pd.DataFrame([row_hit, row_miss])


@pytest.mark.unit
def test_score_inversions_assigns_quality():
    d = Dicast("INV")
    d.load_from_df(_inversion_df())
    d.score_inversions()
    quals = list(d.variants_predict["dicast_qual"])
    assert quals[0] == 1  # all masks satisfied
    assert quals[1] == 0  # fails clip/disco/cov masks


@pytest.mark.unit
def test_score_inversions_long_sv_not_scored():
    df = _inversion_df()
    # Make the otherwise-passing row exceed the length cap.
    df.loc[0, "sv_len"] = 5_000_000
    d = Dicast("INV")
    d.load_from_df(df)
    d.score_inversions()
    assert d.variants_predict.loc[0, "dicast_qual"] == 0


@pytest.mark.unit
def test_score_inversions_chrom_subset():
    d = Dicast("INV")
    d.load_from_df(_inversion_df())
    d.score_inversions(chroms=["chr1"])
    # Only the chr1 row (the hit) is kept.
    assert d.variants_predict.shape[0] == 1
    assert d.variants_predict.loc[0, "dicast_qual"] == 1


@pytest.mark.unit
def test_score_inversions_empty_subset_sets_nan():
    d = Dicast("INV")
    d.load_from_df(_inversion_df())
    d.score_inversions(chroms=["chrZ"])  # matches nothing
    assert d.variants_predict.shape[0] == 0
    assert "dicast_qual" in d.variants_predict.columns


# ---------------------------------------------------------------------------
# score_translocations  (rule-based, no model)
#
# The CLI's default --sv-types list excludes BND and no BND model ships in
# models/, so this scorer is unreachable end-to-end via `dicast.py`. But the
# Dicast class itself still defines score_translocations and the BND branch
# of features_aln_bp unchanged from the lucid dev line, so it remains valid
# to exercise at the unit level.
# ---------------------------------------------------------------------------

def _translocation_df():
    """Row 0 satisfies the translocation rule, row 1 fails it.

    Scored 1 when ALL of:
      * >=2 of four clipreads bins  > 0.2
      * all four cov_mean bins <= 3
      * (>=2 ff/rr > 0.2)  OR  (>=2 rf > 0.3)  OR  (>=2 tx > 0.3)
      * all four mapq_mean bins >= -0.5
      * >=2 of four splitreads bins > 0.1
    """
    bp = ["I", "II", "III", "IV"]
    row_hit = {}
    row_miss = {}
    for b in bp:
        row_hit[f"ill_clipreads_{b}"] = 0.9
        row_miss[f"ill_clipreads_{b}"] = 0.0
        row_hit[f"ill_cov_mean_{b}"] = 1.0
        row_miss[f"ill_cov_mean_{b}"] = 10.0
        # Use the tx (translocation) disco channel for the OR clause.
        row_hit[f"ill_disco_tx_{b}"] = 0.9
        row_miss[f"ill_disco_tx_{b}"] = 0.0
        row_hit[f"ill_disco_ff_{b}"] = 0.0
        row_miss[f"ill_disco_ff_{b}"] = 0.0
        row_hit[f"ill_disco_rr_{b}"] = 0.0
        row_miss[f"ill_disco_rr_{b}"] = 0.0
        row_hit[f"ill_disco_rf_{b}"] = 0.0
        row_miss[f"ill_disco_rf_{b}"] = 0.0
        row_hit[f"ill_mapq_mean_{b}"] = 1.0
        row_miss[f"ill_mapq_mean_{b}"] = -5.0
        row_hit[f"ill_splitreads_{b}"] = 0.9
        row_miss[f"ill_splitreads_{b}"] = 0.0
    row_hit["sv_type"] = "BND"
    row_miss["sv_type"] = "BND"
    row_hit["chrom"] = "chr1"
    row_miss["chrom"] = "chr2"
    return pd.DataFrame([row_hit, row_miss])


@pytest.mark.unit
def test_score_translocations_assigns_quality():
    d = Dicast("BND")
    d.load_from_df(_translocation_df())
    d.score_translocations()
    quals = list(d.variants_predict["dicast_qual"])
    assert quals[0] == 1
    assert quals[1] == 0


@pytest.mark.unit
def test_score_translocations_high_cov_not_scored():
    df = _translocation_df()
    # Bump coverage above the <=3 cap on the otherwise-passing row.
    df.loc[0, "ill_cov_mean_I"] = 50.0
    d = Dicast("BND")
    d.load_from_df(df)
    d.score_translocations()
    assert d.variants_predict.loc[0, "dicast_qual"] == 0


@pytest.mark.unit
def test_score_translocations_empty_subset():
    d = Dicast("BND")
    d.load_from_df(_translocation_df())
    d.score_translocations(chroms=["chrZ"])
    assert d.variants_predict.shape[0] == 0
    assert "dicast_qual" in d.variants_predict.columns


# ---------------------------------------------------------------------------
# to_db / to_df
# ---------------------------------------------------------------------------

_DB_COLUMNS = [
    "single_id", "merged_id", "caller_id", "cohort", "sample", "reference",
    "technology", "caller", "sv_type", "chrom", "chrom_2", "start", "end",
    "sv_len", "filter", "caller_qual", "dicast_qual", "genotype",
    "performed_confirmation", "confirmation_status", "performed_curation",
    "curation_status",
]

_DF_COLUMNS = [
    "id", "cohort", "sample", "reference", "technology", "caller", "sv_type",
    "chrom", "chrom_2", "start", "end", "sv_len", "filter", "qual",
    "dicast_qual", "genotype",
]


def _export_df(columns, n=3):
    return pd.DataFrame({c: list(range(n)) for c in columns})


@pytest.mark.unit
def test_to_db_returns_export_columns():
    d = Dicast("DEL")
    d.variants_predict = _export_df(_DB_COLUMNS)
    out = d.to_db()
    assert list(out.columns) == _DB_COLUMNS
    assert out.shape == (3, len(_DB_COLUMNS))


@pytest.mark.unit
def test_to_df_returns_export_columns():
    d = Dicast("DEL")
    d.variants_predict = _export_df(_DF_COLUMNS)
    out = d.to_df()
    assert list(out.columns) == _DF_COLUMNS
    assert out.shape == (3, len(_DF_COLUMNS))


# ---------------------------------------------------------------------------
# get_feature_importance  (stub model -- no xgboost fitting)
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_feature_importance_maps_and_sorts():
    d = Dicast("DEL")
    n = len(d.features)
    # Ascending importances so the highest is the LAST feature; after sorting
    # descending, that feature must come first.
    importances = np.arange(n, dtype=float)
    d.model = types.SimpleNamespace(feature_importances_=importances)

    out = d.get_feature_importance()

    assert list(out.columns) == ["feature", "importance"]
    assert len(out) == n
    # Sorted descending by importance.
    assert list(out["importance"]) == sorted(importances, reverse=True)
    # The top row maps to the feature that had the largest importance.
    assert out.iloc[0]["feature"] == d.features[-1]
    assert out.iloc[0]["importance"] == importances[-1]
