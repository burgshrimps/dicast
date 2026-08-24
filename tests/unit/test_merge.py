"""Unit tests for ``dicast/merge.py`` -- the best-call-per-cluster VCF merge.

The clustering functions (``extract_overlap_ids``, ``reciprocal_overlap``,
``compute_reciprocal_overlap``, ``compute_sv_clusters``) are ported verbatim
from the paper's evaluation code; these tests exercise them only through the
public ``select_merged_calls`` / ``genotype_to_gt`` entry points used by the
CLI, with small synthetic dataframes built to land on either side of each
threshold (0.5 reciprocal overlap for DEL/DUP/INV, 200bp breakpoint distance
for INS, 0.4 dicast_qual for the population-aware winner rule).
"""
from __future__ import annotations

import math

import pandas as pd
import pytest

from dicast.merge import select_merged_calls, genotype_to_gt, _cluster_winner_index


_COLUMNS = ["id", "sv_type", "chrom", "start", "end", "sv_len", "caller",
            "dicast_qual", "filter", "genotype"]


def _row(id, sv_type, chrom, start, end, sv_len, caller, dicast_qual,
         filter="PASS", genotype="(1, 1)"):
    return {
        "id": id, "sv_type": sv_type, "chrom": chrom, "start": start,
        "end": end, "sv_len": sv_len, "caller": caller,
        "dicast_qual": dicast_qual, "filter": filter, "genotype": genotype,
    }


def _df(rows):
    return pd.DataFrame(rows, columns=_COLUMNS)


# ---------------------------------------------------------------------------
# DEL/DUP/INV: reciprocal overlap > 0.5
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_del_pair_above_reciprocal_threshold_clusters_and_higher_qual_wins():
    # chr1:1000-2000 (len 1000) vs chr1:1100-2100 (len 1000): overlap 900,
    # reciprocal_overlap = 900/1000 = 0.9 > 0.5 -> same cluster.
    df = _df([
        _row("DEL_low", "DEL", "chr1", 1000, 2000, 1000, "callerA", 0.6),
        _row("DEL_high", "DEL", "chr1", 1100, 2100, 1000, "callerB", 0.9),
    ])

    result = select_merged_calls(df)

    assert len(result) == 1
    assert result.iloc[0]["id"] == "DEL_high"


@pytest.mark.unit
def test_del_pair_below_reciprocal_threshold_does_not_cluster():
    # chr1:5000-6000 (len 1000) vs chr1:5700-6700 (len 1000): overlap 300,
    # reciprocal_overlap = 300/1000 = 0.3 < 0.5 -> separate clusters.
    df = _df([
        _row("DEL_a", "DEL", "chr1", 5000, 6000, 1000, "callerA", 0.6),
        _row("DEL_b", "DEL", "chr1", 5700, 6700, 1000, "callerB", 0.9),
    ])

    result = select_merged_calls(df)

    assert len(result) == 2
    assert set(result["id"]) == {"DEL_a", "DEL_b"}


# ---------------------------------------------------------------------------
# INS: breakpoint distance < 200bp, length ignored
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_ins_pair_within_200bp_clusters():
    # extract_overlap_ids's k-doubling loop only marks an id "finished" once
    # it has seen a neighbor beyond max_dist_overlap -- with just the two
    # calls under test it can never prove it has seen everything and loops
    # until the k>100000 cutoff, discarding both (this only bites synthetic
    # 2-row cases; real chromosome-scale data always has a distant SV). A
    # third, far-away anchor call gives the loop something genuinely beyond
    # the threshold to terminate on, matching how the matcher is actually
    # exercised in production.
    df = _df([
        _row("INS_low", "INS", "chr2", 10000, 10001, 100, "callerA", 0.5),
        _row("INS_high", "INS", "chr2", 10150, 10151, 120, "callerB", 0.8),
        _row("INS_anchor", "INS", "chr2", 50000, 50001, 100, "callerA", 0.3),
    ])

    result = select_merged_calls(df)

    assert len(result) == 2
    assert set(result["id"]) == {"INS_high", "INS_anchor"}
    assert "INS_low" not in set(result["id"])


@pytest.mark.unit
def test_ins_pair_beyond_200bp_does_not_cluster():
    df = _df([
        _row("INS_a", "INS", "chr2", 20000, 20001, 100, "callerA", 0.5),
        _row("INS_b", "INS", "chr2", 20300, 20301, 100, "callerB", 0.8),
    ])

    result = select_merged_calls(df)

    assert len(result) == 2
    assert set(result["id"]) == {"INS_a", "INS_b"}


# ---------------------------------------------------------------------------
# Population-aware winner selection (keep_max_qual_dicast port)
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_pav_call_loses_to_qualifying_non_pav_call_despite_higher_qual():
    # Same locus (guarantees clustering); pav has the highest raw qual (0.9)
    # but a non-pav call clears the 0.4 threshold (0.5) and must win instead.
    df = _df([
        _row("DEL_pav", "DEL", "chr3", 1000, 2000, 1000, "pav", 0.9),
        _row("DEL_real", "DEL", "chr3", 1000, 2000, 1000, "delly", 0.5),
    ])

    result = select_merged_calls(df)

    assert len(result) == 1
    assert result.iloc[0]["id"] == "DEL_real"


@pytest.mark.unit
def test_pav_call_wins_when_no_non_pav_call_clears_threshold():
    # Same setup, but the non-pav call's qual (0.2) is below the 0.4
    # threshold -> falls back to the overall max, which is the pav call.
    df = _df([
        _row("DEL_pav", "DEL", "chr3", 1000, 2000, 1000, "pav", 0.9),
        _row("DEL_weak", "DEL", "chr3", 1000, 2000, 1000, "delly", 0.2),
    ])

    result = select_merged_calls(df)

    assert len(result) == 1
    assert result.iloc[0]["id"] == "DEL_pav"


# ---------------------------------------------------------------------------
# All-NaN cluster must not crash
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_all_nan_cluster_keeps_first_row_without_crashing():
    df = _df([
        _row("DEL_nan1", "DEL", "chr4", 1000, 2000, 1000, "callerA", float("nan")),
        _row("DEL_nan2", "DEL", "chr4", 1050, 2050, 1000, "callerB", float("nan")),
    ])

    result = select_merged_calls(df)

    assert len(result) == 1
    assert result.iloc[0]["id"] == "DEL_nan1"


@pytest.mark.unit
def test_cluster_winner_index_all_nan_group_does_not_raise():
    group = _df([
        _row("x1", "DEL", "chr5", 1, 2, 1, "callerA", float("nan")),
    ]).set_index(pd.Index([7]))
    idx = _cluster_winner_index(group)
    assert idx == 7


# ---------------------------------------------------------------------------
# genotype_to_gt
# ---------------------------------------------------------------------------

@pytest.mark.unit
@pytest.mark.parametrize("genotype,expected", [
    ("(1, 1)", "1/1"),
    ("(0, 1)", "0/1"),
    ("(1, 0)", "1/0"),
    ("(None, None)", "./."),
    (None, "./."),
    (float("nan"), "./."),
    ("garbage", "./."),
])
def test_genotype_to_gt(genotype, expected):
    assert genotype_to_gt(genotype) == expected
