"""Unit tests for ``dicast_lib/multi.py`` (cross-sample rescue for the ``multi``
subcommand).

The old ``lucid/dicast`` dev line's cohort mode (``dicast_lib/cohort.py``'s
``Cohort`` class: VCF/CSV/PED-driven ``find_overlapping_variants``,
``_create_vcf_record``, ``update_cohort_information``, ``get_missing_variants``)
does not exist in this repo at all -- it was replaced wholesale by the
``multi`` subcommand's cross-sample rescue (``find_rescue_candidates`` here).
None of the old ``test_cohort.py`` tests carry over; this file is new coverage
for the new mechanism.

``find_rescue_candidates`` takes ``{sample: own_variant_df}`` (one dataframe
per sample, as produced by ``VariantPrep`` from that sample's own caller
VCFs) and, per SV type, clusters variants across samples using bioframe's
``bf.closest`` (DEL/DUP/INV matched by reciprocal overlap > 0.5, INS matched
by breakpoint distance < 200bp) plus a networkx connected-components pass.
For any cluster missing one or more samples, the highest-qual member is
copied into each missing sample's rescue set, with ``id`` prefixed
``rescue_<origin sample>_<origin caller>_`` and ``caller`` rewritten to
``rescue:<origin sample>:<origin caller>``.
"""
from __future__ import annotations

import pandas as pd
import pytest

from dicast_lib.multi import find_rescue_candidates


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Column set find_rescue_candidates / bf.closest actually touch. Mirrors the
# subset of VariantPrep.get_variant_df()'s real columns (see
# dicast_lib/utils.py:caller_vcf_to_dataframe) that the module reads.
_VARIANT_COLS = ["sample", "id", "caller", "sv_type", "chrom", "start", "end", "sv_len", "qual"]


def _variant(sample, id, caller, chrom="chr1", start=1000, end=2000,
             sv_type="DEL", sv_len=1000, qual=50.0):
    """Build one variant row (dict) with every column the module touches."""
    return {
        "sample": sample,
        "id": id,
        "caller": caller,
        "sv_type": sv_type,
        "chrom": chrom,
        "start": start,
        "end": end,
        "sv_len": sv_len,
        "qual": qual,
    }


def _df(rows):
    return pd.DataFrame(rows, columns=_VARIANT_COLS)


# ---------------------------------------------------------------------------
# Basic rescue: variant found by only one sample
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_single_sample_variant_rescued_into_others():
    """A DEL only sample A's caller found is rescued into B and C, tagged
    with the rescue caller/id convention; A itself gets no rescue rows."""
    del_a = _variant("A", "DEL_chr1_1000_2000", "manta", start=1000, end=2000)
    own_variant_dfs = {
        "A": _df([del_a]),
        "B": _df([]),
        "C": _df([]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    assert list(rescued.keys()) == ["A", "B", "C"]
    assert rescued["A"].empty

    for sample in ("B", "C"):
        rows = rescued[sample]
        assert len(rows) == 1
        row = rows.iloc[0]
        assert row["sample"] == sample
        assert row["caller"] == "rescue:A:manta"
        assert row["id"] == "rescue_A_manta_DEL_chr1_1000_2000"
        # The transplanted coordinates/sv_type match the source variant.
        assert row["chrom"] == "chr1"
        assert row["start"] == 1000
        assert row["end"] == 2000
        assert row["sv_type"] == "DEL"


@pytest.mark.unit
def test_variant_found_by_all_samples_not_rescued():
    """When every sample's own calls already overlap (reciprocal overlap >
    0.5), the cluster contains every sample -> nothing gets rescued."""
    del_a = _variant("A", "DEL_A", "manta", start=1000, end=2000)
    del_b = _variant("B", "DEL_B", "delly", start=1020, end=2020)  # recip overlap ~0.96
    own_variant_dfs = {
        "A": _df([del_a]),
        "B": _df([del_b]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    assert rescued["A"].empty
    assert rescued["B"].empty


@pytest.mark.unit
def test_symmetric_rescue_across_two_samples():
    """Two independently-called variants (different SV types so they never
    cluster with each other) each get rescued into the sample that lacks
    them -- rescue runs in both directions within a single call."""
    del_a = _variant("A", "DEL_A", "manta", sv_type="DEL", start=1000, end=2000)
    dup_b = _variant("B", "DUP_B", "delly", sv_type="DUP", start=5000, end=6000)
    own_variant_dfs = {
        "A": _df([del_a]),
        "B": _df([dup_b]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    # A lacks the DUP that B found -> A gets a rescue row sourced from B.
    assert len(rescued["A"]) == 1
    a_row = rescued["A"].iloc[0]
    assert a_row["sample"] == "A"
    assert a_row["sv_type"] == "DUP"
    assert a_row["caller"] == "rescue:B:delly"
    assert a_row["id"] == "rescue_B_delly_DUP_B"

    # B lacks the DEL that A found -> B gets a rescue row sourced from A.
    assert len(rescued["B"]) == 1
    b_row = rescued["B"].iloc[0]
    assert b_row["sample"] == "B"
    assert b_row["sv_type"] == "DEL"
    assert b_row["caller"] == "rescue:A:manta"
    assert b_row["id"] == "rescue_A_manta_DEL_A"


# ---------------------------------------------------------------------------
# INS breakpoint-distance matching
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_ins_within_breakpoint_distance_not_rescued():
    """Two INS calls (one per sample) whose start positions are < 200bp
    apart cluster together -> both samples already 'have' it, no rescue."""
    ins_a = _variant("A", "INS_A", "sniffles", sv_type="INS", start=1000, end=1001, sv_len=300)
    ins_b = _variant("B", "INS_B", "cutesv", sv_type="INS", start=1100, end=1101, sv_len=300)  # 100bp away
    own_variant_dfs = {
        "A": _df([ins_a]),
        "B": _df([ins_b]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    assert rescued["A"].empty
    assert rescued["B"].empty


@pytest.mark.unit
def test_ins_beyond_breakpoint_distance_rescued_both_ways():
    """Two INS calls whose start positions are >= 200bp apart do NOT cluster
    -> each is treated as its own (single-sample) event and gets rescued
    into the other sample."""
    ins_a = _variant("A", "INS_A", "sniffles", sv_type="INS", start=1000, end=1001, sv_len=300)
    ins_b = _variant("B", "INS_B", "cutesv", sv_type="INS", start=1300, end=1301, sv_len=300)  # 300bp away
    own_variant_dfs = {
        "A": _df([ins_a]),
        "B": _df([ins_b]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    assert len(rescued["A"]) == 1
    assert rescued["A"].iloc[0]["id"] == "rescue_B_cutesv_INS_B"
    assert len(rescued["B"]) == 1
    assert rescued["B"].iloc[0]["id"] == "rescue_A_sniffles_INS_A"


# ---------------------------------------------------------------------------
# --pop interaction: the PAV population catalog is added as an extra
# per-sample "caller" (dicast.py / dicast_lib/parsing.py wire
# [sample, 'pav', args.pop_catalog] into every sample's VCF list when --pop
# is set), so identical pav-derived calls land in every sample's own
# dataframe under caller == 'pav'. find_rescue_candidates has no special
# casing for that caller name -- it should cluster/skip exactly like any
# other caller.
# ---------------------------------------------------------------------------


@pytest.mark.unit
def test_pop_catalog_entry_present_in_all_samples_not_rescued():
    """A 'pav' catalog entry that landed in every sample's own dataframe at
    the same coordinates (as --pop adds it identically per sample) is
    already present everywhere -> no rescue, regardless of the special
    caller name."""
    pav_a = _variant("A", "pav.chr1.1000", "pav", start=1000, end=2000, qual=100.0)
    pav_b = _variant("B", "pav.chr1.1000", "pav", start=1000, end=2000, qual=100.0)
    own_variant_dfs = {
        "A": _df([pav_a]),
        "B": _df([pav_b]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    assert rescued["A"].empty
    assert rescued["B"].empty


@pytest.mark.unit
def test_pop_catalog_entry_only_matched_by_one_sample_own_caller_rescued():
    """A pav catalog entry overlaps a variant sample A's own caller found
    but sample B lacks entirely -- B still gets rescued via the normal
    mechanism (matching is by position, not by caller identity), sourced
    from whichever cluster member has the higher qual."""
    del_a_own = _variant("A", "DEL_A_manta", "manta", start=1000, end=2000, qual=40.0)
    pav_a = _variant("A", "pav.chr1.1000", "pav", start=1010, end=2010, qual=100.0)  # overlaps del_a_own
    own_variant_dfs = {
        "A": _df([del_a_own, pav_a]),
        "B": _df([]),
    }

    rescued = find_rescue_candidates(own_variant_dfs)

    # A's own two calls cluster into one group (both from sample A already);
    # B lacks the whole cluster and gets rescued the higher-qual member (pav).
    assert rescued["A"].empty
    assert len(rescued["B"]) == 1
    row = rescued["B"].iloc[0]
    assert row["caller"] == "rescue:A:pav"
    assert row["id"] == "rescue_A_pav_pav.chr1.1000"
