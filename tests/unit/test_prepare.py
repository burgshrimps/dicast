"""Unit tests for ``dicast/prepare.py``.

Covers the ``VariantPrep`` class. The only heavy thing the constructor does is
``pd.read_csv(chrom_sizes, ...)`` over a FAI / chrom-sizes TSV, so each test
writes a tiny FAI into ``tmp_path`` and passes its path. For methods that do
not touch real VCFs, the variant dataframe that drives them is assembled
in-test and assigned directly to ``instance.df_variants``.

This repo's ``VariantPrep`` constructor no longer takes ``sample``, ``vcfs``,
or ``mode`` -- those were dropped along with cohort VCF/CSV mode (replaced by
the ``multi`` subcommand). ``sample``/``vcfs`` are now set via the separate
``read_vcf(vcfs, sample)`` method, and there is no ``mode`` concept or
``filter_variants_cohort`` method at all any more.

``read_variants`` now calls ``dicast.vcf_input.read_caller_vcf`` (the single
place dicast reads a VCF for variant calls) once per ``(caller, path)`` pair
in ``self.vcfs``, so its tests here write real, minimal VCFs via
``tests/fixtures/synthetic_vcf.write_vcf`` rather than staying out of scope.

Expected values are derived from the inputs we construct, not from the methods'
own output.
"""
from __future__ import annotations

import pandas as pd
import pytest

from dicast import prepare
from dicast.vcf_input import VcfInputError
from tests.fixtures import synthetic_vcf as sv


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_fai(tmp_path, contigs):
    """Write a minimal FAI TSV and return its path as a string.

    ``VariantPrep`` reads it with ``index_col=0`` and column names
    ``['size', 'offset', 'linebases', 'linewidth'] -> the first TSV column is
    the contig name (becomes the index) and the second is ``size``.

    ``contigs`` is a list of ``(name, size)`` pairs. Delegates to the shared
    ``synthetic_vcf.write_fai`` fixture writer so the FAI text format has one
    definition across the test suite.
    """
    return sv.write_fai(tmp_path, "ref.fa.fai", contigs)


def _make_prep(tmp_path, contigs=None, sv_types=None, workdir=None,
               sample="sampleA", vcfs=None, call_read_vcf=True):
    """Construct a VariantPrep instance backed by a tiny FAI in tmp_path.

    This repo's constructor only takes cohort/ref/workdir/technology/chroms/
    chrom_sizes/sv_types; sample and vcfs are supplied afterwards via
    ``read_vcf`` (unless ``call_read_vcf=False``, to exercise the
    pre-``read_vcf`` state).
    """
    if contigs is None:
        contigs = [("chr1", 100000), ("chr2", 50000)]
    if sv_types is None:
        sv_types = ["DEL", "DUP", "INS", "INV"]
    if workdir is None:
        workdir = str(tmp_path)
    if vcfs is None:
        vcfs = [["manta", "/path/manta.vcf"]]
    fai = _write_fai(tmp_path, contigs)
    prep = prepare.VariantPrep(
        cohort="cohortA",
        ref="hg38",
        workdir=workdir,
        technology="ill",
        chroms=["chr1", "chr2"],
        chrom_sizes=fai,
        sv_types=sv_types,
    )
    if call_read_vcf:
        prep.read_vcf(vcfs, sample)
    return prep


# ---------------------------------------------------------------------------
# __init__
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_init_stores_scalars_and_lists(tmp_path):
    prep = _make_prep(tmp_path)
    assert prep.cohort == "cohortA"
    assert prep.ref == "hg38"
    assert prep.technology == "ill"
    assert prep.workdir == str(tmp_path)
    # Lists are stored verbatim.
    assert prep.chroms == ["chr1", "chr2"]
    assert prep.sv_types == ["DEL", "DUP", "INS", "INV"]
    # 'mode' was removed along with cohort VCF/CSV mode.
    assert not hasattr(prep, "mode")


@pytest.mark.unit
def test_init_reads_chrom_sizes_into_dataframe(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000), ("chr2", 50000)])
    assert isinstance(prep.chrom_sizes, pd.DataFrame)
    # Contig name became the index; the 'size' column holds the length.
    assert prep.chrom_sizes.loc["chr1", "size"] == 100000
    assert prep.chrom_sizes.loc["chr2", "size"] == 50000
    assert list(prep.chrom_sizes.columns) == ["size", "offset", "linebases", "linewidth"]


@pytest.mark.unit
def test_init_does_not_set_sample_or_vcfs(tmp_path):
    # sample/vcfs are no longer constructor params -- they only appear once
    # read_vcf() is called.
    prep = _make_prep(tmp_path, call_read_vcf=False)
    assert not hasattr(prep, "sample")
    assert not hasattr(prep, "vcfs")


# ---------------------------------------------------------------------------
# read_vcf
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_read_vcf_stores_sample_and_vcfs(tmp_path):
    vcfs = [["manta", "/path/manta.vcf"], ["delly", "/path/delly.vcf"]]
    prep = _make_prep(tmp_path, sample="sampleB", vcfs=vcfs)
    assert prep.sample == "sampleB"
    assert prep.vcfs == vcfs


@pytest.mark.unit
def test_read_vcf_does_not_parse_anything(tmp_path):
    # read_vcf() only stores the raw inputs; parsing happens in the separate
    # read_variants() method.
    prep = _make_prep(tmp_path)
    assert not hasattr(prep, "df_variants")


# ---------------------------------------------------------------------------
# read_variants
#
# Delegates per-file parsing to dicast.vcf_input.read_caller_vcf, so these
# tests write real (minimal) VCFs via synthetic_vcf.write_vcf rather than
# assembling a dataframe by hand. Every record here uses a symbolic DEL ALT
# with END+SVLEN on a single contig ('chr1'), the simplest record shape that
# survives read_caller_vcf's SV contract -- the branch coverage for *how*
# records are normalized (contig mapping, sample policy, SV type precedence,
# ...) lives in tests/unit/test_vcf_input.py, not here.
# ---------------------------------------------------------------------------

def _del_record(pos: int, rec_id: str) -> dict:
    """One minimal, always-kept DEL record dict for synthetic_vcf.write_vcf."""
    return {
        'chrom': 'chr1', 'pos': pos, 'id': rec_id, 'alt': '<DEL>',
        'info': f'SVTYPE=DEL;END={pos + 500};SVLEN=-500',
    }


@pytest.mark.unit
def test_read_variants_concatenates_files_and_keeps_reports(tmp_path):
    manta_vcf = sv.write_vcf(
        tmp_path, 'manta.vcf', [_del_record(1000, 'DELA'), _del_record(2000, 'DELB')])
    delly_vcf = sv.write_vcf(tmp_path, 'delly.vcf', [_del_record(5000, 'DELC')])

    prep = _make_prep(
        tmp_path, contigs=[('chr1', 100000)],
        vcfs=[['manta', manta_vcf], ['delly', delly_vcf]])
    prep.read_variants()

    # One FileReport per input file, in input order, each reporting the
    # right caller and usable count.
    assert [r.caller for r in prep.reports] == ['manta', 'delly']
    assert [r.usable for r in prep.reports] == [2, 1]
    assert all(r.status == 'OK' for r in prep.reports)

    # Both files' rows are concatenated, and internal ids (f"{caller}:{ordinal}")
    # are unique across files.
    assert len(prep.df_variants) == 3
    assert prep.df_variants['id'].is_unique
    assert set(prep.df_variants['id']) == {'manta:0', 'manta:1', 'delly:0'}
    assert set(prep.df_variants['sv_type']) == {'DEL'}


@pytest.mark.unit
def test_read_variants_raises_when_a_file_has_zero_usable_records(tmp_path):
    # A header-only VCF (no records at all) reads fine but has nothing
    # usable -- read_variants must not silently carry it forward.
    empty_vcf = sv.write_vcf(tmp_path, 'empty.vcf', [])

    prep = _make_prep(tmp_path, contigs=[('chr1', 100000)], vcfs=[['manta', empty_vcf]])
    with pytest.raises(VcfInputError) as excinfo:
        prep.read_variants()

    # The error names the offending file and points at 'dicast check'.
    assert empty_vcf in str(excinfo.value)
    assert 'dicast check' in str(excinfo.value)


# ---------------------------------------------------------------------------
# check_out_of_bounds — non-BND branch
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_check_out_of_bounds_non_bnd_well_inside_is_false(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000)])
    # [1000, 2000] padded by the default 50 stays inside [0, 100000] -> not out of bounds.
    # (loc lookups can return a numpy bool, so compare by value, not identity.)
    result = prep.check_out_of_bounds(
        "DEL", "chr1", None, 1000, 2000, prep.chrom_sizes
    )
    assert bool(result) is False


@pytest.mark.unit
def test_check_out_of_bounds_non_bnd_past_contig_end_is_true(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000)])
    # This repo's default padding is 50 (not 100 like the older dev line).
    # end + padding (99970 + 50 = 100020) > 100000 -> out of bounds.
    result = prep.check_out_of_bounds(
        "DEL", "chr1", None, 99000, 99970, prep.chrom_sizes
    )
    assert bool(result) is True


@pytest.mark.unit
def test_check_out_of_bounds_non_bnd_before_contig_start_is_true(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000)])
    # start - padding (40 - 50 = -10) < 0 -> out of bounds.
    result = prep.check_out_of_bounds(
        "DEL", "chr1", None, 40, 2000, prep.chrom_sizes
    )
    assert bool(result) is True


@pytest.mark.unit
def test_check_out_of_bounds_default_padding_is_50(tmp_path):
    # This repo lowered the default padding from 100 (older dev line) to 50.
    # A gap of exactly 50 from the contig end is exactly on the boundary and
    # must NOT be flagged (only end + padding > size trips it).
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000)])
    exactly_at_boundary = prep.check_out_of_bounds(
        "DEL", "chr1", None, 1000, 99950, prep.chrom_sizes
    )
    assert bool(exactly_at_boundary) is False
    one_past_boundary = prep.check_out_of_bounds(
        "DEL", "chr1", None, 1000, 99951, prep.chrom_sizes
    )
    assert bool(one_past_boundary) is True


# ---------------------------------------------------------------------------
# check_out_of_bounds — BND branch (uses chrom + chrom_2)
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_check_out_of_bounds_bnd_both_inside_is_false(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000), ("chr2", 50000)])
    # start=1000 on chr1, end=2000 on chr2, both padded by the default 50 stay inside.
    result = prep.check_out_of_bounds(
        "BND", "chr1", "chr2", 1000, 2000, prep.chrom_sizes
    )
    assert bool(result) is False


@pytest.mark.unit
def test_check_out_of_bounds_bnd_second_chrom_out_is_true(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000), ("chr2", 50000)])
    # end=49960 on chr2: 49960 + 50 = 50010 > 50000 -> chrom2 out of bounds.
    result = prep.check_out_of_bounds(
        "BND", "chr1", "chr2", 1000, 49960, prep.chrom_sizes
    )
    assert bool(result) is True


@pytest.mark.unit
def test_check_out_of_bounds_bnd_first_chrom_out_is_true(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000), ("chr2", 50000)])
    # start=10 on chr1: 10 - 50 = -40 < 0 -> chrom1 out of bounds.
    result = prep.check_out_of_bounds(
        "BND", "chr1", "chr2", 10, 2000, prep.chrom_sizes
    )
    assert bool(result) is True


# ---------------------------------------------------------------------------
# filter_variants
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_filter_variants_drops_outbounds_and_unsupported(tmp_path):
    prep = _make_prep(
        tmp_path,
        contigs=[("chr1", 100000)],
        sv_types=["DEL", "DUP"],
    )
    # Four rows:
    #  - in-bounds DEL  -> kept
    #  - out-of-bounds DEL (runs past chr1 end even with the default padding=50) -> dropped
    #  - in-bounds INS  -> dropped (INS not in sv_types)
    #  - in-bounds DUP  -> kept
    prep.df_variants = pd.DataFrame(
        {
            "sv_type": ["DEL", "DEL", "INS", "DUP"],
            "chrom": ["chr1", "chr1", "chr1", "chr1"],
            "chrom_2": [None, None, None, None],
            "start": [1000, 99000, 3000, 4000],
            "end": [2000, 99990, 3500, 5000],
        }
    )
    prep.filter_variants()
    survivors = prep.df_variants
    # Only the in-bounds DEL and in-bounds DUP remain.
    assert list(survivors["sv_type"]) == ["DEL", "DUP"]
    assert list(survivors["start"]) == [1000, 4000]
    # Index was reset to a contiguous range and 'outbounds' helper col removed.
    assert list(survivors.index) == [0, 1]
    assert "outbounds" not in survivors.columns


@pytest.mark.unit
def test_filter_variants_coerces_string_coords_to_int(tmp_path):
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000)], sv_types=["DEL"])
    # Coordinates arrive as strings; filter_variants casts them via astype(int).
    prep.df_variants = pd.DataFrame(
        {
            "sv_type": ["DEL"],
            "chrom": ["chr1"],
            "chrom_2": [None],
            "start": ["1000"],
            "end": ["2000"],
        }
    )
    prep.filter_variants()
    assert prep.df_variants["start"].dtype == int
    assert prep.df_variants.loc[0, "start"] == 1000
    assert prep.df_variants.loc[0, "end"] == 2000


@pytest.mark.unit
def test_filter_variants_returns_early_on_empty_frame(tmp_path):
    # An empty df_variants (e.g. every input file was itself empty) has no
    # 'chrom'/'sv_type' columns to filter on -- filter_variants must return
    # immediately rather than running .apply()/column indexing on it.
    prep = _make_prep(tmp_path, contigs=[("chr1", 100000)], sv_types=["DEL"])
    prep.df_variants = pd.DataFrame()
    prep.filter_variants()
    assert prep.df_variants.empty


@pytest.mark.unit
def test_filter_variants_raises_when_nothing_left_after_sv_types(tmp_path):
    prep = _make_prep(
        tmp_path,
        contigs=[("chr1", 100000)],
        sv_types=["DUP"],  # the only row below is a DEL, so nothing matches
    )
    prep.df_variants = pd.DataFrame(
        {
            "sv_type": ["DEL"],
            "chrom": ["chr1"],
            "chrom_2": [None],
            "start": [1000],
            "end": [2000],
        }
    )
    with pytest.raises(VcfInputError, match="DUP"):
        prep.filter_variants()


# ---------------------------------------------------------------------------
# filter_variants_cohort -- REMOVED in this repo
# ---------------------------------------------------------------------------
# The older dev line's cohort VCF/CSV mode (and its filter_variants_cohort()
# method, which dropped full-AC and low-qual rows) was removed entirely from
# this repo's VariantPrep in favor of the 'multi' subcommand's cross-sample
# rescue (dicast/multi.py). dicast/prepare.py here defines no such
# method, so both of the older test_filter_variants_cohort_* tests are
# dropped rather than adapted.


# ---------------------------------------------------------------------------
# get_variant_df
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_get_variant_df_returns_assigned_frame(tmp_path):
    prep = _make_prep(tmp_path)
    df = pd.DataFrame({"sv_type": ["DEL"], "chrom": ["chr1"]})
    prep.df_variants = df
    returned = prep.get_variant_df()
    assert returned is df
    pd.testing.assert_frame_equal(returned, df)


# ---------------------------------------------------------------------------
# save_variants
# ---------------------------------------------------------------------------

@pytest.mark.unit
def test_save_variants_writes_tsv(tmp_path):
    prep = _make_prep(tmp_path, workdir=str(tmp_path))
    (tmp_path / "input").mkdir()
    df = pd.DataFrame(
        {
            "sv_type": ["DEL", "DUP"],
            "chrom": ["chr1", "chr2"],
            "start": [1000, 5000],
            "end": [2000, 6000],
        }
    )
    prep.df_variants = df
    prep.save_variants()

    # save_variants writes under workdir/input/ (see the new workdir tree
    # layout built by dicast.parsing.make_workdir_tree).
    expected_path = tmp_path / "input" / "sampleA_hg38.SVs.raw.tsv"
    assert expected_path.exists()

    # Reload and confirm it round-trips to the same frame (TSV, no index col).
    reloaded = pd.read_csv(expected_path, sep="\t")
    pd.testing.assert_frame_equal(reloaded, df)
