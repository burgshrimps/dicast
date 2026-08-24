"""Reusable synthetic reference-track builders for ``ReferenceAnnotator`` tests.

:class:`dicast.collect_reference.ReferenceAnnotator` eagerly reads seven TSV
tracks (repeats, VNTRs, STRs, CpG islands, centromeres, assembly gaps,
alternative haplotypes) and opens one bigWig (GC content) in its constructor.
The helpers here write tiny versions of each into a caller-supplied directory
(typically pytest's ``tmp_path``) and return a ``reference_filenames`` dict with
exactly the eight keys the constructor expects, so a *real* annotator can be
built and driven end-to-end without any production reference data.

The tracks are hand-built so that the single SV call used by the integration
test (a DEL on ``chr1`` spanning ``[1000, 1200)``) overlaps a repeat, a VNTR,
an STR, a CpG island and an assembly gap, and sits a known distance from a
centromere and an alternative haplotype. ``bf.closest(..., return_overlap=True)``
returns ``distance == 0`` for an overlap, so the expected annotation values in
the assertions are derived directly from these coordinates, never from running
the annotator first.

Column layouts mirror what each ``annotate_*`` method reads from its frame:

* repeats   -- ``genoName / genoStart / genoEnd / repClass`` (index col first;
               written via ``index=True`` so ``index_col=0`` round-trips)
* vntrs     -- headerless ``chrom / start / stop / class``
* strs      -- headerless ``chrom / start / stop / len_unit / seq_unit / unknown``
* cpgislands/centromeres/asmb_gaps -- ``chrom / chromStart / chromEnd``
* alt_haps  -- ``tName / tStart / tEnd``
* gc        -- bigWig with one ``chr1`` interval of known value
"""
from __future__ import annotations

import pandas as pd
import pyBigWig

CONTIG = "chr1"
CONTIG_LENGTH = 100_000

# The single call the integration test annotates: chr1 DEL spanning [1000, 1200).
CALL_CHROM = CONTIG
CALL_START = 1000
CALL_END = 1200


def _write_repeats(path) -> None:
    """Repeatmasker-style track read with ``index_col=0``.

    ``annotate_repeats`` keys on ``repClass`` and the ``genoName/genoStart/
    genoEnd`` interval. The first interval (a LINE) overlaps the call window so
    the call gets ``rep_LINE == 0``; the ``Unknown`` row is dropped by the
    method's repClass filter, exercising that branch.
    """
    df = pd.DataFrame(
        {
            "genoName": [CONTIG, CONTIG, CONTIG],
            "genoStart": [1050, 5000, 1100],
            "genoEnd": [1150, 5100, 1120],
            "repClass": ["LINE", "SINE", "Unknown"],
        }
    )
    # index_col=0 on read consumes the first written column; an explicit index
    # column gives a stable round-trip.
    df.index.name = "bin"
    df.to_csv(path, sep="\t", index=True)


def _write_vntrs(path) -> None:
    """Headerless VNTR track: chrom / start / stop / class."""
    rows = [
        (CONTIG, 1180, 1300, "VNTR_A"),  # overlaps the call -> distance 0
        (CONTIG, 8000, 8100, "VNTR_B"),
    ]
    with open(path, "w") as fh:
        for c, s, e, cls in rows:
            fh.write(f"{c}\t{s}\t{e}\t{cls}\n")


def _write_strs(path) -> None:
    """Headerless STR track: chrom / start / stop / len_unit / seq_unit / unknown."""
    rows = [
        (CONTIG, 1190, 1250, 2, "AT", 0),  # overlaps the call -> distance 0
        (CONTIG, 9000, 9050, 3, "CAG", 0),
    ]
    with open(path, "w") as fh:
        for c, s, e, lu, su, u in rows:
            fh.write(f"{c}\t{s}\t{e}\t{lu}\t{su}\t{u}\n")


def _write_chrom_start_end_track(path, rows) -> None:
    """Generic ``chrom / chromStart / chromEnd`` UCSC-style track.

    Used for CpG islands, centromeres and assembly gaps.
    """
    df = pd.DataFrame(rows, columns=["chrom", "chromStart", "chromEnd"])
    df.to_csv(path, sep="\t", index=False)


def _write_alt_haps(path, rows) -> None:
    """``tName / tStart / tEnd`` alt-haplotype track."""
    df = pd.DataFrame(rows, columns=["tName", "tStart", "tEnd"])
    df.to_csv(path, sep="\t", index=False)


def _write_gc_bigwig(path, value: float = 0.42) -> None:
    """Write a bigWig with a single ``chr1`` interval carrying ``value``.

    ``calculate_gc_content`` queries ``bw.stats(chrom, pos-50, pos+50)``. The
    interval below covers [0, CONTIG_LENGTH) so any in-range query returns
    ``value``.
    """
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader([(CONTIG, CONTIG_LENGTH)])
    bw.addEntries([CONTIG], [0], ends=[CONTIG_LENGTH], values=[float(value)])
    bw.close()


def build_reference_files(tmp_path, *, gc_value: float = 0.42) -> dict:
    """Write all eight reference tracks under ``tmp_path``.

    Returns the ``reference_filenames`` dict with exactly the keys
    :meth:`ReferenceAnnotator.__init__` reads. The default coordinates make the
    integration call overlap repeat/VNTR/STR/CpG/assembly-gap tracks (distance
    0) and sit at a known distance from the centromere and alt-haplotype tracks.
    """
    repeats = tmp_path / "repeats.tsv"
    vntrs = tmp_path / "vntrs.tsv"
    strs = tmp_path / "strs.tsv"
    gc = tmp_path / "gc.bw"
    cpg = tmp_path / "cpgislands.tsv"
    centro = tmp_path / "centromeres.tsv"
    gaps = tmp_path / "asmb_gaps.tsv"
    alt = tmp_path / "alt_haps.tsv"

    _write_repeats(repeats)
    _write_vntrs(vntrs)
    _write_strs(strs)
    _write_gc_bigwig(gc, value=gc_value)

    # CpG island overlaps the call -> distance 0.
    _write_chrom_start_end_track(
        cpg, [(CONTIG, 1150, 1250), (CONTIG, 7000, 7100)]
    )
    # Centromere is far away -> a positive distance.
    _write_chrom_start_end_track(
        centro, [(CONTIG, 40_000, 45_000)]
    )
    # Assembly gap overlaps the call -> distance 0.
    _write_chrom_start_end_track(
        gaps, [(CONTIG, 900, 1100), (CONTIG, 6000, 6100)]
    )
    # Alt haplotype is far away -> a positive distance.
    _write_alt_haps(alt, [(CONTIG, 50_000, 55_000)])

    return {
        "repeats_filename": str(repeats),
        "vntrs_filename": str(vntrs),
        "strs_filename": str(strs),
        "gc_filename": str(gc),
        "cpgislands_filename": str(cpg),
        "centromeres_filename": str(centro),
        "asmb_gaps_filename": str(gaps),
        "alt_haps_filename": str(alt),
    }


def make_calls_df():
    """Return a tiny ``df_calls``-shaped frame with the 11 columns load_from_df reads.

    A single non-BND DEL call on ``chr1`` spanning the window the tracks above
    were laid out around. ``chrom_2`` is present (required by the loader and the
    BND split) but unused for a DEL.
    """
    return pd.DataFrame(
        {
            "id": ["call0"],
            "cohort": ["cohortA"],
            "sample": ["sampleA"],
            "reference": ["GRCh38"],
            "technology": ["ONT"],
            "caller": ["sniffles"],
            "sv_type": ["DEL"],
            "chrom": [CALL_CHROM],
            "chrom_2": [CALL_CHROM],
            "start": [CALL_START],
            "end": [CALL_END],
        }
    )


def make_calls_df_with_bnd():
    """Return a ``df_calls``-shaped frame containing both a DEL and a BND.

    The BND triggers ``split_bnd`` (so the ``annotate_*`` methods run their BND
    branches) and the populated ``aggregate_results`` path. The BND breakpoints
    are placed inside the same window the tracks were laid out around so the
    +/-50 split windows still overlap the repeat/VNTR/STR/CpG/gap intervals.
    """
    return pd.DataFrame(
        {
            "id": ["call0", "bnd0"],
            "cohort": ["cohortA", "cohortA"],
            "sample": ["sampleA", "sampleA"],
            "reference": ["GRCh38", "GRCh38"],
            "technology": ["ONT", "ONT"],
            "caller": ["sniffles", "sniffles"],
            "sv_type": ["DEL", "BND"],
            "chrom": [CALL_CHROM, CALL_CHROM],
            "chrom_2": [CALL_CHROM, CALL_CHROM],
            "start": [CALL_START, 1150],
            "end": [CALL_END, 1180],
        }
    )
