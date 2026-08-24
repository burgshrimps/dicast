"""Reusable synthetic VCF builder for the VCF-parser unit tests.

This helper writes a tiny, bgzipped + tabix-indexed VCF into a caller-supplied
directory (typically pytest's ``tmp_path``) and opens it with
:class:`pysam.VariantFile`, so :func:`dicast_lib.utils.caller_vcf_to_dataframe`
can be exercised without any real data.

The header declares every INFO / FORMAT field the parser reads so that field
access never raises ``KeyError``:

* ``caller_vcf_to_dataframe`` reads INFO ``SVTYPE``, ``END``, ``SVLEN`` and
  FORMAT ``GT`` (the latter is read unconditionally for every kept record).

The records are hand-built to collectively cover the parser branches:

* a DEL and an INV on a canonical chrom (with END / SVLEN),
* an INS on a canonical chrom (with SVLEN),
* a DUP on a canonical chrom,
* a BND on a canonical chrom whose ALT carries a mate bracket
  ``N[chrMATE:POS[`` so the regex-based chrom_2 / end extraction fires,
* a DEL on a NON-canonical chrom (dropped by the canonical filter),
* a record with an unsupported SVTYPE (dropped by the sv_types filter),

and the sample column carries a mix of present (``0/1``) and absent (``0/0``)
genotypes so the ``check_genotype`` branch can be exercised.

EXPECTED values used in assertions are derived from what is written here, never
from running the parser first. The two VCF coordinate conventions matter:

* ``rec.start`` is 0-based, so a VCF POS of ``P`` gives ``rec.start == P - 1``
  and the parser stores ``start = rec.start + 1 == P`` (the original POS).
* ``rec.stop`` follows the INFO ``END`` (1-based, inclusive) when present.

NOTE: the older dicast dev line also had a ``sample_vcf_to_dataframe`` parser
(for VCF-based cohort mode, with extra INFO fields CALLER / COHORT_AC /
SUPP_SAMPLES / SUPP_SAMPLES_GT) and a matching ``make_sample_vcf`` builder
here. This repo removed cohort VCF/CSV mode entirely (replaced by the
``multi`` subcommand's cross-sample rescue in ``dicast_lib/multi.py``), and
``dicast_lib/utils.py`` no longer defines ``sample_vcf_to_dataframe`` at all,
so that builder and the INFO fields it alone needed are dropped here.
"""
from __future__ import annotations

import pysam

SAMPLE = "SAMPLE1"

# Canonical chromosomes the parser is told about; ``chrUn`` is deliberately
# omitted so the non-canonical record gets filtered out.
CANONICAL_CHROMS = ["chr1", "chr2", "chr5"]

_HEADER_LINES = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=chr1,length=100000>",
    "##contig=<ID=chr2,length=100000>",
    "##contig=<ID=chr5,length=100000>",
    "##contig=<ID=chrUn,length=100000>",
    '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of SV">',
    '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of SV">',
    '##ALT=<ID=DEL,Description="Deletion">',
    '##ALT=<ID=INS,Description="Insertion">',
    '##ALT=<ID=INV,Description="Inversion">',
    '##ALT=<ID=DUP,Description="Duplication">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
]


# Each record is a dict describing how it is written; the EXPECTED values in the
# tests are read straight off these dicts, so there is a single source of truth.
#
# Fields: chrom, pos (1-based VCF POS), ref_id, alt, svtype, end, svlen, gt.
RECORDS = [
    # 0) DEL on canonical chr1, present in sample (0/1).
    {
        "chrom": "chr1", "pos": 1000, "id": "DEL1", "ref": "N", "alt": "<DEL>",
        "svtype": "DEL", "end": 1500, "svlen": -500, "gt": "0/1",
    },
    # 1) INS on canonical chr2, absent from sample (0/0).
    {
        "chrom": "chr2", "pos": 2000, "id": "INS1", "ref": "N", "alt": "<INS>",
        "svtype": "INS", "end": None, "svlen": 250, "gt": "0/0",
    },
    # 2) INV on canonical chr5, present (1/1).
    {
        "chrom": "chr5", "pos": 3000, "id": "INV1", "ref": "N", "alt": "<INV>",
        "svtype": "INV", "end": 3800, "svlen": 800, "gt": "1/1",
    },
    # 3) DUP on canonical chr1, present (0/1).
    {
        "chrom": "chr1", "pos": 4000, "id": "DUP1", "ref": "N", "alt": "<DUP>",
        "svtype": "DUP", "end": 4600, "svlen": 600, "gt": "0/1",
    },
    # 4) BND on canonical chr1 whose ALT carries a mate bracket to chr5:9000.
    #    The regex pulls chrom_2='chr5' and end='9000' from the ALT string.
    {
        "chrom": "chr1", "pos": 5000, "id": "BND1", "ref": "N",
        "alt": "N[chr5:9000[",
        "svtype": "BND", "end": None, "svlen": None, "gt": "0/1",
    },
    # 5) DEL on NON-canonical chrUn -> dropped by the canonical-chrom filter.
    {
        "chrom": "chrUn", "pos": 6000, "id": "DELUN", "ref": "N", "alt": "<DEL>",
        "svtype": "DEL", "end": 6500, "svlen": -500, "gt": "0/1",
    },
    # 6) Unsupported SVTYPE (CNV) on canonical chr2 -> dropped by sv_types filter.
    {
        "chrom": "chr2", "pos": 7000, "id": "CNV1", "ref": "N", "alt": "<CNV>",
        "svtype": "CNV", "end": 7500, "svlen": 500, "gt": "0/1",
    },
]

# IDs of records that survive both filters (canonical chrom + supported
# svtype). Used to derive the expected row count (order-independent).
KEPT_IDS = ["DEL1", "INS1", "INV1", "DUP1", "BND1"]
# Of the kept records, those whose sample GT contains a '1' allele (i.e. the
# rows that survive check_genotype=True). INS1 is 0/0 so it drops out.
KEPT_PRESENT_IDS = ["DEL1", "INV1", "DUP1", "BND1"]


def record_by_id(rec_id):
    """Return the source RECORDS dict for ``rec_id`` (single source of truth)."""
    for rec in RECORDS:
        if rec["id"] == rec_id:
            return rec
    raise KeyError(rec_id)


def _info_field(rec):
    """Build the INFO column string for one record dict."""
    parts = [f"SVTYPE={rec['svtype']}"]
    if rec["end"] is not None:
        parts.append(f"END={rec['end']}")
    if rec["svlen"] is not None:
        parts.append(f"SVLEN={rec['svlen']}")
    return ";".join(parts)


# Contig order as declared in the header; tabix needs records grouped by this
# order and ascending POS within each contig.
_CONTIG_ORDER = ["chr1", "chr2", "chr5", "chrUn"]


def make_caller_vcf(tmp_path, name="caller.vcf"):
    """Write, bgzip + tabix-index, and open a synthetic VCF for
    :func:`dicast_lib.utils.caller_vcf_to_dataframe`.

    Returns an open :class:`pysam.VariantFile`.
    """
    vcf_path = tmp_path / name
    lines = list(_HEADER_LINES)
    lines.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + SAMPLE
    )
    # Sort by (contig order, POS) so tabix can index; the parser fetches in this
    # coordinate order, not in RECORDS declaration order.
    sorted_records = sorted(
        RECORDS, key=lambda r: (_CONTIG_ORDER.index(r["chrom"]), r["pos"])
    )
    for rec in sorted_records:
        info = _info_field(rec)
        lines.append(
            "\t".join(
                [
                    rec["chrom"],
                    str(rec["pos"]),
                    rec["id"],
                    rec["ref"],
                    rec["alt"],
                    "50",          # QUAL
                    "PASS",        # FILTER
                    info,
                    "GT",          # FORMAT
                    rec["gt"],     # sample genotype
                ]
            )
        )
    vcf_path.write_text("\n".join(lines) + "\n")
    gz = pysam.tabix_index(str(vcf_path), preset="vcf", force=True)
    return pysam.VariantFile(gz)
