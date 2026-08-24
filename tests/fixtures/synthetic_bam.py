"""Reusable synthetic BAM builder for the BAM-dependent unit tests of
:mod:`dicast_lib.collect_illumina`.

This helper writes a tiny, coordinate-sorted, indexed BAM into a caller-supplied
directory (typically pytest's ``tmp_path``) so that
:class:`dicast_lib.collect_illumina.AlignmentAnnotatorIllumina` can be exercised
end-to-end without any real data.

The reads are hand-built to collectively cover the branches the read-based
feature extractor cares about:

* plain matches, soft-clips (op 4), an insertion (op 1), a deletion (op 2)
* forward and reverse strands
* a split read carrying an ``SA`` tag whose supplementary alignment lands inside
  one of the breakpoint bins
* paired mates with FF / RR / RF orientation combinations and assorted
  ``template_length`` signs/magnitudes
* a couple of full-window-covering reads so the coverage signal is positive

The single contig is named ``chr1`` so that ``AlignmentAnnotatorIllumina``'s
baseline samplers (which parse ``int(chrom[3:]) - 1`` -> 0) index it correctly.
The contig is short (5000 bp) so the random baseline sampling is fast and, with
the dense full-length reads added below, lands on a positive median depth.
"""
from __future__ import annotations

import pysam

CONTIG = "chr1"
CONTIG_LENGTH = 5000

# The SV / window the integration test annotates. A DEL spanning [1000, 1400).
SV_START = 1000
SV_END = 1400


def _header(contig: str = CONTIG, length: int = CONTIG_LENGTH) -> pysam.AlignmentHeader:
    return pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": contig, "LN": length}],
        }
    )


def _seq_and_qual(n: int):
    """A trivially valid SEQ/QUAL pair of length ``n`` (all 'A', Phred 30)."""
    return "A" * n, [30] * n


def _new_read(header, name, ref_start, cigar, *, is_reverse=False, mapq=60,
              flag_paired=False, mate_unmapped=True, mate_reverse=False,
              is_read1=True, next_ref_id=0, next_ref_start=0, tlen=0, tags=None):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    # SEQ length must equal the number of query-consuming CIGAR bases
    # (ops M=0, I=1, S=4); D=2 / N=3 consume reference only.
    qlen = sum(length for op, length in cigar if op in (0, 1, 4))
    seq, qual = _seq_and_qual(qlen)
    a.query_sequence = seq
    a.query_qualities = qual
    a.flag = 0
    a.reference_id = 0
    a.reference_start = ref_start
    a.mapping_quality = mapq
    a.cigartuples = cigar
    a.is_reverse = is_reverse
    if flag_paired:
        a.is_paired = True
        a.is_proper_pair = True
        a.is_read1 = is_read1
        a.is_read2 = not is_read1
        a.mate_is_unmapped = mate_unmapped
        a.mate_is_reverse = mate_reverse
        a.next_reference_id = next_ref_id
        a.next_reference_start = next_ref_start
        a.template_length = tlen
    else:
        a.is_paired = False
        a.mate_is_unmapped = True
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
    if tags:
        a.set_tags(tags)
    return a


def make_bam(path, contig: str = CONTIG, length: int = CONTIG_LENGTH):
    """Write a small coordinate-sorted, indexed BAM at ``path``.

    Returns the BAM path as ``str``. A sibling ``<path>.bai`` index is created.
    The reads are described inline; see the module docstring for the rationale.
    """
    path = str(path)
    header = _header(contig, length)
    reads = []

    # --- Dense full-length reads so coverage over the window is clearly > 0
    #     and the random baseline sampling lands on a positive depth. Each
    #     spans [0, 4999), i.e. the whole contig.
    for i in range(8):
        reads.append(_new_read(header, f"dense{i}", 0, [(0, length - 1)]))

    # 1) Plain forward read fully covering the SV start breakpoint region.
    #    300M starting at 950 -> spans 950..1250.
    reads.append(_new_read(header, "plain_fwd_cover", 950, [(0, 300)]))

    # 2) Reverse read, soft-clipped on the left near the start breakpoint.
    #    10S + 250M starting at 1000 -> ref 1000..1250, left soft-clip at 990..999.
    reads.append(_new_read(header, "rev_softclip", 1000, [(4, 10), (0, 250)],
                           is_reverse=True, mapq=55))

    # 3) Read with an insertion near the start breakpoint (op 1).
    #    100M 5I 100M starting at 1020 -> ref 1020..1220.
    reads.append(_new_read(header, "with_insertion", 1020, [(0, 100), (1, 5), (0, 100)]))

    # 4) Read with a deletion near the start breakpoint (op 2).
    #    80M 10D 80M starting at 980 -> ref 980..1150.
    reads.append(_new_read(header, "with_deletion", 980, [(0, 80), (2, 10), (0, 80)]))

    # 5) Split read with a right soft-clip; SA tag points back into bin II
    #    ([sv_start+2, sv_start+52] = [1002, 1052]).
    #    Aligned 150M at 1010, then 20S clipped -> clip span 1160..1179.
    reads.append(_new_read(header, "split_read", 1010, [(0, 150), (4, 20)],
                           tags=[("SA", "chr1,1020,+,20M,60,0", "Z")]))

    # 6) Low mapping-quality read (below the mapq>=20 coverage threshold).
    reads.append(_new_read(header, "low_mapq", 1005, [(0, 180)], mapq=5))

    # --- Paired-end mates around the start breakpoint for orientation buckets.

    # 7) FF discordant pair (both forward). Mate lands in bin III/IV region.
    reads.append(_new_read(header, "pair_ff", 1010, [(0, 100)],
                           is_reverse=False, flag_paired=True, mate_unmapped=False,
                           mate_reverse=False, next_ref_start=1360, tlen=450))

    # 8) RR discordant pair (both reverse).
    reads.append(_new_read(header, "pair_rr", 1015, [(0, 100)],
                           is_reverse=True, flag_paired=True, mate_unmapped=False,
                           mate_reverse=True, next_ref_start=1380, tlen=465))

    # 9) RF discordant pair: read1, forward, reverse mate, positive tlen -> RF.
    reads.append(_new_read(header, "pair_rf", 1020, [(0, 100)],
                           is_reverse=False, flag_paired=True, mate_unmapped=False,
                           mate_reverse=True, is_read1=True, next_ref_start=1370,
                           tlen=350))

    # 10) Normal-ish FR proper pair (read1 forward, reverse mate, but read2 path
    #     not triggered; counts toward all_reads / insert sizes).
    reads.append(_new_read(header, "pair_normal", 1030, [(0, 100)],
                           is_reverse=True, flag_paired=True, mate_unmapped=False,
                           mate_reverse=False, is_read1=False, next_ref_start=1340,
                           tlen=410))

    # Coordinate-sort by reference_start before writing.
    reads.sort(key=lambda r: r.reference_start)

    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for r in reads:
            out.write(r)

    pysam.index(path)
    return path
