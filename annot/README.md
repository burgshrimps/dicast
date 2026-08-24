# hg38 annotation directory

Reference annotation files consumed by dicast's `ReferenceAnnotator`
(`dicast/collect_reference.py`) for feature engineering. All are public;
sources are listed per file below.

## Shipped files

| file                       | source                                                                 |
| -------------------------- | ----------------------------------------------------------------------- |
| `hg38_repeatmasker.tsv`    | UCSC RepeatMasker track. Filtered by `repClass` into transposable-element / low-complexity features (`LINE`, `SINE`, `LTR`, `Retroposon`, `Low_complexity`) and, within `repClass == "Satellite"`, by `repFamily` into satellite subfamilies (`centr` centromeric α-satellite, `telo` telomeric/subtelomeric, `acro` acrocentric rDNA-adjacent). Ambiguous classes (trailing-`?` labels, `Unknown`) are dropped before use. |
| `hg38_strs_chaisson.bed`   | Chaisson et al. STR catalog (1–6 bp motifs, lumped).                     |
| `hg38_vntrs_chaisson.bed`  | Chaisson et al. VNTR catalog.                                            |
| `hg38_cpg_islands.tsv`     | UCSC CpG-island track, used as-is.                                       |
| `hg38_centromeres.tsv`     | UCSC Genome Browser Track; distance-to-centromere feature.              |
| `hg38_asmb_gaps.tsv`       | UCSC Genome Browser Track; distance-to-assembly-gap feature.            |
| `hg38_alt_haps.tsv`        | UCSC Genome Browser Track; distance-to-alternative-haplotype feature.   |
| `hg38_gc_content.bw`       | GC-content BigWig, queried for local GC-content features around each call's breakpoints. |

## Large files: distributed via GitHub Release

`hg38_gc_content.bw` (~1.6 GB) and `hg38_repeatmasker.tsv` (~460 MB) exceed
GitHub's 100 MB file-size limit and are not tracked in this repo (see
`.gitignore`). They're attached to the `annotations-v1` GitHub Release and
fetched on demand by `download_annotations.sh` at the repo root, which
verifies each download against `checksums.md5` in this directory.

`pav_catalog_hg38.vcf.gz` — the PAV population catalog used by the `--pop`
flag — will be published as a release asset the same way once available;
until then `download_annotations.sh` skips it with an informative note.
