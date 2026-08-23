# hg38 annotation directory

Reference annotation files for SV genomic-context analysis. Three source types coexist here:

1. **RepeatMasker / UCSC / Chaisson** — long-standing files used by the existing tGenVar pipeline.
2. **GIAB v3.6 stratifications** — downloaded from https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/ for the genomic-context analysis.

The unified catalog at `04_figures/genomic_context_catalog.tsv` (regenerate with `04_figures/build_annotation_catalog.py`) draws intervals from the files below and tags each row with a `name` (stratum label).

## GIAB v3.6 names and their source BEDs

`name` is the column used in `genomic_context_catalog.tsv`; `source BED` is the file on disk in this directory.

### Homopolymers & STRs (motif-split)

| catalog `name`     | source BED(s)                                                                                          | what it means                                                                                                    |
| ------------------ | ------------------------------------------------------------------------------------------------------ | ---------------------------------------------------------------------------------------------------------------- |
| `giab_homopolymer` | `GRCh38_AllHomopolymers_ge7bp_imperfectge11bp_slop5.bed.gz`                                            | Homopolymer runs (1-bp motif): perfect runs ≥ 7 bp **or** imperfect runs ≥ 11 bp.                                |
| `giab_STR_di`      | `GRCh38_SimpleRepeat_diTR_10to49_slop5.bed.gz` ∪ `..._50to149_slop5.bed.gz` ∪ `..._ge150_slop5.bed.gz` | Dinucleotide STRs (2-bp motif, e.g. CACA). Union of three array-length bins, merged via `bioframe.merge`.         |
| `giab_STR_tri`     | `GRCh38_SimpleRepeat_triTR_14to49_slop5.bed.gz` ∪ `..._50to149_slop5.bed.gz` ∪ `..._ge150_slop5.bed.gz` | Trinucleotide STRs (3-bp motif, e.g. CAG). Union of three array-length bins.                                     |
| `giab_STR_quad`    | `GRCh38_SimpleRepeat_quadTR_19to49_slop5.bed.gz` ∪ `..._50to149_slop5.bed.gz` ∪ `..._ge150_slop5.bed.gz` | Tetranucleotide STRs (4-bp motif). Union of three array-length bins.                                             |

### Tandem repeats by array length (disjoint, all motifs)

| catalog `name`      | source BED                                                    | what it means                                                   |
| ------------------- | ------------------------------------------------------------- | --------------------------------------------------------------- |
| `giab_TR_le50`      | `GRCh38_AllTandemRepeats_le50bp_slop5.bed.gz`                 | Any tandem-repeat array ≤ 50 bp long.                           |
| `giab_TR_51to200`   | `GRCh38_AllTandemRepeats_51to200bp_slop5.bed.gz`              | Array 51 – 200 bp long.                                         |
| `giab_TR_201to10kb` | `GRCh38_AllTandemRepeats_201to10000bp_slop5.bed.gz`           | Array 201 bp – 10 kb long.                                      |
| `giab_TR_ge10kb`    | `GRCh38_AllTandemRepeats_ge10001bp_slop5.bed.gz`              | Array ≥ 10 kb long (big-VNTR territory).                        |
| `giab_TR_any`       | `GRCh38_AllTandemRepeats.bed.gz`                              | Umbrella: any tandem repeat, any motif, any length.             |

### Satellites & segdups

| catalog `name`       | source BED                                          | what it means                                                                      |
| -------------------- | --------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `giab_satellites`    | `GRCh38_satellites_slop5.bed.gz`                    | GIAB-curated human satellite catalog (α-satellite, HSat1/2/3, etc.).               |
| `giab_segdups`       | `GRCh38_segdups.bed.gz`                             | Segmental duplications (UCSC `genomicSuperDups`, GIAB-versioned).                  |

### Baseline

| catalog `name` | source BED                                    | what it means                                                                          |
| -------------- | --------------------------------------------- | -------------------------------------------------------------------------------------- |
| `giab_easy`    | `GRCh38_notinalldifficultregions.bed.gz`      | Complement of GIAB's union of all "difficult" strata — the easy part of the genome.   |

All GIAB BEDs carry the `_slop5` suffix where applicable, meaning GIAB pre-extended the intervals by 5 bp on each side before publishing; using strict overlap against these BEDs effectively gives a ±5 bp margin for free.

## RepeatMasker / UCSC / Chaisson sources

Existing files; used by the catalog for the non-GIAB strata.

### Transposable elements & low complexity (`repClass` split)

| catalog `name`        | source file            | filter (`repClass`)     |
| --------------------- | ---------------------- | ----------------------- |
| `rep_LINE`            | `hg38_repeatmasker.tsv`| `LINE`                  |
| `rep_SINE`            | `hg38_repeatmasker.tsv`| `SINE`                  |
| `rep_LTR`             | `hg38_repeatmasker.tsv`| `LTR`                   |
| `rep_Retroposon`      | `hg38_repeatmasker.tsv`| `Retroposon`            |
| `rep_Low_complexity`  | `hg38_repeatmasker.tsv`| `Low_complexity`        |

### Satellite subfamilies (`repClass == "Satellite"`, `repFamily` split)

GIAB's `giab_satellites` is a single umbrella; RepeatMasker splits satellite DNA into biologically meaningful subfamilies, added here for finer-grained stratification.

| catalog `name`                 | source file            | filter                                   | what it captures                                       |
| ------------------------------ | ---------------------- | ---------------------------------------- | ------------------------------------------------------ |
| `rep_satellite_centromeric`    | `hg38_repeatmasker.tsv`| `repClass=="Satellite" & repFamily=="centr"` | Centromeric α-satellite (ALR) arrays.                   |
| `rep_satellite_telomeric`      | `hg38_repeatmasker.tsv`| `repClass=="Satellite" & repFamily=="telo"`  | Telomeric / subtelomeric repeat arrays.                 |
| `rep_satellite_acrocentric`    | `hg38_repeatmasker.tsv`| `repClass=="Satellite" & repFamily=="acro"`  | Acrocentric (rDNA-adjacent) satellites on chrs 13/14/15/21/22. Very sparse (~60 intervals). |

### UCSC / Chaisson

| catalog `name` | source file                | notes                                                                                                                                |
| -------------- | -------------------------- | ------------------------------------------------------------------------------------------------------------------------------------ |
| `cpg_islands`  | `hg38_cpg_islands.tsv`     | UCSC CpG-island track, used as-is.                                                                                                   |
| `rep_STR`      | `hg38_strs_chaisson.bed`   | Chaisson et al. STR catalog (1–6 bp motifs, lumped). Parallel to GIAB motif-resolved `giab_STR_{di,tri,quad}`.                       |
| `rep_VNTR`     | `hg38_vntrs_chaisson.bed`  | Chaisson et al. VNTR catalog. Parallel to GIAB tandem-repeat length bins (different curation).                                       |

## Other files (not consumed by the catalog)

Present in the dir but not referenced by `genomic_context_catalog.tsv`:

- `hg38_centromeres.tsv`, `hg38_asmb_gaps.tsv`, `hg38_alt_haps.tsv` — used by `dicast`'s `ReferenceAnnotator` for feature engineering, not stratification.
- `hg38_genomicSuperDups.tsv` — legacy UCSC segdup track; superseded by `giab_segdups` (same underlying track, GIAB-versioned).
- `hg38_genes.tsv`, `hg38_orphanet.tsv` — gene-level annotation, unrelated to context stratification.
- `hg38_gc_content.bw` — GC-content BigWig for dicast features.
- `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai` — FASTA index.

## GIAB strata we did NOT download

The GIAB v3.6 directory has many more files we skipped to keep the catalog readable:

- **LowComplexity/**: per-motif-length homopolymer splits (`homopolymer_4to6`, `_7to11`, `_ge12`, `_ge21` + AT/GC variants); imperfect-homopolymer splits; `AllTandemRepeatsandHomopolymers_slop5` umbrella; complement `notinAll*` files. The `AllTandemRepeats_ge101bp` bin is obsoleted by the new disjoint length bins above.
- **Mappability/**, **OtherDifficult/** (MHC/KIR/VDJ/gaps/assembly-collapse FP regions), **XY/**, **GCcontent/**, **Functional/**, **FunctionalTechnicallyDifficult/** — potentially useful for future expansions; see Tier-1 shortlist in the session notes.
- **Ancestry/** — population-specific, not relevant to our 9-patient cohort.
- **GenomeSpecific/** — HG002/HG003/HG004-specific; our cohort is not HG002.
