# dicast

Confidence scoring for structural-variant calls from short-read (Illumina) sequencing.

[![CI](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml/badge.svg)](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml)
[![License: GPL-3.0](https://img.shields.io/badge/License-GPL--3.0-blue.svg)](LICENSE)
[![Python 3.10](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

## What is dicast

dicast is a command-line tool that scores the confidence of structural-variant
(SV) calls made by other SV callers (delly, manta, lumpy, cnvnator, gridss,
or any other caller that emits a standard SV VCF). Given a coordinate-sorted,
indexed Illumina BAM plus one or more per-caller VCFs, it extracts alignment
and reference-context features around each call's breakpoints and scores it
with a pretrained XGBoost model, one model per SV type (DEL, DUP, INS;
INV is scored with a fixed rule rather than a trained model). Every call
gets a `dicast_qual` score in `[0, 1]`, written as a column in a combined TSV
and, if you re-emit the input VCFs, as a `DQ` INFO tag. Scoring is
caller-agnostic and cohort-agnostic — the shipped models are trained
artifacts, but any directory of matching XGBoost JSON files works with
`--models`.

## Requirements

- Linux or macOS
- [conda](https://docs.conda.io/) or [mamba](https://mamba.readthedocs.io/)
- ~2.3 GB free disk space if you download the full reference annotations
  (the quickstart demo below needs none of that)
- Your own data: a coordinate-sorted, indexed hg38 BAM, at least one SV
  caller's VCF for that sample, and a `.fai` index for the hg38 FASTA you
  aligned against

dicast currently supports hg38 only — the shipped reference annotations and
trained models are all hg38-specific.

## Installation

```bash
git clone https://github.com/burgshrimps/dicast.git
cd dicast

conda env create -f environment.yml
conda activate dicast
```

(If you use mamba instead of conda, swap in `mamba env create -f
environment.yml` / `mamba activate dicast`.)

```bash
bash download_annotations.sh
```

This fetches the two hg38 annotation files too large for git (~2.1 GB
combined: GC content and RepeatMasker). **The quickstart below works without
this step** — it uses a small, pre-sliced set of annotations under
`test_data/annot/`.

## Quickstart

A small, self-contained chr21-only demo (real HG002 reads and calls) lets you
run dicast end-to-end with zero downloads:

```bash
python3 dicast.py call \
    --sample demo --workdir demo_out \
    --fai test_data/hg38.fa.fai \
    --bam test_data/demo.bam \
    --vcfs delly=test_data/demo_delly.vcf.gz \
    --annot-dir test_data/annot \
    --chrom chr21 --sv_types DEL INS --threads 2
```

This writes several files to `demo_out/`; the one you want is
`demo_hg38.SVs.dicast.tsv`, one row per call with a `dicast_qual` score:

```
id            cohort  sample  reference  technology  caller  sv_type  chrom  chrom_2  start     end       sv_len  filter  qual  dicast_qual  genotype
DEL00000001   none    demo    hg38       ill         delly   DEL      chr21  NA       9101395   9101696   301     PASS    NA    0.102        (1, 1)
DEL00000002   none    demo    hg38       ill         delly   DEL      chr21  NA       9616899   9618664   1765    PASS    NA    0.005        (1, 1)
DEL00000003   none    demo    hg38       ill         delly   DEL      chr21  NA       13278442  13278523  81      PASS    NA    0.637        (1, 1)
```

It also re-emits `test_data/demo_delly.vcf.gz` as
`demo_out/demo_delly.dicast.vcf`, the same VCF with a `DQ` INFO field added
per record (e.g. `DQ=0.637`). See `test_data/README.md` for exactly how this
demo dataset was built and why it uses artificially downsampled coverage at
the variant loci.

## Usage on your own data

```bash
python3 dicast.py call \
    --sample SAMPLE_NAME \
    --workdir WORKDIR \
    --fai /path/to/reference.fa.fai \
    --bam /path/to/sample.bam \
    --vcfs delly=/path/to/delly.vcf.gz manta=/path/to/manta.vcf.gz \
    --threads 8
```

See `dicast.sh` for a fuller worked example (five callers).

Notes:

- `--vcfs` takes one or more `caller=path` entries. `caller` is a free-text
  label — it does not have to be a caller dicast was trained on, but the
  label you pick here is what shows up in the output `caller` column and,
  for `add_info_tag_to_vcf`, which VCF gets tagged with which subset of
  scores. `--vcfs` is required; the caller name and file path are joined
  with `=`, e.g. `delly=formatted_variants.vcf.gz`.
- `--annot-dir` defaults to the `annot/` directory shipped in this repo,
  and `--models` defaults to the `models/` directory shipped in this repo
  — you only need to pass either flag to point at a different set. You can
  also override individual annotation files (`--repeats`, `--cgis`,
  `--centromeres`, `--gaps`, `--althaps`, `--vntrs`, `--strs`, `--gc`)
  without replacing the whole directory.
  Individual flags always win over `--annot-dir`; both `--annot-dir` and
  `--models` require hg38 if you're using the shipped resources.
- `--chrom` restricts feature extraction and prediction to specific
  chromosomes (default: all standard chromosomes through chrX). `--sv_types`
  restricts which SV types are processed (default: `DEL DUP INS INV`).
- `--threads` controls parallelism across the per-chromosome,
  per-SV-type alignment feature collection jobs.

## Outputs

For a `call` run with `--sample SAMPLE --ref REF` (`REF` defaults to
`hg38`), `--workdir` fills up with intermediate and final files named
`SAMPLE_REF.SVs.*`. The one you care about is
**`SAMPLE_REF.SVs.dicast.tsv`**, one row per input call:

| column | meaning |
|---|---|
| `id` | dicast-internal call id (e.g. `DEL00000001`) |
| `cohort` | `--cohort` value (default `none`) |
| `sample` | `--sample` value |
| `reference` | `--ref` value |
| `technology` | `--technology` value (default `ill`) |
| `caller` | the `caller` label from `--vcfs`, or `rescue:<origin sample>:<origin caller>` for a cross-sample rescue row (multi-sample mode only) |
| `sv_type` | `DEL`, `DUP`, `INS`, or `INV` |
| `chrom`, `chrom_2`, `start`, `end` | breakpoint coordinates |
| `sv_len` | SV length (bp) |
| `filter` | the FILTER field from the source VCF record |
| `qual` | the caller-reported QUAL, if present |
| `dicast_qual` | dicast's confidence score, `[0, 1]` (`NA` if the model couldn't score the call) |
| `genotype` | genotype tuple, e.g. `(1, 1)` |

If you re-emit the input VCFs (which `dicast.py call` always does), each
caller's VCF gets rewritten into the `--workdir` as `<name>.dicast.vcf`
(`.vcf`/`.vcf.gz` suffix replaced), with a new INFO field:

```
##INFO=<ID=DQ,Number=1,Type=String,Description="Dicast Quality Score">
```

`DQ` carries the same `dicast_qual` value as the TSV; a record present in the
VCF but missing from the scored TSV (filtered out upstream) gets `DQ=-1`.

## Population mode

`--pop` adds the PAV structural-variant population catalog as an extra
"caller" (labeled `pav`) and switches DEL/INS scoring to the
population-aware models (`dicast_DEL_pop.json` / `dicast_INS_pop.json`),
which were additionally trained with common population variants as
positives, favouring recall on common variants. DUP and INV have no
population-aware model and always fall back to the normal model / rule.

```bash
python3 dicast.py call ... --pop
```

`--pop-catalog` overrides the catalog VCF path (default:
`pav_catalog_hg38.vcf.gz` under `--annot-dir`). The catalog file itself will
be published alongside the annotation release; until then, `--pop` requires
you to supply your own via `--pop-catalog`.

## Multi-sample mode

`multi` scores several samples together (e.g. a trio) with **cross-sample
rescue**: for every cluster of matching calls found across the samples'
own callers, any sample missing from that cluster gets the call added to
its own list and scored against its own BAM too. This recovers calls a
sample's callers missed but that show up as real signal once you actually
look at that sample's reads.

```bash
python3 dicast.py multi \
    --bams MOTHER=/path/mother.bam FATHER=/path/father.bam CHILD=/path/child.bam \
    --vcfs MOTHER:delly=/path/mother_delly.vcf.gz MOTHER:manta=/path/mother_manta.vcf.gz \
           FATHER:delly=/path/father_delly.vcf.gz FATHER:manta=/path/father_manta.vcf.gz \
           CHILD:delly=/path/child_delly.vcf.gz  CHILD:manta=/path/child_manta.vcf.gz \
    --workdir WORKDIR --fai /path/to/reference.fa.fai
```

`--bams` takes `sample=bam_file` entries (at least two samples required);
`--vcfs` takes `sample:caller=vcf_file` entries. See `dicast_multi.sh` for a
full worked example. In the output TSV, rescued rows are distinguishable by
their `caller` value: `rescue:<origin sample>:<origin caller>` (e.g. a
deletion found only by the mother's `delly` call, rescued into the child,
shows up as `rescue:MOTHER:delly` in the child's output).

## Models

Shipped in `models/`, one XGBoost JSON per SV type (plus population-aware
variants for DEL/INS):

| model | SV type | trained on |
|---|---|---|
| `dicast_DEL.json` | DEL | tGenVar cohort (8 samples), delly/manta/lumpy/cnvnator/gridss calls |
| `dicast_DUP.json` | DUP | tGenVar cohort (7 samples), delly/manta/lumpy/cnvnator/gridss calls |
| `dicast_INS.json` | INS | tGenVar cohort (8 samples), delly/manta/gridss calls |
| `dicast_DEL_pop.json` | DEL (population-aware) | tGenVar cohort + PAV population calls as additional positives |
| `dicast_INS_pop.json` | INS (population-aware) | tGenVar cohort + PAV population calls as additional positives |

INV has no trained model — it's scored with a fixed rule based on clipped
reads, discordant-pair signal, coverage, mapping quality, and split reads
(see `Dicast.score_inversions` in `dicast_lib/model.py`).

Each model has a sibling `*_metadata.json` (training date, feature list,
cohorts, callers, and positive/negative call counts). Swap in your own
models by pointing `--models` at a directory with the same
`dicast_<SV_TYPE>[_pop].json` naming convention; `dicast_lib/model.py` is
the training entry point.

## Reference annotations

`annot/` ships the eight hg38 annotation files dicast's feature engineering
depends on:

| file | source |
|---|---|
| `hg38_repeatmasker.tsv` | UCSC RepeatMasker track (transposable-element and low-complexity classes, plus satellite subfamilies) |
| `hg38_strs_chaisson.bed` | Chaisson et al. short tandem repeat catalog |
| `hg38_vntrs_chaisson.bed` | Chaisson et al. VNTR catalog |
| `hg38_cpg_islands.tsv` | UCSC CpG island track |
| `hg38_centromeres.tsv` | UCSC Genome Browser centromere track |
| `hg38_asmb_gaps.tsv` | UCSC Genome Browser assembly-gap track |
| `hg38_alt_haps.tsv` | UCSC Genome Browser alternative-haplotype track |
| `hg38_gc_content.bw` | GC-content BigWig |

Six of these are tracked in git under `annot/`. `hg38_gc_content.bw`
(~1.6 GB) and `hg38_repeatmasker.tsv` (~460 MB) exceed GitHub's 100 MB file
limit and are fetched on demand from a GitHub Release by
`bash download_annotations.sh`, which verifies each against
`annot/checksums.md5`. Full source details are in `annot/README.md`.
`--annot-dir` defaults to this repo's `annot/`; individual `--repeats` /
`--cgis` / `--centromeres` / `--gaps` / `--althaps` / `--vntrs` / `--strs` /
`--gc` flags override one file at a time.

## Full CLI reference

<details>
<summary><code>python3 dicast.py call --help</code></summary>

```
usage: dicast.py call [-h] [--cohort COHORT] --sample SAMPLE
                      [--chrom CHROM [CHROM ...]] [--ref REF]
                      [--technology TECHNOLOGY] --workdir WORKDIR --fai FAI
                      [--annot-dir ANNOT_DIR] [--repeats REPEATS]
                      [--cgis CGIS] [--centromeres CENTROMERES] [--gaps GAPS]
                      [--althaps ALTHAPS] [--vntrs VNTRS] [--strs STRS]
                      [--gc GC] --bam BAM --vcfs [VCFS ...] [--models MODELS]
                      [--threads THREADS] [--pop] [--pop-catalog POP_CATALOG]
                      [--benchmark BENCHMARK]
                      [--sv_types SV_TYPES [SV_TYPES ...]]

options:
  -h, --help            show this help message and exit
  --cohort COHORT       Cohort name
  --sample SAMPLE       Sample name
  --chrom CHROM [CHROM ...]
                        Chromosomes
  --ref REF             Reference genome name
  --technology TECHNOLOGY
                        Sequencing technology name
  --workdir WORKDIR     Working and output directory
  --fai FAI             FAI file of the reference genome
  --annot-dir ANNOT_DIR
                        Directory with default hg38 annotation files
  --repeats REPEATS     TSV file with repeats annotated by repeatmasker
  --cgis CGIS           TSV file with CpG island annotations
  --centromeres CENTROMERES
                        TSV file with centromere annotations
  --gaps GAPS           TSV file with assembly gap annotations
  --althaps ALTHAPS     TSV file with alternative haplotype annotations
  --vntrs VNTRS         BED file with VNTR regions from Chaisson
  --strs STRS           BED file with STR regions from Chaisson
  --gc GC               BIGWIG file with GC content
  --bam BAM             BAM file
  --vcfs [VCFS ...]     List of VCF files. Needs to be in the format
                        method=vcf_file
  --models MODELS       Directory with trained models
  --threads THREADS     Number of threads to use
  --pop                 Add the PAV population catalog as an additional caller
                        and prefer population-aware models
  --pop-catalog POP_CATALOG
                        VCF file with the PAV population catalog
  --benchmark BENCHMARK
                        Path to write a per-stage TSV with wall-time, CPU-time
                        and peak RSS (feature_collection, prediction, total).
                        If unset, no benchmark is written.
  --sv_types SV_TYPES [SV_TYPES ...]
                        Restrict feature extraction and prediction to these SV
                        types. Default: DEL DUP INS INV.
```

</details>

<details>
<summary><code>python3 dicast.py multi --help</code></summary>

```
usage: dicast.py multi [-h] [--cohort COHORT] --bams [BAMS ...] --vcfs
                       [VCFS ...] [--chrom CHROM [CHROM ...]] [--ref REF]
                       [--technology TECHNOLOGY] --workdir WORKDIR --fai FAI
                       [--annot-dir ANNOT_DIR] [--repeats REPEATS]
                       [--cgis CGIS] [--centromeres CENTROMERES] [--gaps GAPS]
                       [--althaps ALTHAPS] [--vntrs VNTRS] [--strs STRS]
                       [--gc GC] [--models MODELS] [--threads THREADS] [--pop]
                       [--pop-catalog POP_CATALOG] [--benchmark BENCHMARK]
                       [--sv_types SV_TYPES [SV_TYPES ...]]

options:
  -h, --help            show this help message and exit
  --cohort COHORT       Cohort name
  --bams [BAMS ...]     Per-sample BAM files. Format: sample=bam_file
  --vcfs [VCFS ...]     Per-sample, per-caller VCF files. Format:
                        sample:caller=vcf_file
  --chrom CHROM [CHROM ...]
                        Chromosomes
  --ref REF             Reference genome name
  --technology TECHNOLOGY
                        Sequencing technology name
  --workdir WORKDIR     Working and output directory
  --fai FAI             FAI file of the reference genome
  --annot-dir ANNOT_DIR
                        Directory with default hg38 annotation files
  --repeats REPEATS     TSV file with repeats annotated by repeatmasker
  --cgis CGIS           TSV file with CpG island annotations
  --centromeres CENTROMERES
                        TSV file with centromere annotations
  --gaps GAPS           TSV file with assembly gap annotations
  --althaps ALTHAPS     TSV file with alternative haplotype annotations
  --vntrs VNTRS         BED file with VNTR regions from Chaisson
  --strs STRS           BED file with STR regions from Chaisson
  --gc GC               BIGWIG file with GC content
  --models MODELS       Directory with trained models
  --threads THREADS     Number of threads to use
  --pop                 Add the PAV population catalog as an additional caller
                        (for every sample) and prefer population-aware models
  --pop-catalog POP_CATALOG
                        VCF file with the PAV population catalog
  --benchmark BENCHMARK
                        Path to write a per-stage TSV with wall-time, CPU-time
                        and peak RSS (feature_collection, prediction, total).
                        If unset, no benchmark is written.
  --sv_types SV_TYPES [SV_TYPES ...]
                        Restrict feature extraction and prediction to these SV
                        types. Default: DEL DUP INS INV.
```

</details>

## Reproducing the paper benchmark

`dicast_benchmark_hg002.sh` reproduces the runtime / peak-memory benchmark
reported in the paper on HG002 (30x WGS Illumina) — see that script for the
input paths to configure.

## Citation

The dicast paper has been accepted in principle at *Genome Biology*
(DOI to be assigned on publication). Until then, please cite this repository.

```bibtex
@article{dicast,
  title   = {dicast: TBD},
  author  = {TBD},
  journal = {Genome Biology},
  year    = {2026},
  note    = {Accepted in principle. DOI TBD.}
}
```

## License

GPL-3.0. See [LICENSE](LICENSE).
