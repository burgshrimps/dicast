![dicast](docs/img/dicast-banner.png)

**dicast**: a machine learning method for accurate detection of structural
variants from short-read sequencing data

---

[![CI](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml/badge.svg)](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

**dicast** is a tool for scoring the confidence of structural variant (SV)
calls made by other SV callers (delly, manta, lumpy, cnvnator, gridss, or any
caller that emits a standard SV VCF). Given a coordinate-sorted, indexed
Illumina BAM plus one or more per-caller VCFs, it extracts alignment and
reference-context features around each call's breakpoints and scores it with
a pretrained XGBoost model, one model per SV type (DEL, DUP, INS; INV is
scored with a fixed rule). Every call gets a `dicast_qual` score in `[0, 1]`,
written as a column in a combined TSV and as a `DQ` INFO tag in re-emitted
VCFs. dicast is caller-agnostic (caller names are free-text labels), can
score several samples together with cross-sample rescue (e.g. trios), and
has a population mode that folds a PAV population catalog into the scoring.
The shipped annotations and models are hg38-specific.

## Installation

With conda (or mamba — swap the command accordingly):

```bash
git clone https://github.com/burgshrimps/dicast.git
cd dicast
conda env create -f environment.yml
conda activate dicast
```

## Usage

*On first use dicast asks where to store the hg38 annotation files it
downloads once (~2.1 GB; press Enter to accept the suggested location, the
repo's `annot/` folder) and reuses them from then on. `--models` defaults to
the `models/` directory shipped in this repo.*

**Single sample:**

```bash
dicast call \
    --sample SAMPLE_NAME \
    --workdir WORKDIR \
    --fai /path/to/reference.fa.fai \
    --bam /path/to/sample.bam \
    --vcfs delly=/path/to/delly.vcf.gz manta=/path/to/manta.vcf.gz \
    --threads 8
```

`--vcfs` takes one or more `caller=path` entries; the label you pick is what
shows up in the output `caller` column and decides which re-emitted VCF gets
which scores.

**Multiple samples (e.g. a trio) with cross-sample rescue** — every sample is
additionally scored on calls that only the *other* samples' callers found
(see [Rescued calls](#rescued-calls)); at least two samples are required:

```bash
dicast multi \
    --bams MOTHER=/path/mother.bam FATHER=/path/father.bam CHILD=/path/child.bam \
    --vcfs MOTHER:delly=/path/mother_delly.vcf.gz MOTHER:manta=/path/mother_manta.vcf.gz \
           FATHER:delly=/path/father_delly.vcf.gz FATHER:manta=/path/father_manta.vcf.gz \
           CHILD:delly=/path/child_delly.vcf.gz  CHILD:manta=/path/child_manta.vcf.gz \
    --workdir WORKDIR --fai /path/to/reference.fa.fai
```

`--bams` takes `sample=bam_file` entries and `--vcfs` takes
`sample:caller=vcf_file` entries.

**Population mode** (`call` and `multi`): `--pop` adds the PAV
structural-variant population catalog as an extra "caller" (labeled `pav`)
and switches DEL/INS scoring to the population-aware models, which favour
recall on common variants. DUP and INV always fall back to the normal
model / rule:

```bash
dicast call ... --pop
```

`--pop-catalog` overrides the catalog VCF path (default:
`pav_catalog_hg38.vcf.gz` under `--annot-dir`). The catalog file itself will
be published alongside the annotation release; until then, `--pop` requires
you to supply your own via `--pop-catalog`.

Other useful flags: `--chrom` restricts work to specific chromosomes
(default: all standard chromosomes through chrX), `--sv_types` to specific
SV types (default: `DEL DUP INS INV`), and `--threads` parallelises the
per-chromosome, per-SV-type feature collection. For detailed descriptions of
all parameters and their defaults, run `dicast call --help` and
`dicast multi --help`.

## Interpreting dicast output

A run fills `--workdir` with intermediate and final files named
`SAMPLE_REF.SVs.*` (`REF` defaults to `hg38`). The two you care about are the
scores TSV and the DQ-tagged VCFs.

### The scores TSV

**`SAMPLE_REF.SVs.dicast.tsv`** has one row per input call, for example:

```
id            cohort  sample  reference  technology  caller  sv_type  chrom  chrom_2  start     end       sv_len  filter  qual  dicast_qual  genotype
DEL00000001   none    demo    hg38       ill         delly   DEL      chr21  NA       9101395   9101696   301     PASS    NA    0.102        (1, 1)
DEL00000002   none    demo    hg38       ill         delly   DEL      chr21  NA       9616899   9618664   1765    PASS    NA    0.005        (1, 1)
DEL00000003   none    demo    hg38       ill         delly   DEL      chr21  NA       13278442  13278523  81      PASS    NA    0.637        (1, 1)
```

| column | meaning |
|---|---|
| `id` | dicast-internal call id (e.g. `DEL00000001`) |
| `cohort` | `--cohort` value (default `none`) |
| `sample` | `--sample` value |
| `reference` | `--ref` value |
| `technology` | `--technology` value (default `ill`) |
| `caller` | the `caller` label from `--vcfs`, or `rescue:<origin sample>:<origin caller>` for a rescued row |
| `sv_type` | `DEL`, `DUP`, `INS`, or `INV` |
| `chrom`, `chrom_2`, `start`, `end` | breakpoint coordinates |
| `sv_len` | SV length (bp) |
| `filter` | the FILTER field from the source VCF record |
| `qual` | the caller-reported QUAL, if present |
| `dicast_qual` | dicast's confidence score, `[0, 1]` (`NA` if the model couldn't score the call) |
| `genotype` | genotype tuple, e.g. `(1, 1)` |

### The DQ-tagged VCFs

Each caller's input VCF is re-emitted into `--workdir` as
`<name>.dicast.vcf` (`.vcf`/`.vcf.gz` suffix replaced), with a new INFO
field:

```
##INFO=<ID=DQ,Number=1,Type=String,Description="Dicast Quality Score">
```

`DQ` carries the same `dicast_qual` value as the TSV; a record present in the
VCF but missing from the scored TSV (filtered out upstream) gets `DQ=-1`.

### Rescued calls

In `multi` mode, calls found in one sample but missing from another sample's
own callers are transplanted into that sample's candidate set and scored
against its own BAM — recovering calls the sample's callers missed but that
show up as real signal once you actually look at that sample's reads.
Rescued rows are distinguishable by their `caller` value:
`rescue:<origin sample>:<origin caller>` (e.g. a deletion found only by the
mother's delly call, rescued into the child, shows up as
`rescue:MOTHER:delly` in the child's output).

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
(see `Dicast.score_inversions` in `dicast/model.py`).

Swap in your own models by pointing `--models` at a directory with the same
`dicast_<SV_TYPE>[_pop].json` naming convention; `dicast/model.py` is
the training entry point.

## Reference annotations

dicast's feature engineering depends on eight hg38 annotation files:

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

None of these are tracked in git — they are release assets, and dicast
downloads the ones it needs automatically on first use (~2.1 GB total,
checksum-verified, one time), asking once where to store them. Set the
`DICAST_DATA_DIR` environment variable to fix the location, run
`dicast-fetch-annotations` to prefetch everything (e.g. on a cluster ahead of
offline jobs), or point `--annot-dir` (or the individual flags) at your own
copies.
