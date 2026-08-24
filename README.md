![dicast](docs/img/dicast-banner.png)

# A machine learning method for accurate detection of structural variants from short-read sequencing data

[![CI](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml/badge.svg)](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

**dicast** is a tool for scoring the confidence of structural variant (SV)
calls made by other SV callers (delly, manta, lumpy, cnvnator, gridss, or any
caller that emits a standard SV VCF). Given a coordinate-sorted, indexed
Illumina BAM plus one or more per-caller VCFs, it extracts alignment and
reference-context features around each call's breakpoints and scores it with
a pretrained XGBoost model, one model per SV type (DEL, DUP, INS). Every call gets a `dicast_qual` score in `[0, 1]`,
written as a column in a combined TSV and as a `DQ` INFO tag in re-emitted
VCFs. dicast is caller-agnostic (caller names are free-text labels), can
score several samples together with cross-sample rescue (e.g. trios), and
has a population mode that folds a PAV population catalog into the scoring.
The shipped annotations and models are hg38-specific.

## Installation

With conda (or mamba; swap the command accordingly):

```bash
git clone https://github.com/burgshrimps/dicast.git
cd dicast
conda env create -f environment.yml
conda activate dicast
```

Or with pip (Python 3.10+):

```bash
pip install git+https://github.com/burgshrimps/dicast.git
```

## Usage

*On first use dicast asks where to store the hg38 annotation files it
downloads once (~2.1 GB; press Enter to accept the suggested location, the
repo's `annot/` folder) and reuses them from then on.*

**Single sample:**

```bash
dicast call \
    --sample SAMPLE_NAME \
    --workdir WORKDIR \
    --fai /path/to/reference.fa.fai \
    --bam /path/to/sample.bam \
    --vcfs delly=/path/to/delly.vcf.gz manta=/path/to/manta.vcf.gz \
    --threads 24
```

`--vcfs` takes one or more `caller=path` entries: the input SV calls that
dicast scores. We recommend using the unfiltered calls of `manta`, `delly`,
`lumpy`, `gridss`, and `cnvnator`.

**Multiple samples:**

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
recall on common variants. DUP always falls back to the normal model:

```bash
dicast call ... --pop
```

`--pop-catalog` overrides the catalog VCF path (default:
`pav_catalog_hg38.vcf.gz` under `--annot-dir`). The catalog file itself will
be published alongside the annotation release; until then, `--pop` requires
you to supply your own via `--pop-catalog`.

## Output

A `call` run fills `--workdir` with a fixed tree of intermediate and final
files, all named `SAMPLE_REF.SVs.*` (`REF` defaults to `hg38`):

```
WORKDIR/
├── input/SAMPLE_REF.SVs.raw.tsv                      parsed input calls
├── features/
│   ├── ref/SAMPLE_REF.SVs.ref.tsv                    reference features
│   ├── aln/SAMPLE_REF.SVs.aln.ill.CHROM.SVTYPE.tsv   alignment feature shards
│   └── SAMPLE_REF.SVs.annot.tsv                      combined feature matrix
└── output/
    ├── SAMPLE_REF.SVs.dicast.tsv                     scores (the deliverable)
    ├── SAMPLE_CALLER.dicast.vcf                      per input caller, DQ-tagged
    └── SAMPLE_REF.SVs.dicast.merged.vcf              merged best-per-cluster VCF
```

`multi` builds the exact same tree per sample, under `WORKDIR/SAMPLE/...`.
The three files under `output/` are the ones you care about: the scores TSV,
the DQ-tagged per-caller VCFs, and the merged VCF.

### The scores TSV

**`output/SAMPLE_REF.SVs.dicast.tsv`** has one row per input call, for example:

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
| `sv_type` | `DEL`, `DUP`, or `INS` |
| `chrom`, `chrom_2`, `start`, `end` | breakpoint coordinates |
| `sv_len` | SV length (bp) |
| `filter` | the FILTER field from the source VCF record |
| `qual` | the caller-reported QUAL, if present |
| `dicast_qual` | dicast's confidence score, `[0, 1]` (`NA` if the model couldn't score the call) |
| `genotype` | genotype tuple, e.g. `(1, 1)` |

### The DQ-tagged VCFs

Each caller's input VCF is re-emitted into `output/` as
**`SAMPLE_CALLER.dicast.vcf`** (named after the sample and the `caller`
label from `--vcfs`, not the input filename, so two callers whose input
files happen to share a basename never collide), with a new INFO field:

```
##INFO=<ID=DQ,Number=1,Type=String,Description="Dicast Quality Score">
```

`DQ` carries the same `dicast_qual` value as the TSV; a record present in the
VCF but missing from the scored TSV (filtered out upstream) gets `DQ=-1`.

### Merged VCF

**`output/SAMPLE_REF.SVs.dicast.merged.vcf`** collapses the scores TSV down
to one call per real-world SV event: calls across all callers (and, in
`multi` mode, rescued calls) are clustered per SV type (DEL/DUP by
>50% reciprocal breakpoint overlap, INS by <200bp breakpoint distance), and
only the highest-`dicast_qual` call in each cluster survives, genotype
included. Selection is population-aware: a `pav` population-catalog call
only wins its cluster if no non-population caller call in that cluster
reaches `dicast_qual >= 0.4`. It's a fresh, minimal VCF (not a merge of the
input VCFs' headers) with a `CALLER` INFO tag naming which caller produced
the winning call and a `DQ` INFO tag carrying its `dicast_qual`.

### Rescued calls

In `multi` mode, calls found in one sample but missing from another sample's
own callers are transplanted into that sample's candidate set and scored
against its own BAM, recovering calls the sample's callers missed but that
show up as real signal once you actually look at that sample's reads.
Rescued rows are distinguishable by their `caller` value:
`rescue:<origin sample>:<origin caller>` (e.g. a deletion found only by the
mother's delly call, rescued into the child, shows up as
`rescue:MOTHER:delly` in the child's output).
