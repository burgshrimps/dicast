![dicast](docs/img/dicast-banner.png)

# A machine learning method for accurate detection of structural variants from short-read sequencing data

[![CI](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml/badge.svg)](https://github.com/burgshrimps/dicast/actions/workflows/ci.yml)
[![License: GPL-3.0](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-3.10-blue.svg)](environment.yml)

**dicast** is a machine learning method for determining whether a structural
variant (SV) call is a real variant or a likely false-positive artefact. For
every input variant (we recommend combining several SV callers, see below),
dicast builds an internal representation from over 100 features describing
the genomic context and the alignment signal around the variant's
breakpoints, and scores it with a pretrained XGBoost model, one per SV type
(DEL, DUP, INS). Every variant receives an easily interpretable confidence score between 0
and 1. We show
that this approach outperforms individual SV callers as well as commonly
used consensus approaches.

*See also: to visually inspect the read-level evidence behind individual SV
calls, check out dicast's companion tool
[cuban](https://github.com/burgshrimps/cuban).*

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
    --workdir WORKDIR --fai /path/to/reference.fa.fai \
    --threads 24
```

`--bams` takes `sample=bam_file` entries and `--vcfs` takes
`sample:caller=vcf_file` entries.

In multi-sample mode dicast not only scores the variants present in each
sample's input VCFs but also variants that occur in the other samples of the
run. For example, if based on the input VCFs a variant supposedly occurs only
in the child and not the parents, dicast checks the same region in the
parents' sequencing data for signal supporting an SV call.

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
    ├── SAMPLE_REF.SVs.dicast.tsv                     TSV with all variant scores
    ├── SAMPLE_CALLER.dicast.vcf                      per input caller, DQ-tagged
    └── SAMPLE_REF.SVs.dicast.merged.vcf              merged best-per-cluster VCF
```

`multi` builds the exact same tree per sample, under `WORKDIR/SAMPLE/...`.

The scores TSV lists all input variants, each assigned an individual dicast
quality score. The per-caller VCF files are the input VCFs re-emitted with an
additional INFO tag `DQ` carrying the dicast quality score. The merged VCF
first builds an overlap graph of SV calls likely representing the same
variant, then uses the dicast score to pick one representative variant per
cluster: essentially a deduplicated set of scored calls based on the input
VCFs.
