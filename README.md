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
`lumpy`, `gridss`, and `cnvnator`. Before a first run against new VCFs,
check them with [`dicast check`](#checking-inputs-dicast-check) -- see
[Input VCF requirements](#input-vcf-requirements) for what dicast expects
from a VCF and what it does automatically vs. what needs a fix first.

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

## Input VCF requirements

dicast is caller-agnostic: it needs `chrom, pos, id, sv_type, end, sv_len,
qual, filter, genotype` for each record and normalizes any reasonable VCF
dialect into that shape in one place before scoring. Specifically, it
auto-fixes:

- **SV type**, from `INFO/SVTYPE`, falling back to `INFO/SIMPLE_TYPE`
  (gridss's `simple-event-annotation.R` output), a symbolic ALT
  (`<DUP:TANDEM>`, `<DEL:ME:ALU>`, ...), or -- for sequence-resolved
  callers -- the REF/ALT length difference.
- **End position and length**, from `INFO/END` (or the REF allele of a
  sequence-resolved deletion), falling back to `INFO/SVLEN`; insertion
  lengths from `INFO/SVLEN` or the inserted sequence.
- **Contig naming**, mapping `chr1` <-> `1` per contig against your `--fai`.
- **The sample column**, used regardless of its name when a VCF has exactly
  one; for a multi-sample (cohort) VCF, `--sample` picks the column and only
  records where that sample carries an alt allele are kept.
- **Record IDs**, replaced internally with a globally unique
  `caller:ordinal` id; the VCF's own ID (if any) is preserved separately
  (`vcf_id` in `SVs.raw.tsv`/`SVs.dicast.tsv`, and as the VCF ID column of
  the merged output VCF -- see [Output](#output)).

What is **dropped** (and counted, never silently): inversions (`INV`),
breakends/translocations (`BND`, including delly `TRA`), copy-number-only
calls (`CNV`), any other SV type outside `{DEL, DUP, INS}`, and
sequence-resolved indels below 50 bp. A non-hg38 or non-`chr`-named FAI
produces a loud warning rather than an error. Run `dicast check` (below) to
see exactly what a given file will hit before running `call`/`multi` on it.

## Checking inputs: `dicast check`

`dicast check` reads a set of VCFs the same way `call`/`multi` would and
reports what it understood, without touching a workdir, models, annotations
or a BAM:

```bash
dicast check \
    --vcfs delly=/path/delly.vcf.gz manta=/path/manta.vcf.gz \
    --fai /path/to/reference.fa.fai \
    --sample SAMPLE_NAME
```

```
CALLER  STATUS           READ  USABLE  PATH
-------------------------------------------
delly   OK                 20      20  demo_delly.vcf.gz
    kept:    DEL=10, DUP=0, INS=10
    contigs: as-is (1/1 mapped)
    sample:  demo
```

Add `--out report.tsv` for the same numbers as a fixed-column TSV (one row
per file), and `--chrom` to restrict which chromosomes count (default:
`chr1..chr22, chrX`, same as `call`/`multi`). It exits 1 if any file ends up
with zero usable records -- the same failure mode `call`/`multi` would hit
on that input.

## Breakend-only VCFs (gridss, delly TRA)

dicast does not pair breakends itself, so a file whose records are all
`SVTYPE=BND` (raw gridss output, or delly's `TRA` calls) has every record
dropped and ends up with zero usable records. `dicast check` reports this
explicitly with a hint. Fixes:

- **gridss**: annotate the raw VCF with
  [`scripts/gridss_simple_event_annotation.R`](scripts/gridss_simple_event_annotation.R)
  (gridss's own `simple-event-annotation.R`, adapted to take file arguments
  and to work with current Bioconductor; needs R with
  `VariantAnnotation` and `StructuralVariantAnnotation`):

  ```bash
  Rscript scripts/gridss_simple_event_annotation.R gridss_raw.vcf gridss_simple.vcf
  ```

  This adds `INFO/SIMPLE_TYPE` (DEL/DUP/INS/...) and `INFO/SVLEN`, which
  dicast reads in preference to `SVTYPE=BND`.
- **delly**: use its `DEL`/`DUP`/`INS` calls (delly's default `call` mode
  already produces these) rather than filtering to `TRA`.

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

The scores TSV lists all input variants, each identified by dicast's
internal `id` (`caller:ordinal`, unique across every input file) and, where
the input VCF had one, the original `vcf_id`, and assigned an individual
dicast quality score. The per-caller VCF files are the input VCFs
re-emitted with an additional INFO tag `DQ` carrying the dicast quality
score, keyed to each record by `id`. The merged VCF first builds an overlap
graph of SV calls likely representing the same variant, then uses the
dicast score to pick one representative variant per cluster: essentially a
deduplicated set of scored calls based on the input VCFs. Each merged
record's VCF ID column carries the original caller ID (`vcf_id`) when the
input VCF had one, and always carries `INFO/DICAST_ID` (dicast's internal
`id`) so a record can be traced back to its row in the scores TSV either
way.

**To consider a call a true positive, we recommend the following score
thresholds:**

| SV type | threshold |
|---|---|
| DEL | 0.45 |
| INS | 0.30 |
| DUP | 0.40 |

## Use with AI coding agents

The repository ships an [Agent Skill](https://agentskills.io) that teaches
Claude Code and Codex how to run dicast, check and fix input VCFs with
`dicast check`, and interpret the output. Inside a clone of the repo both
tools pick it up automatically (Claude Code from `.claude/skills/dicast`,
Codex from `.agents/skills/dicast`), so you can ask for example *"check
whether these delly and manta VCFs are ready for dicast"* or *"run dicast on
this BAM and these caller VCFs and tell me which calls look real"*. To use
it from other projects, copy or symlink `.agents/skills/dicast` into
`~/.claude/skills/` (Claude Code) or `~/.codex/skills/` (Codex).
