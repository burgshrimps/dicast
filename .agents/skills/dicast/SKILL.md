---
name: dicast
description: Score structural variant (SV) calls from short-read BAM + caller VCFs with dicast, and check/fix input VCFs before scoring. Use when asked to run, score, filter, or triage SV calls (DEL, DUP, INS) from callers such as delly, manta, lumpy/smoove, gridss, cnvnator, or dysgu; to check whether an SV VCF is in a format dicast understands (dicast check); or to fix a VCF dicast rejects or reports zero usable records for (non-chr contigs, non-hg38 reference, a multi-sample cohort VCF, gridss/delly breakend-only output). dicast is not an SV caller -- it re-scores calls already made by other tools -- and it does not pair breakends or score INV/CNV/TRA calls.
---

# dicast: score structural variant calls

dicast takes SV calls already made by one or more callers plus the BAM they
were called from, and assigns each call a confidence score (`dicast_qual`,
0-1) using a pretrained XGBoost model per SV type (DEL, DUP, INS). It is
caller-agnostic and works best combining several callers' raw, unfiltered
output for the same sample.

## 1. Preflight

```bash
dicast --help
```

If the command is missing, install it:

```bash
git clone https://github.com/burgshrimps/dicast.git && cd dicast
conda env create -f environment.yml && conda activate dicast
# or: pip install git+https://github.com/burgshrimps/dicast.git
```

On first `call`/`multi` run, dicast asks where to store the hg38 annotation
files it downloads once (~2.1 GB; press Enter to accept the suggested
location). `--annot-dir` points it at an existing copy instead.

Also needed before running anything:

- The BAM (`--bam` / `--bams`) needs an index (`.bai`/`.csi`) next to it.
- The BAM and the `--fai` must both be `chr`-named hg38/GRCh38
  (`chr1..chr22, chrX, ...`); dicast's alignment-feature collection indexes
  the BAM by chromosome name and assumes that layout. Anything else is
  unsupported -- see "Fixing inputs" below for the input-VCF side of this
  (contig naming, non-hg38 FAI), which `dicast check` catches without
  touching the BAM at all.

## 2. Check inputs first

Before running `call`/`multi`, always run `dicast check` on the input VCFs.
It reads every file the same way `call`/`multi` would, but only reports what
it understood -- no workdir, models, annotations or BAM required:

```bash
dicast check \
    --vcfs delly=/path/delly.vcf.gz manta=/path/manta.vcf.gz \
    --fai /path/to/reference.fa.fai \
    --sample SAMPLE_NAME
```

Add `--out report.tsv` to also get the same numbers as a fixed-column TSV
(one row per file), and `--chrom` to restrict which contigs count (defaults
to `chr1..chr22, chrX`, same as `call`/`multi`).

Example output, one summary row per file plus indented detail:

```
CALLER  STATUS           READ  USABLE  PATH
-------------------------------------------
delly   OK                 20      20  demo_delly.vcf.gz
    kept:    DEL=10, DUP=0, INS=10
    contigs: as-is (1/1 mapped)
    sample:  demo
gridss  ZERO USABLE       186       0  gridss_raw.vcf.gz
    kept:    DEL=0, DUP=0, INS=0
    dropped: BND=186
    contigs: as-is (1/1 mapped)
    sample:  caller.bam
    hint: this file looks breakend-only (every record dropped as BND). dicast does not pair breakends itself -- for gridss, annotate the raw VCF with scripts/gridss_simple_event_annotation.R first; for delly, use its DEL/DUP/INS calls instead of TRA/BND.
```

How to read it:

- **STATUS** is `OK`, `ZERO USABLE` (file parsed but nothing survived --
  `call`/`multi` would exit 1 on this exact input), or `ERROR` (the file
  itself, or the requested `--sample`, could not be resolved at all; see the
  printed `error:` line).
- **kept** breaks usable records down by SV type (DEL/DUP/INS -- the only
  types dicast's models cover).
- **dropped** breaks unusable records down by reason (`contig, no_svtype,
  INV, BND, CNV, other_svtype, small_variant, missing_end, bnd_mate,
  absent_in_sample`) -- see
  [reference/drop-reasons.md](reference/drop-reasons.md) for what each one
  means and how to fix it.
- **contigs** shows how contig names in the file were mapped onto the FAI's
  names: `as-is`, `+chr` (a `chr` prefix was added), or `-chr` (one was
  stripped), plus how many of the file's distinct contigs mapped at all.
- **sample** is the sample column dicast actually read genotypes from (only
  shown when the file has samples); a `warning:` line fires if it silently
  picked the file's one sample column over a different `--sample` you asked
  for.
- A `hint:` line fires specifically when a file is breakend-only (every
  record dropped as `BND`) -- see "Fixing inputs" below.

Fix whatever `dicast check` flags, re-run it, and only move on to
`call`/`multi` once every file you plan to use is `OK` (or an intentionally
partial `ZERO USABLE` file you are dropping from the run).

## 3. Fixing inputs

**Wrong contig naming (`dropped: contig`, or `STATUS ERROR` naming the
file's contigs):** dicast auto-maps `chr1` <-> `1` per contig, so this only
fires when a file mixes styles inconsistently or uses names the FAI doesn't
have at all (e.g. GRCh37 `MT` against an hg38 FAI's `chrM`). Rename with
`bcftools`:

```bash
# one line per old->new name, e.g. "1\tchr1"
bcftools annotate --rename-chrs chr_map.txt in.vcf.gz -Oz -o out.vcf.gz
```

**Cohort/multi-sample VCF, wrong sample picked or `STATUS ERROR` listing
sample names:** dicast's `--sample` only matters when a file has more than
one sample column; if the wrong one gets picked, or `--sample` doesn't match
any column, extract just that sample first:

```bash
bcftools view -s SAMPLE_NAME in.vcf.gz -Oz -o out.vcf.gz
```

**gridss raw output (`STATUS ZERO USABLE`, `dropped: BND`, breakend-only
hint):** gridss's raw VCF is all breakend (`BND`) records; dicast does not
pair breakends itself. Run the event-simplification script dicast ships
(`scripts/gridss_simple_event_annotation.R`, adapted from gridss's own
`example/simple-event-annotation.R`; needs R with the Bioconductor packages
VariantAnnotation and StructuralVariantAnnotation, e.g. via conda
`bioconductor-variantannotation bioconductor-structuralvariantannotation r-stringr`).
It adds `INFO/SIMPLE_TYPE` (DEL/DUP/INS/...) and `INFO/SVLEN`, which dicast
reads in preference to `SVTYPE=BND`:

```bash
Rscript scripts/gridss_simple_event_annotation.R gridss_raw.vcf gridss_simple.vcf
```

**cnvnator output:** cnvnator's native output isn't VCF at all; convert with
the script it ships before handing dicast the result:

```bash
cnvnator2VCF.pl -prefix study -reference hg38 calls.txt > cnvnator.vcf
```

**delly TRA/BND-only output (`STATUS ZERO USABLE`, `dropped: BND`,
breakend-only hint):** delly's translocation (`TRA`) calls are breakend
pairs like gridss's, and dicast drops them the same way. Request delly's
`DEL`/`DUP`/`INS` calls instead (delly's default `call` mode already
produces these alongside `TRA`; just don't filter to `SVTYPE=TRA` before
handing the file to dicast).

**Everything else (manta, lumpy/smoove, dysgu, and any other
sequence-resolved or symbolic-ALT caller):** these work as-is, no
pre-processing needed. dicast reads `INFO/SVTYPE`/`INFO/END`/`INFO/SVLEN`,
symbolic ALTs (`<DUP:TANDEM>`, `<DEL:ME:ALU>`, ...), and sequence-resolved
REF/ALT alleles (`>=50 bp` length difference) as sources of SV type and
length, in that order.

**Non-hg38 input:** unsupported. `dicast check` (via the same FAI check
`call`/`multi` runs) warns loudly -- `chr1..chr22`/`chrX` missing from the
FAI, or `chr1`'s length not matching hg38/GRCh38 -- but keeps going rather
than erroring, so do not ignore that warning; re-align/re-call against hg38
first.

## 4. Commands

Single sample (`call`), after `dicast check` is clean for every input:

```bash
dicast call \
    --sample SAMPLE_NAME \
    --workdir WORKDIR \
    --fai /path/to/reference.fa.fai \
    --bam /path/to/sample.bam \
    --vcfs delly=/path/delly.vcf.gz manta=/path/manta.vcf.gz \
    --threads 24
```

Multiple samples with cross-sample rescue (`multi`, e.g. a trio -- checks
each sample's other relatives' BAM for signal at variants only one of them
called):

```bash
dicast multi \
    --bams MOTHER=/path/mother.bam FATHER=/path/father.bam CHILD=/path/child.bam \
    --vcfs MOTHER:delly=/path/mother_delly.vcf.gz MOTHER:manta=/path/mother_manta.vcf.gz \
           FATHER:delly=/path/father_delly.vcf.gz FATHER:manta=/path/father_manta.vcf.gz \
           CHILD:delly=/path/child_delly.vcf.gz  CHILD:manta=/path/child_manta.vcf.gz \
    --workdir WORKDIR --fai /path/to/reference.fa.fai \
    --threads 24
```

Both accept `--chrom` (default: all of `chr1..chr22, chrX`) and `--sv_types`
(default: `DEL DUP INS`) to restrict the run, and `--pop` to add the PAV
population catalog as an extra pseudo-caller and switch to
population-aware models for DEL/INS (downloaded automatically like the
other annotations, unless `--pop-catalog` points at your own copy).

## 5. Outputs and thresholds

`call` fills `--workdir` (and `multi` fills `--workdir/SAMPLE`) with, most
relevantly: `output/SAMPLE_REF.SVs.dicast.tsv` (every input variant with its
`dicast_qual` score), one `output/SAMPLE_CALLER.dicast.vcf` per input caller
(same VCF, `INFO/DQ` added), and `output/SAMPLE_REF.SVs.dicast.merged.vcf`
(one deduplicated, best-scoring call per cluster across callers). Run
`dicast check --out report.tsv` beforehand if you also want a machine-
readable record of what was dropped from the inputs.

Recommended thresholds for calling a scored variant a true positive:

| SV type | threshold |
|---|---|
| DEL | 0.45 |
| INS | 0.30 |
| DUP | 0.40 |

To visually inspect the read-level evidence behind a specific call (does
the coverage/read-pair/alignment signal actually support it?), pair dicast
with its companion tool [cuban](https://github.com/burgshrimps/cuban) -- it
renders one PNG per call from a merged VCF plus the sample's BAM:
`cuban --vcf output/SAMPLE_REF.SVs.dicast.merged.vcf --sample SAMPLE_NAME:/path/sample.bam --outdir plots/`.
