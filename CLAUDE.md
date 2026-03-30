# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Dicast is a Python tool for structural variant (SV) confidence scoring using XGBoost models. It extracts genomic features from reference files and BAM alignments, then predicts SV quality scores.

## Environment Setup

- Conda environment: `conda activate svlearn` (must be active before running anything)
- Python 3.10, dependencies defined in `environment.yml`

## Running

```bash
# Single sample
python3 dicast.py call --bam <bam> --vcf <vcf> [options]

# Cohort mode
python3 dicast.py cohort --vcf <vcf> [options]
```

See `dicast.sh` and `dicast_cohort.sh` for full usage examples.

## Architecture

- `dicast.py` — main entry point (~20K lines, orchestrates the full pipeline)
- `dicast_lib/` — core modules: parsing, model, prepare, collect_reference, collect_illumina, cohort, evaluate, utils
- `configs/` — YAML feature configuration files
- `models/deployment/` — trained XGBoost models (JSON format, not pickle)

Two modes: **call** (single-sample feature extraction + prediction) and **cohort** (multi-sample analysis with family filtering).

Pipeline stages (call mode): variant preparation -> reference feature collection -> alignment feature collection (parallelized with joblib) -> feature combination -> XGBoost prediction -> VCF annotation with DQ info tag.

## Key Gotchas

- Output format is CSV (recently changed from TSV)
- Models must be JSON format in `models/deployment/`
- Sample IDs can be integers — handle type coercion carefully
- URL-encoded sample lists in VCF cohort fields (`%2C` -> `,`) need decoding
- Supported SV types: DEL, DUP, INS, INV (BND/translocations have special handling)
- The notebook `eva.ipynb` references old import paths (`lib.evaluate` instead of `dicast_lib.evaluate`)

## Obsidian

- company: lucid
- project: dicast
- tag: #lucid
- todoist_project: Lucid Genomics
- todoist_section: Dicast
