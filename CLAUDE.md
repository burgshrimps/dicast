# CLAUDE.md

Dicast: SV confidence scoring using XGBoost. Extracts genomic features from reference files and BAM alignments, predicts quality scores.

## Environment

- `conda activate svlearn` (Python 3.10, deps in `environment.yml`)

## Running

```bash
python3 dicast.py call --sample <s> --bam <bam> --vcfs caller=path ... --models <dir> --workdir <dir>
python3 dicast.py cohort --cohort <name> --csv <csv> --bams <bam ...> --ped <ped> --models <dir> --workdir <dir> --fai <fai>
```

Full examples: `dicast.sh`, `dicast_cohort.sh`

## Architecture

- `dicast.py` — entry point, orchestrates call and cohort modes
- `dicast_lib/` — parsing, prepare, collect_reference, collect_illumina, model, model_exome, cohort, evaluate, utils
- `configs/` — YAML feature toggle files
- `models/deployment/` — production XGBoost models (.json) + exome models (.pkl)

## Key Gotchas

- Output format is CSV (recently changed from TSV)
- Models must be JSON format in `models/deployment/`
- Sample IDs can be integers — handle type coercion carefully
- URL-encoded sample lists in VCF cohort fields (`%2C` -> `,`) need decoding
- Supported SV types: DEL, DUP, INS, INV (BND/translocations have special handling)
- The notebook `eva.ipynb` references old import paths (`lib.evaluate` instead of `dicast_lib.evaluate`)

## Detailed Documentation

- [docs/architecture.md](docs/architecture.md) — pipeline stages, data flow, class relationships
- [docs/cli-reference.md](docs/cli-reference.md) — complete argument reference for call and cohort
- [docs/features.md](docs/features.md) — feature taxonomy, bin system, config YAML format
- [docs/models.md](docs/models.md) — model formats, deployment structure, INV/BND heuristic scoring
- [docs/cohort-mode.md](docs/cohort-mode.md) — cohort pipeline, PED, family filtering, overlap detection
- [docs/development.md](docs/development.md) — training, evaluation, notebooks, dependencies

## Git & Versioning

- **Workflow**: feature branch → PR → review → merge to `main`
- **Semantic versioning**: tag `main` after each merge. Bump patch (`0.1.1`) for fixes/cleanup, minor (`0.2.0`) for new features, major (`1.0.0`) for breaking changes.
- **Tagging**: `git tag X.Y.Z && git push origin --tags`. No `v` prefix.

## Obsidian

- company: lucid
- project: dicast
- tag: #lucid
- todoist_project: Lucid Genomics
- todoist_section: Dicast
