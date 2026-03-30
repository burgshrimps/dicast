# Development Guide

## Training Workflow

Source: `Dicast.train()` in `dicast_lib/model.py`

### Requirements

Training requires a DataFrame with:
- All feature columns (variant + reference + alignment, see [features.md](features.md))
- A `confirmation_status` column (1 = true positive, 0 = false positive)
- Columns: `chrom`, `cohort`, `sample`, `reference`, `caller`

### Steps

```python
from dicast_lib.model import Dicast

dicast = Dicast('DEL', feature_config_file='configs/config.yaml')
dicast.load_from_df(training_df)        # or load_from_csv()
dicast.impute_missing_values()
dicast.train(
    model_type='XGBoost',
    model_params={'n_estimators': 100, 'random_state': 42, 'n_jobs': 1},
    chroms=['chr1', 'chr2', ...]        # optional: restrict to specific chromosomes
)
dicast.save('models/deployment/dicast_DEL.json')
```

`save()` writes both the model JSON and a `*_metadata.json` with training details (see [models.md](models.md)).

### Feature Importance

```python
importance_df = dicast.get_feature_importance()
# Returns DataFrame with 'feature' and 'importance' columns, sorted descending
```

## Evaluation (`Eva` class)

Source: `dicast_lib/evaluate.py`

The `Eva` class compares dicast predictions against a benchmark (ground truth) VCF and other SV caller methods.

### Initialization

```python
from dicast_lib.evaluate import Eva

params = {
    'sample': 'HG002',
    'ref': 'hg38',
    'dicast': '/path/to/dicast_predictions.tsv',
    'dicast_ref_annot': '/path/to/ref_annotated.tsv',
    'benchmark': '/path/to/benchmark.tsv',
    'curation_root': '/path/to/curation/',
    'curation_date': '20240115',
    'vcf': {
        'ill': {
            'delly': '/path/to/delly.vcf',
            'manta': '/path/to/manta.vcf',
            ...
        }
    }
}
eva = Eva(params, params_ref)
```

### Key Methods

- `read_method_variants()` -- reads all caller VCFs
- `read_benchmark_variants()` -- reads ground truth TSV
- `read_dicast_variants()` -- reads dicast prediction TSV
- Overlap computation uses bioframe with configurable thresholds:
  - `max_dist_overlap = 500` (bp distance for INS)
  - `min_size_overlap = 0.7` (reciprocal overlap for DEL/DUP/INV)

### Notebook

`eva.ipynb` is the evaluation notebook. **Known issue**: it references the old import path `lib.evaluate` instead of `dicast_lib.evaluate`. Update the import before using:

```python
# Old (broken):
from lib.evaluate import Eva
# New (correct):
from dicast_lib.evaluate import Eva
```

## Configuration Variants

See [features.md](features.md) for full config documentation. The config variants are useful for ablation experiments:

| Config | Use Case |
|---|---|
| `configs/config.yaml` | Production: all features |
| `configs/config_noref_novar.yaml` | Ablation: alignment + GC only |
| `configs/config_noref_novar_onlycov.yaml` | Ablation: coverage + GC only |

Pass to `Dicast()` constructor: `Dicast(sv_type, feature_config_file='configs/config.yaml')`

Note: the scoring functions in `dicast.py` (`score_variants()`) do **not** pass a config file, so they use the full feature set by default.

## Legacy README

The `README.md` documents 5 legacy modes (prepare, train, predict, test, manual curation) that no longer exist in the CLI. The codebase has been refactored to use `call` and `cohort` subcommands. The README is retained for historical reference but is not current.

## Key Dependencies

From `environment.yml`:

| Package | Version | Purpose |
|---|---|---|
| xgboost | 1.7.4 | Model training and prediction |
| pysam | 0.18.0 | BAM/VCF reading |
| bioframe | 0.3.2 | Genomic interval operations (overlap detection) |
| networkx | 2.8.4 | Graph-based variant clustering |
| pyBigWig | 0.3.18 | GC content BigWig file reading |
| pandas | 1.4.0 | Data manipulation throughout |
| numpy | 1.21.5 | Numerical operations |
| scikit-learn | 1.0.2 | Metrics (precision, recall, ROC) |
| joblib | 1.1.0 | Parallel alignment feature collection |
| vcfpy | (via pip) | VCF writing with DQ INFO tag |
| plotly | 5.5.0 | Evaluation visualizations |
| pyyaml | 6.0 | Config file parsing |

## Utility Functions

Source: `dicast_lib/utils.py`

Key functions used across the codebase:
- `replace_filename(filename, sample, ref)` -- substitutes SAMPLE/REF placeholders in file paths
- `caller_vcf_to_dataframe()` -- parses a single-caller VCF into a DataFrame
- `sample_vcf_to_dataframe()` -- parses a sample-level VCF (cohort mode)
- `parse_vcf()` -- VCF parsing for evaluation
- `mad()` -- median absolute deviation calculation (used in baseline statistics)
- `read_parameters()` -- JSON parameter file reader (legacy)
