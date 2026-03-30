# Models

Source: `dicast_lib/model.py`, `dicast_lib/model_exome.py`

## Deployment Models

Located in `models/deployment/`:

### WGS Models (XGBoost JSON)

| File | SV Type |
|---|---|
| `dicast_DEL.json` | Deletions |
| `dicast_DUP.json` | Duplications |
| `dicast_INS.json` | Insertions |

Loaded via:
```python
model = xgb.XGBClassifier()
model._Booster = xgb.Booster()
model._Booster.load_model(model_filename)
```

### Exome Models (Pickle)

| File | SV Type |
|---|---|
| `dicast_exome_DEL.pkl` | Deletions |
| `dicast_exome_DUP.pkl` | Duplications |

Loaded via `pickle.load()`. Uses `DicastExome` class which has fewer alignment features (no coverage, no `ill_disco_tx`).

### Metadata Files

Each model has a companion `*_metadata.json`:

```json
{
    "date": "2024-01-15 10:30:00",
    "model_type": "XGBoost",
    "model_params": {
        "n_estimators": 100,
        "random_state": 42,
        "n_jobs": 1
    },
    "features": ["sv_len", "rep_LINE", ...],
    "sv_type": "DEL",
    "reference": "hg38",
    "cohorts": ["tgenvar"],
    "samples": ["HG002", "HG003", ...],
    "callers": ["delly", "manta", "lumpy", ...],
    "chroms_train": ["chr1", "chr2", ...],
    "number_variants": 12345,
    "number_variants_positive": 6789,
    "number_variants_negative": 5556
}
```

## Heuristic Scoring (No Model)

### Inversions (`score_inversions`)

INV variants have no trained model. Instead, a rule-based heuristic assigns binary 0/1 scores. A variant scores 1 (pass) if **all** conditions are met:

| Condition | Threshold |
|---|---|
| Clipped reads | `ill_clipreads_{I,II,III,IV} > 0.2` in >= 3 of 4 bins |
| Discordant reads (relaxed) | `ill_disco_{ff,rr}_{I,II,III,IV} > 0.1` in >= 3 of 8 columns |
| Discordant reads (strict) | `ill_disco_{ff,rr}_{I,II,III,IV} > 0.2` in >= 2 of 8 columns |
| Coverage | `ill_cov_mean_{I,II,III,IV} <= 3.5` in all 4 bins |
| Length | `sv_len < 3,000,000` |

### Translocations/BND (`score_translocations`)

BND variants also use a heuristic. Scores 1 if **all** conditions are met:

| Condition | Threshold |
|---|---|
| Clipped reads | `> 0.2` in >= 2 of 4 bins |
| Coverage | `<= 3` in all 4 bins |
| Discordant (inv pattern) | `ill_disco_{ff,rr} > 0.2` in >= 2 of 8 columns |
| Discordant (dup pattern) | `ill_disco_rf > 0.3` in >= 2 of 4 columns |
| Discordant (tra pattern) | `ill_disco_tx > 0.3` in >= 2 of 4 columns |
| Mapping quality | `ill_mapq_mean >= -0.5` in all 4 bins |
| Split reads | `> 0.1` in >= 2 of 4 bins |

Discordant conditions use OR logic: at least one of inv/dup/tra patterns must match.

## Model Prediction Flow

```python
dicast = Dicast(sv_type)                     # Initialize with SV type
dicast.load_from_csv(variant_features_file)  # Load annotated variants
dicast.impute_missing_values()               # Handle NAs (see features.md)
dicast.load(model_filename)                  # Load XGBoost JSON model
dicast.predict()                             # predict_proba[:, 1] -> dicast_qual
df_output = dicast.to_df()                   # Export prediction columns
```

Output columns from `to_df()`:
```
id, cohort, sample, reference, technology, caller, sv_type,
chrom, chrom_2, start, end, sv_len, filter, qual, dicast_qual, genotype
```

## XGBoost Version Compatibility

The code includes a compatibility check for xgboost 1.7.4:
```python
if not hasattr(self.model, 'n_classes_'):
    self.model.n_classes_ = 2
```

This is needed because older model files may not have the `n_classes_` attribute set on the loaded classifier.

## Training (`Dicast.train`)

Training requires a DataFrame with a `confirmation_status` column (binary: 1 = true positive, 0 = false positive).

```python
dicast = Dicast(sv_type, feature_config_file='configs/config.yaml')
dicast.load_from_df(training_df)
dicast.impute_missing_values()
dicast.train(model_type='XGBoost', model_params={...}, chroms=[...])
dicast.save('models/deployment/dicast_DEL.json')  # Saves model + metadata JSON
```

## Development Models

`models/development/` contains historical model versions (not used at runtime). Naming pattern: `dicast_{SVTYPE}_{YYYYMMDD}.json`, some with `_pop` suffix for population-level models.
