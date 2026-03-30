# Feature System

Source: `dicast_lib/model.py` (feature definitions), `dicast_lib/collect_illumina.py` (bin system), `configs/*.yaml` (toggles)

## Feature Categories

### 1. Variant Features (1 feature)

| Feature | Description |
|---|---|
| `sv_len` | Length of the structural variant |

### 2. Reference Features (23 features)

Collected by `ReferenceAnnotator` in `dicast_lib/collect_reference.py`.

| Feature | Source File | Description |
|---|---|---|
| `rep_LINE` | RepeatMasker TSV | LINE repeat overlap |
| `rep_SINE` | RepeatMasker TSV | SINE repeat overlap |
| `rep_LTR` | RepeatMasker TSV | LTR repeat overlap |
| `rep_DNA` | RepeatMasker TSV | DNA repeat overlap |
| `rep_Simple_repeat` | RepeatMasker TSV | Simple repeat overlap |
| `rep_Satellite` | RepeatMasker TSV | Satellite repeat overlap |
| `rep_Low_complexity` | RepeatMasker TSV | Low complexity region overlap |
| `rep_Retroposon` | RepeatMasker TSV | Retroposon overlap |
| `rep_snRNA` | RepeatMasker TSV | snRNA overlap |
| `rep_tRNA` | RepeatMasker TSV | tRNA overlap |
| `rep_srpRNA` | RepeatMasker TSV | srpRNA overlap |
| `rep_rRNA` | RepeatMasker TSV | rRNA overlap |
| `rep_RC` | RepeatMasker TSV | Rolling circle repeat overlap |
| `rep_scRNA` | RepeatMasker TSV | scRNA overlap |
| `rep_RNA` | RepeatMasker TSV | RNA overlap |
| `rep_VNTR` | Chaisson BED | VNTR region overlap |
| `rep_STR` | Chaisson BED | STR region overlap |
| `cpg_islands` | TSV | CpG island overlap |
| `centromeres` | TSV | Centromere overlap |
| `asmb_gaps` | TSV | Assembly gap overlap |
| `alt_haps` | TSV | Alternative haplotype overlap |
| `GC_content_left` | BigWig | GC content left of variant |
| `GC_content_right` | BigWig | GC content right of variant |

### 3. Alignment Features (variable count, SV-type dependent)

Collected by `AlignmentAnnotatorIllumina` in `dicast_lib/collect_illumina.py`.

## Bin System

The bin system defines genomic windows around and within SVs where alignment features are measured. All bins are 50bp wide (offset range [-52, -2] or [+2, +52] relative to a position).

### Breakpoint Bins

**DEL, DUP, INV** -- 4 bins:
```
  I: [start-52, start-2]     (upstream of start breakpoint)
  II: [start+2, start+52]    (downstream of start breakpoint)
  III: [end-52, end-2]        (upstream of end breakpoint)
  IV: [end+2, end+52]         (downstream of end breakpoint)
```

**INS** -- 2 bins (no end breakpoint):
```
  I: [start-52, start-2]
  II: [start+2, start+52]
```

**BND** -- 4 bins across two chromosomes:
```
  I: [chrom:start-52, chrom:start-2]
  II: [chrom:start+2, chrom:start+52]
  III: [chrom_2:end-52, chrom_2:end-2]
  IV: [chrom_2:end+2, chrom_2:end+52]
```

### Body Bins (DEL, DUP, INV only)

The SV body (between start+52 and end-52) is divided into 4 quarters:
```
  IIa:  [start+52, body_I]     (1st quarter)
  IIb:  [body_I, body_II]      (2nd quarter)
  IIIb: [body_II, body_III]    (3rd quarter)
  IIIa: [body_III, end-52]     (4th quarter)
```

INS and BND have no body bins.

### Connection Features

Connections measure relationships between pairs of breakpoint bins. Pairs follow this pattern:
```
  I-II, I-III, I-IV, II-III, II-IV, III-IV
```

INS has no connection features (only 2 bins, and connections require the full set).

## Alignment Feature Metrics

### Breakpoint Metrics (per bin)

| Metric | Description |
|---|---|
| `ill_cov_mean` | Mean log2 fold-change coverage vs baseline |
| `ill_cov_std` | Std of log2 fold-change coverage vs baseline |
| `ill_isize_mean` | Mean insert size deviation from baseline |
| `ill_isize_std` | Std of insert size deviation from baseline |
| `ill_mapq_mean` | Mean mapping quality deviation from baseline |
| `ill_mapq_std` | Std of mapping quality deviation from baseline |
| `ill_clipreads` | Fraction of clipped reads |
| `ill_splitreads` | Fraction of split reads |
| `ill_disco_ff` | Discordant pairs (forward-forward) |
| `ill_disco_rr` | Discordant pairs (reverse-reverse) |
| `ill_disco_rf` | Discordant pairs (reverse-forward) |
| `ill_disco_tx` | Discordant pairs (translocation) |

### Body Metrics (per body bin)

| Metric | Description |
|---|---|
| `ill_cov_mean` | Mean log2 fold-change coverage vs baseline |
| `ill_cov_std` | Std of log2 fold-change coverage vs baseline |

### Connection Metrics (per bin pair)

| Metric | Description |
|---|---|
| `ill_splitreads` | Split reads spanning both bins |
| `ill_disco_ff` | Discordant pairs (forward-forward) between bins |
| `ill_disco_rr` | Discordant pairs (reverse-reverse) between bins |
| `ill_disco_rf` | Discordant pairs (reverse-forward) between bins |

## Feature Naming Convention

Features are named by combining metric + bin suffix:

- **Breakpoint**: `{metric}_{bin}` -- e.g., `ill_cov_mean_I`, `ill_splitreads_IV`
- **Body**: `{metric}_{body_bin}` -- e.g., `ill_cov_mean_IIa`, `ill_cov_std_IIIa`
- **Connection**: `{metric}_{binA}_{binB}` -- e.g., `ill_splitreads_I_II`, `ill_disco_rf_II_IV`

## Feature Counts by SV Type

| SV Type | BP Bins | Body Bins | Connections | Total Alignment Features |
|---|---|---|---|---|
| DEL | 12 x 4 = 48 | 2 x 4 = 8 | 4 x 6 = 24 | 80 |
| DUP | 12 x 4 = 48 | 2 x 4 = 8 | 4 x 6 = 24 | 80 |
| INV | 12 x 4 = 48 | 2 x 4 = 8 | 4 x 6 = 24 | 80 |
| INS | 12 x 2 = 24 | 0 | 0 | 24 |

Total features = 1 (variant) + 23 (reference) + alignment = 104 (DEL/DUP/INV) or 48 (INS).

## Exome Differences (`DicastExome`)

Source: `dicast_lib/model_exome.py`

Exome models use fewer alignment features:
- **No coverage features** in breakpoint bins (`ill_cov_mean`, `ill_cov_std` removed)
- **No body bins** at all (empty for all SV types)
- **No `ill_disco_tx`** in breakpoint features
- Connection features remain the same

This results in 9 breakpoint metrics (vs 12 for WGS).

## Configuration Files

Source: `configs/*.yaml`

Config YAML has 5 sections with binary (0/1) toggles:

```yaml
variant:
  sv_len: 1
reference:
  rep_LINE: 1
  # ... (23 features)
alignment:
  ill_cov_mean: 1
  # ... (12 metrics)
alignment_body:
  ill_cov_mean: 1
  ill_cov_std: 1
alignment_conn:
  ill_splitreads: 1
  ill_disco_ff: 1
  ill_disco_rr: 1
  ill_disco_rf: 1
```

### Config Variants

| File | Description |
|---|---|
| `configs/config.yaml` | All features enabled (default) |
| `configs/config_noref_novar.yaml` | Alignment features + GC content only (no variant/reference) |
| `configs/config_noref_novar_onlycov.yaml` | Coverage + GC content only |

Config variants are useful for ablation studies and feature importance analysis.

## Imputation Rules

Source: `Dicast.impute_missing_values()` in `dicast_lib/model.py`

| Feature | Imputation | Reason |
|---|---|---|
| `GC_content_left` | Median fill | BigWig may not cover all positions |
| `GC_content_right` | Median fill | BigWig may not cover all positions |
| `ill_cov_mean_{I,II,III,IV}` | Fill with `cov_thr` (5) | Coverage exceeded threshold during feature extraction |
| `ill_cov_std_{I,II,III,IV}` | Fill with 1 | Coverage exceeded threshold during feature extraction |

The `cov_thr = 5` represents the log2 threshold where feature collection was aborted (extremely high coverage).
