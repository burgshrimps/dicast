# Cohort Mode

Source: `dicast_lib/cohort.py`, cohort section of `dicast.py`

## Purpose

Regenotype structural variants across multiple samples in a cohort. Detects variants that individual callers missed in some samples by re-evaluating alignment evidence using the Dicast model.

## Input Formats

### CSV File (`--csv`)

Required columns:

| Column | Description |
|---|---|
| `ID` | Variant identifier |
| `SAMPLE` | Sample name |
| `TYPE` | SV type (DEL, DUP, INS, INV) |
| `CHR` | Chromosome |
| `START` | Start position |
| `END` | End position |
| `SIZE` | SV length |
| `FILTER` | Filter status |
| `QUAL` | Quality score (dicast quality from previous run) |
| `GT` | Genotype (e.g., `0/1`) |
| `COHORT_AC` | Cohort allele count |
| `COHORT_SC` | Cohort sample count |
| `COHORT_SUP_SAMPLES` | Supporting samples as list of dicts |

The `COHORT_SUP_SAMPLES` column uses this format:
```
[{'id': 'HG002', 'gt': '0/1'}, {'id': 'HG003', 'gt': '0/1'}]
```

This is parsed via `ast.literal_eval()` in `Cohort._parse_cohort_samples()`.

### BAM Files (`--bams`)

Space-separated paths. Sample name is extracted from the path:
```python
bam_dict = {bam_path.split('/')[-2]: bam_path for bam_path in arguments.bams}
```

The parent directory name of each BAM must match the sample name in the CSV.

### PED File (`--ped`)

TSV with header row:

| Name | Family | Mother | Father |
|---|---|---|---|
| HG002 | GIAB | HG004 | HG003 |
| HG003 | GIAB | 0 | 0 |
| HG004 | GIAB | 0 | 0 |

Used to build `family_dict`: maps each sample to other samples in the same family.

## Pipeline Flow

### 1. Single-Sample Shortcut

If only 1 sample in the CSV, the file is simply copied to `{csv_path}.regenotyped.csv` with no processing.

### 2. Variant Preparation

```python
VP = VariantPrep(cohort, ref, workdir, technology, chroms, fai, sv_types, mode='cohort')
VP.read_csv(df)
VP.reformat_csv(cohort, technology, ref)
```

`reformat_csv()` maps CSV columns to internal format:
```
ID -> id, SAMPLE -> sample, TYPE -> sv_type, CHR -> chrom, START -> start,
END -> end, SIZE -> sv_len, FILTER -> filter, QUAL -> qual, GT -> genotype,
COHORT_AC -> cohort_ac, COHORT_SC -> cohort_sc, COHORT_SUP_SAMPLES -> cohort_samples
```

Adds: `cohort`, `technology`, `reference`, `caller='dicast'`, `chrom_2=NaN`.

### 3. Filtering

Two-stage filtering:

1. **Standard filters** (`filter_variants()`):
   - Out-of-bounds variants (50bp padding from chromosome ends)
   - Unsupported SV types (keeps DEL, DUP, INS, INV only)

2. **Cohort-specific filters** (`filter_variants_cohort(num_samples)`):
   - Variants present in all samples (`cohort_sc >= num_samples`) -- nothing to regenotype
   - Variants with `qual < 0.1` -- too low quality to consider

Both filtered and unfiltered DataFrames are saved:
- `{workdir}/{cohort}_{ref}.SVs.raw.tsv` (filtered)
- `{workdir}/{cohort}_{ref}.SVs.raw.unfiltered.tsv` (unfiltered, used for overlap detection)

### 4. Missing Variant Detection (`Cohort.get_missing_variants()`)

For each sample, finds variants where the sample is **not** in `cohort_samples`:

```python
mask = ~cohort_df['cohort_samples'].apply(
    lambda x: any(item.get('id') == sample for item in self._parse_cohort_samples(x))
)
```

With `--filter-fam`: additionally requires that at least one family member is in `cohort_samples`. If a sample has no family members, no variants are considered missing.

Missing variants are saved as per-sample TSV files: `{workdir}/{sample}_{ref}.SVs.raw.tsv`

### 5. Per-Sample Feature Collection + Scoring

Same as call mode (see [architecture.md](architecture.md)):
1. Reference features -> `{sample}_{ref}.SVs.ref.tsv`
2. Alignment features (parallelized) -> `{sample}_{ref}.SVs.aln.ill.{chrom}.{svtype}.tsv`
3. Feature combination -> `{sample}_{ref}.SVs.annot.tsv`
4. Scoring -> `{sample}_{ref}.SVs.dicast.tsv`

### 6. Cohort Update (`Cohort.update_cohort_information()`)

Loads all per-sample dicast predictions. For variants scoring above threshold:
- Adds sample to `cohort_samples`
- Increments `cohort_ac` and `cohort_sc`
- Assigns genotype `0/1`

**Threshold**: `dicast_thr = 0.4`

### 7. Overlap Detection (`Cohort.find_overlapping_variants()`)

Processes each sample separately. Combines the sample's original variants with newly regenotyped variants, then clusters overlapping ones.

**Overlap criteria by SV type:**

| SV Type | Method | Threshold |
|---|---|---|
| DEL, DUP, INV | Reciprocal overlap (via bioframe `closest`, k=10) | >= 0.5 (50%) |
| INS | Breakpoint distance (via bioframe `closest`, k=10) | < 200bp |

**Clustering**: Uses networkx connected components on the overlap graph. Each overlapping pair becomes an edge; connected components form clusters.

**Best variant selection** (per cluster, per sample):
1. Highest `cohort_ac`
2. Tiebreak: highest `dicast_qual`

All non-best variants in a cluster are added to the sample's blacklist.

### 8. CSV Update (`Cohort.update_csv_file()`)

Two passes:

1. **Update existing rows**: for variants in `variant_cohort_map`, update `COHORT_AC`, `COHORT_SC`, `COHORT_SUP_SAMPLES`
2. **Add new rows**: for regenotyped variants above threshold that aren't already in the CSV for that sample, create new rows with:
   - `QUAL` = dicast_qual
   - `FILTER` = `['PASS']`
   - `GT` = `0/1`
   - `ID` = `{SAMPLE}.MERGED.{TYPE}.{CHR_num}.{START}.{SIZE}`
   - `NUM_SUPP_CALLERS` = 1
   - `DICAST` = True, all other callers (DELLY, MANTA, LUMPY, GRIDSS, CNVNATOR, SNIFFLES) = False

Finally, for each row, the current sample is removed from `COHORT_SUP_SAMPLES` (supporting samples should not include self).

Output: `{original_csv}.regenotyped.csv`, sorted by SAMPLE and ID.
