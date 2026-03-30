# Architecture

## Entry Point

`dicast.py` dispatches on `arguments.command`:
- `call` -- single-sample feature extraction + prediction
- `cohort` -- multi-sample regenotyping

All imports come from `dicast_lib/`: `parsing`, `prepare`, `collect_reference`, `collect_illumina`, `model`, `model_exome`, `cohort`, `evaluate`, `utils`.

## Call Mode Pipeline

### 1. Variant Preparation (`VariantPrep`)

```
dicast_lib/prepare.py
```

- Reads VCF files from multiple callers via pysam (format: `caller=path`)
- Each VCF is parsed by `caller_vcf_to_dataframe()` from `utils.py`
- Merges all caller DataFrames into one unified DataFrame
- Filters: removes out-of-bounds variants (checked against FAI chromosome sizes with 50bp padding) and unsupported SV types
- Output: `{workdir}/{sample}_{ref}.SVs.raw.tsv`

Constructor signature (single mode):
```python
VariantPrep(cohort, sample, ref, workdir, technology, chroms, fai, sv_types)
```

### 2. Reference Feature Collection (`ReferenceAnnotator`)

```
dicast_lib/collect_reference.py
```

Loads 8 reference annotation files and annotates each variant:

| Method | Reference File | Output Features |
|---|---|---|
| `annotate_repeats()` | RepeatMasker TSV | `rep_LINE`, `rep_SINE`, `rep_LTR`, `rep_DNA`, `rep_Simple_repeat`, `rep_Satellite`, `rep_Low_complexity`, `rep_Retroposon`, `rep_snRNA`, `rep_tRNA`, `rep_srpRNA`, `rep_rRNA`, `rep_RC`, `rep_scRNA`, `rep_RNA` |
| `annotate_vntrs()` | Chaisson BED | `rep_VNTR` |
| `annotate_strs()` | Chaisson BED | `rep_STR` |
| `annotate_cpg_islands()` | TSV | `cpg_islands` |
| `annotate_centromeres()` | TSV | `centromeres` |
| `annotate_asmb_gaps()` | TSV | `asmb_gaps` |
| `annotate_alt_haps()` | TSV | `alt_haps` |
| `annotate_gc_content()` | BigWig | `GC_content_left`, `GC_content_right` |

BND variants are split before annotation (`split_bnd()`), then results are aggregated (`aggregate_results()`).

Output: `{workdir}/{sample}_{ref}.SVs.ref.tsv`

### 3. Alignment Feature Collection (`AlignmentAnnotatorIllumina`)

```
dicast_lib/collect_illumina.py
```

Parallelized with joblib across `(chrom x sv_type)` combinations using `--threads`.

For each combination:
1. `calculate_coverage_baseline()` -- samples 1000 random 1000bp regions on the chromosome
2. `calculate_insertsize_baseline()` -- samples insert sizes from the chromosome
3. `calculate_mapping_quality_baseline()` -- samples mapping qualities
4. `annotate_coverage()` -- log2 fold-change coverage in each bin vs baseline
5. `annotate_read_based_features()` -- clip reads, split reads, discordant pairs per bin

Output: `{workdir}/{sample}_{ref}.SVs.aln.ill.{chrom}.{svtype}.tsv` (one file per chrom x svtype)

### 4. Feature Combination

`combine_feature_files()` in `dicast.py`:
- Merges raw + ref + all alignment TSVs on `id` column
- Output: `{workdir}/{sample}_{ref}.SVs.annot.tsv`

### 5. Scoring (`Dicast` / `DicastExome`)

```
dicast_lib/model.py, dicast_lib/model_exome.py
```

`score_variants()` in `dicast.py` iterates over SV types:
- **DEL, DUP, INS**: loads XGBoost model from `{models}/dicast_{SVTYPE}.json`, calls `impute_missing_values()` then `predict()`
- **INV**: no model -- uses `score_inversions()` heuristic (binary 0/1)
- **Exome mode**: uses `DicastExome` class with pickle models (`.pkl`), fewer alignment features (no coverage)

Output: `{workdir}/{sample}_{ref}.SVs.dicast.tsv`

### 6. VCF Annotation

`add_info_tag_to_vcf()` in `dicast.py`:
- Reads the dicast TSV, matches variants by caller
- Adds `DQ` INFO tag (dicast quality score) to each input VCF
- Output: `{original_vcf_path}/*.dicast.vcf`

## Cohort Mode Pipeline

See [cohort-mode.md](cohort-mode.md) for full details. Key differences from call mode:

1. Reads a CSV file instead of VCFs (`VariantPrep` in cohort mode via `read_csv()` + `reformat_csv()`)
2. Additional cohort-specific filtering: removes variants present in all samples, filters `qual < 0.1`
3. Creates a `Cohort` object to determine missing variants per sample
4. Runs reference + alignment + scoring per sample (same as call mode)
5. Post-scoring: updates cohort info, detects overlapping variants, writes `.regenotyped.csv`
6. Single-sample shortcut: if only 1 sample, copies CSV to `.regenotyped.csv` without processing

## Class Dependency Diagram

```
dicast.py
  |-- VariantPrep          (dicast_lib/prepare.py)
  |-- ReferenceAnnotator    (dicast_lib/collect_reference.py)
  |-- AlignmentAnnotatorIllumina  (dicast_lib/collect_illumina.py)
  |-- Dicast               (dicast_lib/model.py)
  |-- DicastExome           (dicast_lib/model_exome.py)
  |-- Cohort               (dicast_lib/cohort.py)

Eva                        (dicast_lib/evaluate.py)  -- used by eva.ipynb, not by dicast.py
```

## Data Flow Summary

```
VCFs (per caller)
  --> raw.tsv (unified variant table)
    --> ref.tsv (+ reference features)
    --> aln.ill.{chrom}.{svtype}.tsv (+ alignment features, parallelized)
      --> annot.tsv (merged features)
        --> dicast.tsv (+ quality scores)
          --> *.dicast.vcf (annotated input VCFs)
```

## Global Constants

- `chroms`: chr1-chr22 + chrX (defined at module level in `dicast.py`)
- `sv_types`: DEL, DUP, INS, INV
- Chromosome filtering: `--chrom` argument restricts processing to specified chromosomes
