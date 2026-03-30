# CLI Reference

Source: `dicast_lib/parsing.py`

```
python3 dicast.py {call,cohort} [options]
```

## call -- Single-Sample Feature Extraction + Prediction

### Required Arguments

| Argument | Description |
|---|---|
| `--sample` | Sample name |
| `--workdir` | Working and output directory |
| `--models` | Directory containing trained models |
| `--bam` | BAM file path |
| `--vcfs` | Space-separated `caller=path` pairs (e.g., `delly=/path/to.vcf manta=/path/to.vcf`) |

### Optional Arguments

| Argument | Default | Description |
|---|---|---|
| `--cohort` | None | Cohort name |
| `--chrom` | `all` | Chromosome(s), space-separated (e.g., `--chrom chr1 chr2`) |
| `--ref` | `hg38` | Reference genome name |
| `--technology` | `ill` | Sequencing technology name |
| `--threads` | `1` | Number of parallel jobs for alignment feature collection |
| `--exome` | False | Flag for exome sequencing data |
| `--exome_regions` | None | BED file with enrichment kit regions (**required** if `--exome` is set) |

### Reference Files

All optional, but needed for reference feature collection:

| Argument | Format | Description |
|---|---|---|
| `--fai` | FAI | Reference genome index (chromosome sizes) |
| `--repeats` | TSV | RepeatMasker annotations |
| `--cgis` | TSV | CpG island annotations |
| `--centromeres` | TSV | Centromere annotations |
| `--gaps` | TSV | Assembly gap annotations |
| `--althaps` | TSV | Alternative haplotype annotations |
| `--vntrs` | BED | VNTR regions (Chaisson) |
| `--strs` | BED | STR regions (Chaisson) |
| `--gc` | BigWig | GC content |

## cohort -- Multi-Sample Cohort Analysis

### Required Arguments

| Argument | Description |
|---|---|
| `--cohort` | Cohort name |
| `--workdir` | Working and output directory |
| `--fai` | FAI file of the reference genome |
| `--models` | Directory with trained models |

### Cohort-Specific Arguments

| Argument | Default | Description |
|---|---|---|
| `--csv` | None | CSV file with variants (see [cohort-mode.md](cohort-mode.md) for format) |
| `--bams` | None | Space-separated BAM file paths. Sample name is extracted from `path.split('/')[-2]` |
| `--ped` | None | PED file with family information (TSV: Name, Family, Mother, Father) |
| `--filter-fam` | False | Restrict regenotyping to family members only |

### Shared Optional Arguments

Same as call mode: `--chrom`, `--ref`, `--technology`, `--threads`, `--exome`, `--exome_regions`, and all reference files.

## Validation Rules

- `--exome` requires `--exome_regions` to be set (for both `call` and `cohort`)
- In cohort mode: number of BAM files must match number of unique samples in the CSV
- VCF paths in `--vcfs` are split on `=` (format: `caller=path`)

## Full Examples

See `dicast.sh` (call mode) and `dicast_cohort.sh` (cohort mode) for complete working examples.
