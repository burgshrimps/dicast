# dicast — Structural-variant confidence scoring

**dicast** is a standalone command-line tool that scores the confidence of
structural-variant (SV) calls. Given any aligned BAM, the matching per-caller
VCFs (delly, manta, lumpy, cnvnator, gridss, …), and the genomic-context
annotations of a reference build, it emits a per-call quality score in `[0, 1]`
under the `qual_dicast` column (or as a `DQ` INFO tag in a re-emitted VCF).

dicast is **caller-agnostic and cohort-agnostic** — it does not depend on the
tGenVar paper. The trained models shipped here are paper artifacts, but you
can swap them for your own (any XGBoost JSON / pkl model produced by
`02_dicast/dicast_lib/model.py`) and dicast will use them transparently.

This is the third of five sub-archives that accompany the paper:

```
00_ground_truth/        SV ground-truth construction
01_variant_database/    SQLite DB + dicast model training
02_dicast/              ← you are here
03_cuban/               visualization library
04_figures/             paper figures
```

## Layout

```
.
├── README.md
├── dicast.py                   tool entry point — `python3 dicast.py {call,cohort}`
├── dicast.sh                   single-sample call example (edit paths, then run)
├── dicast_cohort.sh            cohort call example
├── dicast_benchmark_hg002.sh   reproduce the runtime/memory benchmark
├── environment.yml             conda environment
├── configs/                    feature-toggle YAMLs (config.yaml is the default)
├── dicast_lib/                 Python library (parsing, prepare, model, evaluate, …)
├── docs/                       per-area docs (architecture, CLI, features, models, cohort, dev)
└── models/                     ← trained XGBoost models for hg38 (see below)
```

## Shipped models

This archive ships the **five trained models referenced from the paper figure
scripts in `04_figures/`**. All are XGBoost JSON files trained on the tGenVar
9-patient cohort against the consensus ground truth from `00_ground_truth/`.

| File | SV type | Variant | Notes |
|---|---|---|---|
| `models/dicast_DEL.json`     | DEL | normal | trained on tGenVar consensus truth |
| `models/dicast_INS.json`     | INS | normal | trained on tGenVar consensus truth |
| `models/dicast_DUP.json`     | DUP | normal | trained on tGenVar consensus truth |
| `models/dicast_DEL_pop.json` | DEL | population | trained additionally with HGSVC PAV calls as positives — favours common-variant recall |
| `models/dicast_INS_pop.json` | INS | population | population-trained INS model (no DUP pop variant — HGSVC PAV does not provide DUP truth) |

Each model has a sibling `*_metadata.json` listing training cohort
(`SAMPLE_001..SAMPLE_009` placeholders), reference build, callers, features,
and class counts.

To use a different model, point `--models` at any directory containing one
JSON per SV type with the same naming pattern, or train your own:
`02_dicast/dicast_lib/model.py` is the training entry point and
`01_variant_database/scripts/03_train.py` is the driver used in this paper.

## Dependencies

- Python 3.10 with: pandas, numpy, scipy, scikit-learn, xgboost, networkx,
  matplotlib, seaborn, pyyaml, tqdm, joblib, pysam, pybedtools, pyBigWig,
  bioframe, vcfpy.
- bcftools, bgzip, tabix.

`environment.yml` ships a minimal conda spec; recreate with
`conda env create -f environment.yml`.

## Quickstart — score a single sample

1. Set `REF_DIR` to your local hg38 reference root (used to resolve the FASTA
   index `${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai`).
   Annotation TSVs ship with this archive under `annot/` and need no
   external setup.
2. Edit `dicast.sh` to point `BAM` and `VCF_DIR` at your inputs.
3. Run:
   ```
   bash dicast.sh
   ```
   Per-call quality scores land in `${WORKDIR}/<SAMPLE>.dicast.csv`. With
   `--vcf` (add the flag in `dicast.sh`), dicast also re-emits the input VCFs
   with a `DQ` INFO tag.

For cohort mode (joint-called CSV + pedigree + family filtering) see
`dicast_cohort.sh`. For the runtime / peak-memory benchmark see
`dicast_benchmark_hg002.sh`.

## Reference annotations

This archive ships the eight hg38 annotation files dicast consumes under
`annot/`. The shell scripts default to `${DICAST_DIR}/annot` so no extra
setup is needed; override `ANNOT` (in `dicast.sh`) or `ANNOT_DIR`
(in `dicast_benchmark_hg002.sh`) to point at your own copy.

- `hg38_repeatmasker.tsv`, `hg38_vntrs_chaisson.bed`, `hg38_strs_chaisson.bed`
- `hg38_cpg_islands.tsv`, `hg38_centromeres.tsv`
- `hg38_asmb_gaps.tsv`, `hg38_alt_haps.tsv`
- `hg38_gc_content.bw`

All eight are public — sources are documented in `annot/README.md` and
recipes for regenerating them live in `01_variant_database/`.
