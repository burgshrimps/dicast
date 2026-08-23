#!/usr/bin/env bash
# Single-sample dicast example.
# Edit the variables below to point at your own annotation, BAM, and VCF inputs,
# then run: bash dicast.sh
#
# Required environment: see environment.yml.

set -euo pipefail

# --- Inputs (edit these) ---
SAMPLE=${SAMPLE:-HG002}
REF=${REF:-hg38}
TECH=${TECH:-ill}
WORKDIR=${WORKDIR:-./workdir}

# Reference FASTA index (provide your own); annotations ship in ./annot.
DICAST_DIR=${DICAST_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
REF_DIR=${REF_DIR:-/path/to/reference/hg38}
FAI=${FAI:-${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai}
ANNOT=${ANNOT:-${DICAST_DIR}/annot}

# Aligned BAM + per-caller VCFs.
BAM=${BAM:-/path/to/aligned/${SAMPLE}.bam}
VCF_DIR=${VCF_DIR:-/path/to/vcfs/${SAMPLE}}

# Model directory shipped with this archive.
MODELS=${MODELS:-./models}

mkdir -p "${WORKDIR}"

python3 dicast.py call \
    --cohort tgenvar \
    --sample "${SAMPLE}" \
    --ref "${REF}" \
    --technology "${TECH}" \
    --workdir "${WORKDIR}" \
    --fai "${FAI}" \
    --repeats     "${ANNOT}/hg38_repeatmasker.tsv" \
    --cgis        "${ANNOT}/hg38_cpg_islands.tsv" \
    --centromeres "${ANNOT}/hg38_centromeres.tsv" \
    --gaps        "${ANNOT}/hg38_asmb_gaps.tsv" \
    --althaps     "${ANNOT}/hg38_alt_haps.tsv" \
    --vntrs       "${ANNOT}/hg38_vntrs_chaisson.bed" \
    --strs        "${ANNOT}/hg38_strs_chaisson.bed" \
    --gc          "${ANNOT}/hg38_gc_content.bw" \
    --bam "${BAM}" \
    --vcfs delly=${VCF_DIR}/delly/formatted_variants.vcf.gz \
           manta=${VCF_DIR}/manta/formatted_variants.vcf.gz \
           lumpy=${VCF_DIR}/lumpy/formatted_variants.vcf.gz \
           cnvnator=${VCF_DIR}/cnvnator/formatted_variants.vcf.gz \
           gridss=${VCF_DIR}/gridss/formatted_variants.vcf.gz \
    --models "${MODELS}" \
    --threads "${THREADS:-30}"
