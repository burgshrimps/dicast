#!/usr/bin/env bash
# Multi-sample dicast example: score a trio, with cross-sample rescue.
#
# Each sample is scored against its own callers' variants PLUS variants that
# were only found by the OTHER supplied samples' callers (e.g. a deletion the
# child's callers found but the mother's callers missed still gets evaluated
# against the mother's BAM). Rescued rows are flagged in the output TSV via
# their 'caller' column, e.g. 'rescue:CHILD:delly'.
#
# Edit the variables below to point at your own annotation, BAM, and VCF
# inputs, then run: bash dicast_multi.sh
#
# Required environment: see environment.yml.

set -euo pipefail

# --- Inputs (edit these) ---
MOTHER=${MOTHER:-HG004}
FATHER=${FATHER:-HG003}
CHILD=${CHILD:-HG002}
REF=${REF:-hg38}
TECH=${TECH:-ill}
WORKDIR=${WORKDIR:-./workdir}

# Reference FASTA index (provide your own); annotations ship in ./annot and
# are picked up automatically from there (override with ANNOT_DIR if needed).
DICAST_DIR=${DICAST_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
REF_DIR=${REF_DIR:-/path/to/reference/hg38}
FAI=${FAI:-${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai}
ANNOT_DIR=${ANNOT_DIR:-${DICAST_DIR}/annot}

# Aligned BAM + per-sample, per-caller VCFs.
BAM_DIR=${BAM_DIR:-/path/to/aligned}
VCF_DIR=${VCF_DIR:-/path/to/vcfs}

# Model directory shipped with this archive.
MODELS=${MODELS:-./models}

mkdir -p "${WORKDIR}"

python3 dicast.py multi \
    --cohort trio \
    --ref "${REF}" \
    --technology "${TECH}" \
    --workdir "${WORKDIR}" \
    --fai "${FAI}" \
    --annot-dir "${ANNOT_DIR}" \
    --bams "${MOTHER}=${BAM_DIR}/${MOTHER}.bam" \
           "${FATHER}=${BAM_DIR}/${FATHER}.bam" \
           "${CHILD}=${BAM_DIR}/${CHILD}.bam" \
    --vcfs "${MOTHER}:delly=${VCF_DIR}/${MOTHER}/delly/formatted_variants.vcf.gz" \
           "${MOTHER}:manta=${VCF_DIR}/${MOTHER}/manta/formatted_variants.vcf.gz" \
           "${FATHER}:delly=${VCF_DIR}/${FATHER}/delly/formatted_variants.vcf.gz" \
           "${FATHER}:manta=${VCF_DIR}/${FATHER}/manta/formatted_variants.vcf.gz" \
           "${CHILD}:delly=${VCF_DIR}/${CHILD}/delly/formatted_variants.vcf.gz" \
           "${CHILD}:manta=${VCF_DIR}/${CHILD}/manta/formatted_variants.vcf.gz" \
    --models "${MODELS}" \
    --threads "${THREADS:-30}"
