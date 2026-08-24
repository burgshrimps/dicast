#!/usr/bin/env bash
# Reproduce the dicast runtime / peak-memory benchmark on HG002 (30x WGS Illumina).
# This is the script that produced the timing supplementary table for the paper.
# Adapt the path variables at the top to point at your local copy of the HG002
# inputs.

set -euo pipefail

PYTHON=${PYTHON:-python3}

# HG002 reference + inputs. Defaults are /path/to/... placeholders; override
# either inline or via the environment.
REF_FA=${REF_FA:-/path/to/reference/hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa}
REF_FAI=${REF_FAI:-${REF_FA}.fai}
# Annotations ship under ./annot in this archive.
DICAST_DIR=${DICAST_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
ANNOT_DIR=${ANNOT_DIR:-${DICAST_DIR}/annot}

BAM=${BAM:-/path/to/HG002/ill/bam/GRCh38.bwa_mem.pe.sorted.mdup.recal.bam}
VCF_DIR=${VCF_DIR:-/path/to/HG002/ill/vcf}

# Shipped models + scratch dir (DICAST_DIR set above).
MODELS_DIR=${MODELS_DIR:-${DICAST_DIR}/models}
WORKDIR=${WORKDIR:-${DICAST_DIR}/tmp_bench_hg002}

# Where to write the benchmark TSV + /usr/bin/time -v cross-check + hardware log.
BENCH_DIR=${BENCH_DIR:-${DICAST_DIR}/bench_output}
BENCH_TSV=${BENCH_DIR}/dicast_hg002_timing.tsv
TIME_LOG=${BENCH_DIR}/dicast_hg002_time_v.txt
HW_LOG=${BENCH_DIR}/dicast_hg002_hardware.txt

mkdir -p "${WORKDIR}" "${BENCH_DIR}"

# Record hardware for the supplementary-table footer.
{
  echo "# Host: $(hostname)"
  echo "# Date: $(date -Iseconds)"
  echo
  echo "## lscpu"
  lscpu
  echo
  echo "## meminfo (MemTotal)"
  grep -E '^MemTotal' /proc/meminfo
  echo
  echo "## uname"
  uname -a
  echo
  echo "## python"
  "${PYTHON}" --version
  echo "## python exe"
  echo "${PYTHON}"
} > "${HW_LOG}"

cd "${DICAST_DIR}"

# /usr/bin/time -v gives an independent cross-check on peak RSS and wall-time
# (compare against the in-process --benchmark TSV).
/usr/bin/time -v -o "${TIME_LOG}" \
    "${PYTHON}" dicast.py call \
        --cohort 1kg \
        --sample HG002 \
        --ref hg38 \
        --technology ill \
        --workdir "${WORKDIR}" \
        --fai "${REF_FAI}" \
        --annot-dir "${ANNOT_DIR}" \
        --bam "${BAM}" \
        --vcfs delly=${VCF_DIR}/delly/formatted_variants.vcf.gz \
               manta=${VCF_DIR}/manta/formatted_variants.vcf.gz \
               lumpy=${VCF_DIR}/lumpy/formatted_variants.vcf.gz \
               cnvnator=${VCF_DIR}/cnvnator/formatted_variants.vcf.gz \
               gridss=${VCF_DIR}/gridss/formatted_variants.vcf.gz \
        --models "${MODELS_DIR}" \
        --sv_types DEL DUP INS \
        --threads "${THREADS:-30}" \
        --benchmark "${BENCH_TSV}"

echo
echo "dicast benchmark TSV: ${BENCH_TSV}"
echo "/usr/bin/time -v log: ${TIME_LOG}"
echo "Hardware log:         ${HW_LOG}"
