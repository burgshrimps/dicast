#!/usr/bin/env bash
# Fetches the large hg38 annotation files that exceed GitHub's 100MB limit
# and are therefore attached to the "annotations-v1" GitHub Release instead
# of being tracked in this repo. Safe to re-run: files that already exist
# with a matching checksum are skipped.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANNOT_DIR="${SCRIPT_DIR}/annot"
CHECKSUM_FILE="${ANNOT_DIR}/checksums.md5"
RELEASE_URL="https://github.com/burgshrimps/dicast/releases/download/annotations-v1"

# filename -> URL, one pair per line. Edit here to point at a mirror.
FILE_NAMES=(
    "hg38_gc_content.bw"
    "hg38_repeatmasker.tsv"
)
FILE_URLS=(
    "${RELEASE_URL}/hg38_gc_content.bw"
    "${RELEASE_URL}/hg38_repeatmasker.tsv"
)

# Optional: published later by the maintainer, used by dicast's --pop flag.
# Left out of the required set above so a 404 here never fails the run.
OPTIONAL_FILE="pav_catalog_hg38.vcf.gz"
OPTIONAL_URL="${RELEASE_URL}/${OPTIONAL_FILE}"

echo "dicast annotation downloader"
echo "Target directory: ${ANNOT_DIR}"
echo "Required downloads: ~2.1 GB total (hg38_gc_content.bw ~1.6 GB, hg38_repeatmasker.tsv ~0.5 GB)"
echo

checksum_of() {
    local file="$1"
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "${file}" | awk '{print $1}'
    else
        md5 -r "${file}" | awk '{print $1}'
    fi
}

expected_checksum() {
    local name="$1"
    if [[ -f "${CHECKSUM_FILE}" ]]; then
        awk -v n="${name}" '$2 == n {print $1}' "${CHECKSUM_FILE}"
    fi
}

download() {
    local url="$1"
    local dest_part="$2"
    if command -v curl >/dev/null 2>&1; then
        curl -L -C - --fail -o "${dest_part}" "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget -c -O "${dest_part}" "${url}"
    else
        echo "  Error: neither curl nor wget is available." >&2
        return 1
    fi
}

fetch_file() {
    local name="$1"
    local url="$2"
    local required="$3"
    local dest="${ANNOT_DIR}/${name}"
    local dest_part="${dest}.part"
    local expected
    expected="$(expected_checksum "${name}")"

    echo "== ${name} =="

    if [[ -z "${expected}" ]]; then
        echo "  Warning: no checksum entry for ${name} in $(basename "${CHECKSUM_FILE}"); skipping verification."
    fi

    if [[ -f "${dest}" ]]; then
        if [[ -n "${expected}" ]]; then
            local actual
            actual="$(checksum_of "${dest}")"
            if [[ "${actual}" == "${expected}" ]]; then
                echo "  Already present and checksum OK, skipping."
                return 0
            fi
            echo "  Existing file failed checksum, re-downloading."
        else
            echo "  Already present, skipping (unverified)."
            return 0
        fi
    fi

    echo "  Downloading from ${url}"
    if ! download "${url}" "${dest_part}"; then
        rm -f "${dest_part}"
        if [[ "${required}" == "yes" ]]; then
            echo "  Download failed. Fetch it manually from: ${url}"
            echo "  and place it at: ${dest}"
            return 1
        else
            echo "  Not available yet (expected while the catalog is unpublished): ${url}"
            return 0
        fi
    fi

    if [[ -n "${expected}" ]]; then
        local actual
        actual="$(checksum_of "${dest_part}")"
        if [[ "${actual}" != "${expected}" ]]; then
            echo "  Checksum mismatch (expected ${expected}, got ${actual})."
            rm -f "${dest_part}"
            echo "  Fetch it manually from: ${url}"
            echo "  and place it at: ${dest}"
            return 1
        fi
    fi

    mv "${dest_part}" "${dest}"
    echo "  Done."
}

mkdir -p "${ANNOT_DIR}"

status=0
ok_count=0
total_required=${#FILE_NAMES[@]}

for i in "${!FILE_NAMES[@]}"; do
    if fetch_file "${FILE_NAMES[$i]}" "${FILE_URLS[$i]}" "yes"; then
        ok_count=$((ok_count + 1))
    else
        status=1
    fi
    echo
done

fetch_file "${OPTIONAL_FILE}" "${OPTIONAL_URL}" "no"
echo

echo "Summary: ${ok_count}/${total_required} required annotation files ready in ${ANNOT_DIR}"
if [[ ${status} -ne 0 ]]; then
    echo "One or more required downloads failed - see messages above for manual URLs."
    exit 1
fi
echo "All required annotation files are in place."
