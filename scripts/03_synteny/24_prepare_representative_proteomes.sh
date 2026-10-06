#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=prepare_proteomes
#SBATCH --output=07_annotation/logs/prepare_proteomes_%j.out
#SBATCH --error=07_annotation/logs/prepare_proteomes_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
RUN_ROOT="${PROJECT_DIR}/07_annotation/braker4_runs"
OUT_DIR="${PROJECT_DIR}/08_comparative_inputs/proteomes"
MANIFEST="${PROJECT_DIR}/08_comparative_inputs/proteome_manifest.tsv"

SPECIES=(VPAN VSCU VANA VARV VPER VSER VTRI VVER PMAJ)

mkdir -p "${OUT_DIR}"
mkdir -p "$(dirname "${MANIFEST}")"

printf "species_code\tsource_file\tstandardised_file\tprotein_count\tduplicate_headers\tstatus\n" \
    > "${MANIFEST}"

for CODE in "${SPECIES[@]}"
do
    RESULTS_DIR="${RUN_ROOT}/${CODE}/output/${CODE}/results"

    SOURCE=$(find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.longest.aa.gz" -o -name "braker.longest.aa" \) \
        -size +0c \
        -print -quit 2>/dev/null || true)

    TARGET="${OUT_DIR}/${CODE}.longest.aa.fa"

    if [[ -z "${SOURCE}" ]]; then
        printf "%s\tNA\t%s\t0\tNA\tFAIL_MISSING_SOURCE\n" \
            "${CODE}" "${TARGET}" >> "${MANIFEST}"
        continue
    fi

    if [[ "${SOURCE}" == *.gz ]]; then
        gzip -cd "${SOURCE}" > "${TARGET}"
    else
        cp "${SOURCE}" "${TARGET}"
    fi

    if [[ ! -s "${TARGET}" ]]; then
        printf "%s\t%s\t%s\t0\tNA\tFAIL_EMPTY_TARGET\n" \
            "${CODE}" "${SOURCE}" "${TARGET}" >> "${MANIFEST}"
        continue
    fi

    COUNT=$(grep -c '^>' "${TARGET}")

    DUPLICATES=$(
        grep '^>' "${TARGET}" |
        sed 's/^>//' |
        awk '{print $1}' |
        sort |
        uniq -d |
        wc -l
    )

    if [[ "${DUPLICATES}" -eq 0 ]]; then
        STATUS="PASS"
    else
        STATUS="FAIL_DUPLICATE_HEADERS"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${SOURCE}" \
        "${TARGET}" \
        "${COUNT}" \
        "${DUPLICATES}" \
        "${STATUS}" \
        >> "${MANIFEST}"
done

column -t -s $'\t' "${MANIFEST}"
