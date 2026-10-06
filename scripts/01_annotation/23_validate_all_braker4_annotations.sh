#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=validate_braker4
#SBATCH --output=07_annotation/logs/validate_braker4_%j.out
#SBATCH --error=07_annotation/logs/validate_braker4_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
RUN_ROOT="${PROJECT_DIR}/07_annotation/braker4_runs"
QC_DIR="${PROJECT_DIR}/07_annotation/qc"

mkdir -p "${QC_DIR}"

SUMMARY="${QC_DIR}/annotation_summary.tsv"
MISSING="${QC_DIR}/annotation_missing_files.tsv"
FEATURES="${QC_DIR}/annotation_feature_counts.tsv"

SPECIES=(VPAN VSCU VANA VARV VPER VSER VTRI VVER PMAJ)

printf "species_code\tresults_dir\tgff3\tgtf\tall_proteins\tcoding_sequences\tlongest_proteins\tlongest_gtf\tstatus\n" \
    > "${SUMMARY}"

printf "species_code\tmissing_file\n" > "${MISSING}"

printf "species_code\tgenes\tmRNA\ttranscripts\tall_proteins\tcoding_sequences\tlongest_proteins\tavg_transcripts_per_gene\n" \
    > "${FEATURES}"

open_stream() {
    local FILE="$1"

    if [[ "${FILE}" == *.gz ]]; then
        gzip -cd "${FILE}"
    else
        cat "${FILE}"
    fi
}

find_result() {
    local DIR="$1"
    shift

    local NAME
    local FILE=""

    for NAME in "$@"
    do
        FILE=$(find "${DIR}" \
            -maxdepth 1 \
            -type f \
            -name "${NAME}" \
            -size +0c \
            -print -quit 2>/dev/null || true)

        if [[ -n "${FILE}" ]]; then
            printf "%s\n" "${FILE}"
            return 0
        fi
    done

    return 1
}

for CODE in "${SPECIES[@]}"
do
    RESULTS_DIR="${RUN_ROOT}/${CODE}/output/${CODE}/results"

    echo
    echo "============================================================"
    echo "Validating ${CODE}"
    echo "============================================================"

    GFF3=$(find_result "${RESULTS_DIR}" "braker.gff3.gz" "braker.gff3" || true)
    GTF=$(find_result "${RESULTS_DIR}" "braker.gtf.gz" "braker.gtf" || true)
    AA=$(find_result "${RESULTS_DIR}" "braker.aa.gz" "braker.aa" || true)
    CDS=$(find_result "${RESULTS_DIR}" "braker.codingseq.gz" "braker.codingseq" || true)
    LONGEST_AA=$(find_result "${RESULTS_DIR}" "braker.longest.aa.gz" "braker.longest.aa" || true)
    LONGEST_GTF=$(find_result "${RESULTS_DIR}" "braker.longest.gtf.gz" "braker.longest.gtf" || true)

    STATUS="PASS"

    declare -A FILES=(
        [braker_gff3]="${GFF3}"
        [braker_gtf]="${GTF}"
        [braker_aa]="${AA}"
        [braker_codingseq]="${CDS}"
        [braker_longest_aa]="${LONGEST_AA}"
        [braker_longest_gtf]="${LONGEST_GTF}"
    )

    for LABEL in \
        braker_gff3 \
        braker_gtf \
        braker_aa \
        braker_codingseq \
        braker_longest_aa \
        braker_longest_gtf
    do
        FILE="${FILES[${LABEL}]}"

        if [[ -z "${FILE}" || ! -s "${FILE}" ]]; then
            printf "%s\t%s\n" "${CODE}" "${LABEL}" >> "${MISSING}"
            STATUS="FAIL"
        elif [[ "${FILE}" == *.gz ]]; then
            if ! gzip -t "${FILE}"; then
                printf "%s\t%s_corrupt\n" "${CODE}" "${LABEL}" >> "${MISSING}"
                STATUS="FAIL"
            fi
        fi
    done

    GENES=0
    MRNA=0
    TRANSCRIPTS=0
    ALL_PROTEINS=0
    CDS_COUNT=0
    LONGEST_COUNT=0
    AVG="NA"

    if [[ -n "${GFF3}" && -s "${GFF3}" ]]; then
        GENES=$(
            open_stream "${GFF3}" |
            awk -F'\t' '
                $0 !~ /^#/ && NF >= 3 && $3 == "gene" {
                    n++
                }
                END {
                    print n+0
                }
            '
        )

        MRNA=$(
            open_stream "${GFF3}" |
            awk -F'\t' '
                $0 !~ /^#/ && NF >= 3 && $3 == "mRNA" {
                    n++
                }
                END {
                    print n+0
                }
            '
        )

        TRANSCRIPTS=$(
            open_stream "${GFF3}" |
            awk -F'\t' '
                $0 !~ /^#/ && NF >= 3 &&
                ($3 == "transcript" || $3 == "mRNA") {
                    n++
                }
                END {
                    print n+0
                }
            '
        )
    fi

    if [[ -n "${AA}" && -s "${AA}" ]]; then
        ALL_PROTEINS=$(open_stream "${AA}" | grep -c '^>' || true)
    fi

    if [[ -n "${CDS}" && -s "${CDS}" ]]; then
        CDS_COUNT=$(open_stream "${CDS}" | grep -c '^>' || true)
    fi

    if [[ -n "${LONGEST_AA}" && -s "${LONGEST_AA}" ]]; then
        LONGEST_COUNT=$(open_stream "${LONGEST_AA}" | grep -c '^>' || true)
    fi

    if [[ "${GENES}" -gt 0 ]]; then
        AVG=$(
            awk -v transcripts="${ALL_PROTEINS}" -v genes="${GENES}" \
                'BEGIN {printf "%.3f", transcripts/genes}'
        )
    fi

    if [[ "${ALL_PROTEINS}" -ne "${CDS_COUNT}" ]]; then
        echo "WARNING: ${CODE} protein and CDS counts differ." >&2
        STATUS="CHECK"
    fi

    if [[ "${LONGEST_COUNT}" -gt "${ALL_PROTEINS}" ]]; then
        echo "WARNING: ${CODE} longest-protein count exceeds all proteins." >&2
        STATUS="CHECK"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${RESULTS_DIR}" \
        "${GFF3:-NA}" \
        "${GTF:-NA}" \
        "${AA:-NA}" \
        "${CDS:-NA}" \
        "${LONGEST_AA:-NA}" \
        "${LONGEST_GTF:-NA}" \
        "${STATUS}" \
        >> "${SUMMARY}"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${GENES}" \
        "${MRNA}" \
        "${TRANSCRIPTS}" \
        "${ALL_PROTEINS}" \
        "${CDS_COUNT}" \
        "${LONGEST_COUNT}" \
        "${AVG}" \
        >> "${FEATURES}"

    echo "Genes:               ${GENES}"
    echo "mRNA features:       ${MRNA}"
    echo "Transcript features: ${TRANSCRIPTS}"
    echo "All proteins:        ${ALL_PROTEINS}"
    echo "Coding sequences:    ${CDS_COUNT}"
    echo "Longest proteins:    ${LONGEST_COUNT}"
    echo "Status:               ${STATUS}"

    unset FILES
done

echo
echo "============================================================"
echo "Annotation validation finished"
echo "============================================================"
echo "Summary:"
echo "${SUMMARY}"
echo
echo "Feature counts:"
echo "${FEATURES}"
echo
echo "Missing files:"
echo "${MISSING}"

column -t -s $'\t' "${FEATURES}"
