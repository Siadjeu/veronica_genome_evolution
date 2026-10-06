#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --job-name=prepare_braker4
#SBATCH --output=07_annotation/logs/prepare_braker4_%j.out
#SBATCH --error=07_annotation/logs/prepare_braker4_%j.err

set -euo pipefail

# ============================================================
# Fixed project directory
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
cd "${PROJECT_DIR}"

GENOME_SOURCE_DIR="${PROJECT_DIR}/02_genomes/public"
REPORT_DIR="${PROJECT_DIR}/00_metadata/sequence_reports"

OUTDIR="${PROJECT_DIR}/07_annotation/inputs/genomes"
MAPDIR="${PROJECT_DIR}/07_annotation/inputs/header_maps"
IDDIR="${PROJECT_DIR}/07_annotation/inputs/primary_assembly_ids"
QCDIR="${PROJECT_DIR}/07_annotation/qc"
LOGDIR="${PROJECT_DIR}/07_annotation/logs"
TMPDIR="${PROJECT_DIR}/07_annotation/tmp/genome_preparation"

QC_TABLE="${QCDIR}/braker4_genome_input_qc.tsv"

mkdir -p \
    "${OUTDIR}" \
    "${MAPDIR}" \
    "${IDDIR}" \
    "${QCDIR}" \
    "${LOGDIR}" \
    "${TMPDIR}"

# ============================================================
# Activate tools
# ============================================================

CONDA_SH="${CONDA_SH:-${HOME}/miniforge3/etc/profile.d/conda.sh}"
source "${CONDA_SH}"
conda activate genome_download

for TOOL in jq seqkit gzip awk sha256sum
do
    if ! command -v "${TOOL}" >/dev/null 2>&1
    then
        echo "ERROR: Required tool unavailable: ${TOOL}" >&2
        exit 1
    fi
done

# ============================================================
# Species included in the primary annotation set
# ============================================================

CODES=(
    VPAN
    VSCU
    VANA
    VARV
    VPER
    VSER
    VTRI
    VVER
    PMAJ
)

printf "species_code\tsource_genome\tsource_sequences\tselected_nuclear_sequences\toutput_genome\toutput_sequences\toutput_bp\tsha256\tstatus\n" \
    > "${QC_TABLE}"

# ============================================================
# Process genomes
# ============================================================

for CODE in "${CODES[@]}"
do
    SOURCE="${GENOME_SOURCE_DIR}/${CODE}.genome.fa.gz"
    REPORT="${REPORT_DIR}/${CODE}.sequence_report.jsonl"

    IDFILE="${IDDIR}/${CODE}.primary_nuclear_ids.txt"
    EXTRACTED="${TMPDIR}/${CODE}.selected.fa"
    OUTPUT="${OUTDIR}/${CODE}.nuclear.fa"
    HEADER_MAP="${MAPDIR}/${CODE}.header_map.tsv"

    echo
    echo "============================================================"
    echo "Preparing ${CODE}"
    echo "Source: ${SOURCE}"
    echo "Output: ${OUTPUT}"
    echo "============================================================"

    if [[ ! -s "${SOURCE}" ]]
    then
        echo "ERROR: Missing source genome for ${CODE}: ${SOURCE}" >&2
        exit 1
    fi

    gzip -t "${SOURCE}"

    SOURCE_COUNT=$(zgrep -c '^>' "${SOURCE}" || true)

    if [[ "${SOURCE_COUNT}" -eq 0 ]]
    then
        echo "ERROR: No FASTA records in ${SOURCE}" >&2
        exit 1
    fi

    rm -f \
        "${IDFILE}" \
        "${EXTRACTED}" \
        "${OUTPUT}" \
        "${HEADER_MAP}"

    # --------------------------------------------------------
    # VPAN is the locally assembled nuclear genome
    # --------------------------------------------------------

    if [[ "${CODE}" == "VPAN" ]]
    then
        gzip -cd "${SOURCE}" > "${EXTRACTED}"
        SELECTED_COUNT="${SOURCE_COUNT}"

    # --------------------------------------------------------
    # NCBI genomes: retain only Primary Assembly records
    #
    # This keeps:
    #   assembled-molecule
    #   unlocalized-scaffold
    #   unplaced-scaffold
    #
    # It excludes non-nuclear assembly units.
    # --------------------------------------------------------

    else
        if [[ ! -s "${REPORT}" ]]
        then
            echo "ERROR: Missing sequence report for ${CODE}: ${REPORT}" >&2
            exit 1
        fi

        if ! jq -e . "${REPORT}" >/dev/null 2>&1
        then
            echo "ERROR: Invalid JSONL sequence report: ${REPORT}" >&2
            exit 1
        fi

        jq -r '
            select(.assemblyUnit == "Primary Assembly")
            |
            (
                .genbankAccession
                // .refseqAccession
                // .sequenceName
                // empty
            )
        ' "${REPORT}" |
        sed '/^[[:space:]]*$/d' |
        sort -u > "${IDFILE}"

        SELECTED_COUNT=$(wc -l < "${IDFILE}")

        if [[ "${SELECTED_COUNT}" -eq 0 ]]
        then
            echo "ERROR: No Primary Assembly sequence IDs found for ${CODE}" >&2
            exit 1
        fi

        seqkit grep \
            --pattern-file "${IDFILE}" \
            "${SOURCE}" \
            > "${EXTRACTED}"
    fi

    if [[ ! -s "${EXTRACTED}" ]]
    then
        echo "ERROR: No sequences extracted for ${CODE}" >&2
        exit 1
    fi

    EXTRACTED_COUNT=$(grep -c '^>' "${EXTRACTED}" || true)

    if [[ "${EXTRACTED_COUNT}" -ne "${SELECTED_COUNT}" ]]
    then
        echo "ERROR: Sequence extraction mismatch for ${CODE}" >&2
        echo "Expected:  ${SELECTED_COUNT}" >&2
        echo "Extracted: ${EXTRACTED_COUNT}" >&2

        if [[ -s "${IDFILE}" ]]
        then
            echo "First requested identifiers:" >&2
            head -n 20 "${IDFILE}" >&2
        fi

        echo "First source FASTA headers:" >&2
        zgrep '^>' "${SOURCE}" | head -n 20 >&2

        exit 1
    fi

    # --------------------------------------------------------
    # Rename all sequences deterministically
    # --------------------------------------------------------

    awk \
        -v code="${CODE}" \
        -v mapfile="${HEADER_MAP}" '
        BEGIN {
            OFS="\t"
        }

        /^>/ {
            original=$0
            sub(/^>/, "", original)

            split(original, parts, /[[:space:]]+/)
            original_id=parts[1]

            record_number++
            new_id=sprintf("%s_seq%06d", code, record_number)

            print original_id, new_id >> mapfile
            print ">" new_id
            next
        }

        {
            gsub(/[[:space:]]/, "")
            print toupper($0)
        }
    ' "${EXTRACTED}" > "${OUTPUT}"

    rm -f "${EXTRACTED}"

    # --------------------------------------------------------
    # Validate prepared genome
    # --------------------------------------------------------

    if [[ ! -s "${OUTPUT}" ]]
    then
        echo "ERROR: Empty prepared genome for ${CODE}: ${OUTPUT}" >&2
        exit 1
    fi

    OUTPUT_COUNT=$(grep -c '^>' "${OUTPUT}" || true)

    if [[ "${OUTPUT_COUNT}" -ne "${SELECTED_COUNT}" ]]
    then
        echo "ERROR: Header-renaming count mismatch for ${CODE}" >&2
        echo "Expected: ${SELECTED_COUNT}" >&2
        echo "Output:   ${OUTPUT_COUNT}" >&2
        exit 1
    fi

    if grep '^>' "${OUTPUT}" | grep -q '[[:space:]]'
    then
        echo "ERROR: Whitespace remains in FASTA headers for ${CODE}" >&2
        exit 1
    fi

    if grep -v '^>' "${OUTPUT}" | grep -q '[^ACGTRYSWKMBDHVNacgtryswkmbdhvn]'
    then
        echo "ERROR: Unexpected nucleotide characters in ${OUTPUT}" >&2
        exit 1
    fi

    OUTPUT_BP=$(
        seqkit stats -T "${OUTPUT}" |
        awk 'NR == 2 {print $5}'
    )

    HASH=$(sha256sum "${OUTPUT}" | awk '{print $1}')

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tREADY\n" \
        "${CODE}" \
        "${SOURCE}" \
        "${SOURCE_COUNT}" \
        "${SELECTED_COUNT}" \
        "${OUTPUT}" \
        "${OUTPUT_COUNT}" \
        "${OUTPUT_BP}" \
        "${HASH}" \
        >> "${QC_TABLE}"

    echo "Prepared successfully:"
    echo "  Source sequences:   ${SOURCE_COUNT}"
    echo "  Nuclear sequences:  ${OUTPUT_COUNT}"
    echo "  Nuclear bases:      ${OUTPUT_BP}"
    echo "  Output:             ${OUTPUT}"
done

# ============================================================
# Final global validation
# ============================================================

EXPECTED_GENOMES=${#CODES[@]}

READY_GENOMES=$(
    awk -F'\t' '
        NR > 1 && $9 == "READY" {
            count++
        }

        END {
            print count+0
        }
    ' "${QC_TABLE}"
)

FILES_PRESENT=0

for CODE in "${CODES[@]}"
do
    OUTPUT="${OUTDIR}/${CODE}.nuclear.fa"

    if [[ -s "${OUTPUT}" ]]
    then
        FILES_PRESENT=$((FILES_PRESENT + 1))
    else
        echo "ERROR: Expected prepared genome is missing: ${OUTPUT}" >&2
        exit 1
    fi
done

if [[ "${READY_GENOMES}" -ne "${EXPECTED_GENOMES}" ]]
then
    echo "ERROR: Expected ${EXPECTED_GENOMES} READY genomes, found ${READY_GENOMES}" >&2
    exit 1
fi

if [[ "${FILES_PRESENT}" -ne "${EXPECTED_GENOMES}" ]]
then
    echo "ERROR: Expected ${EXPECTED_GENOMES} genome files, found ${FILES_PRESENT}" >&2
    exit 1
fi

echo
echo "============================================================"
echo "BRAKER4 genome preparation completed successfully"
echo "============================================================"
echo "Prepared genomes: ${FILES_PRESENT}/${EXPECTED_GENOMES}"
echo
column -t -s $'\t' "${QC_TABLE}" || cat "${QC_TABLE}"
