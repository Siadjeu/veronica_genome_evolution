#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=annotation_table
#SBATCH --output=07_annotation/logs/annotation_table_%j.out
#SBATCH --error=07_annotation/logs/annotation_table_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
RUN_ROOT="${PROJECT_DIR}/07_annotation/braker4_runs"
QC_DIR="${PROJECT_DIR}/07_annotation/qc"

FEATURE_FILE="${QC_DIR}/annotation_feature_counts.tsv"

OUTPUT="${QC_DIR}/manuscript_annotation_statistics_final.tsv"
BUSCO_CHECK="${QC_DIR}/busco_summary_files.tsv"
RAW_SUMMARIES="${QC_DIR}/busco_raw_summaries.txt"

SPECIES=(VPAN VSCU VANA VARV VPER VSER VTRI VVER PMAJ)

mkdir -p "${QC_DIR}"

test -s "${FEATURE_FILE}" || {
    echo "ERROR: Missing annotation feature table:" >&2
    echo "${FEATURE_FILE}" >&2
    exit 1
}

declare -A SCIENTIFIC_NAME
declare -A EVIDENCE

SCIENTIFIC_NAME[VPAN]="Veronica panormitana"
SCIENTIFIC_NAME[VSCU]="Veronica scutellata"
SCIENTIFIC_NAME[VANA]="Veronica anagallis-aquatica"
SCIENTIFIC_NAME[VARV]="Veronica arvensis"
SCIENTIFIC_NAME[VPER]="Veronica persica"
SCIENTIFIC_NAME[VSER]="Veronica serpyllifolia"
SCIENTIFIC_NAME[VTRI]="Veronica triloba"
SCIENTIFIC_NAME[VVER]="Veronica verna"
SCIENTIFIC_NAME[PMAJ]="Plantago major"

EVIDENCE[VPAN]="RNA-seq + proteins"
EVIDENCE[VSCU]="Proteins"
EVIDENCE[VANA]="Proteins"
EVIDENCE[VARV]="Proteins"
EVIDENCE[VPER]="Proteins"
EVIDENCE[VSER]="Proteins"
EVIDENCE[VTRI]="Proteins"
EVIDENCE[VVER]="Proteins"
EVIDENCE[PMAJ]="Proteins"

# ============================================================
# Parse BUSCO summary
#
# Accepts standard BUSCO strings such as:
# C:98.7%[S:17.8%,D:80.9%],F:0.3%,M:1.0%,n:1614
# ============================================================

parse_busco_summary() {
    local SUMMARY_FILE="$1"

    python - "${SUMMARY_FILE}" <<'PY'
import re
import sys
from pathlib import Path

path = Path(sys.argv[1])

if not path.is_file() or path.stat().st_size == 0:
    print("NA\tNA\tNA\tNA\tNA\tNA")
    raise SystemExit(0)

text = path.read_text(errors="replace")

pattern = re.compile(
    r"C:\s*([0-9.]+)%"
    r"\s*\[\s*S:\s*([0-9.]+)%"
    r"\s*,\s*D:\s*([0-9.]+)%\s*\]"
    r"\s*,\s*F:\s*([0-9.]+)%"
    r"\s*,\s*M:\s*([0-9.]+)%"
    r"\s*,\s*n:\s*([0-9]+)"
)

matches = pattern.findall(text)

if not matches:
    print("NA\tNA\tNA\tNA\tNA\tNA")
else:
    complete, single, duplicated, fragmented, missing, total = matches[-1]
    print(
        "\t".join(
            [
                complete,
                single,
                duplicated,
                fragmented,
                missing,
                total,
            ]
        )
    )
PY
}

# ============================================================
# Output headers
# ============================================================

printf "species_code\tscientific_name\tevidence\tpredicted_genes\ttranscripts\tall_proteins\tcoding_sequences\trepresentative_proteins\ttranscripts_per_gene\tprotein_busco_complete_pct\tprotein_busco_single_copy_pct\tprotein_busco_duplicated_pct\tprotein_busco_fragmented_pct\tprotein_busco_missing_pct\tprotein_busco_n\tgenome_busco_complete_pct\tgenome_busco_single_copy_pct\tgenome_busco_duplicated_pct\tgenome_busco_fragmented_pct\tgenome_busco_missing_pct\tgenome_busco_n\n" \
    > "${OUTPUT}"

printf "species_code\tprotein_busco_summary\tgenome_busco_summary\tstatus\n" \
    > "${BUSCO_CHECK}"

: > "${RAW_SUMMARIES}"

# ============================================================
# Process each species
# ============================================================

for CODE in "${SPECIES[@]}"
do
    QC_RESULTS="${RUN_ROOT}/${CODE}/output/${CODE}/results/quality_control"

    PROTEIN_BUSCO="${QC_RESULTS}/busco_proteins_short_summary.txt"
    GENOME_BUSCO="${QC_RESULTS}/busco_genome_short_summary.txt"
    GENERAL_BUSCO="${QC_RESULTS}/busco_summary.txt"

    STATUS="PASS"

    if [[ ! -s "${PROTEIN_BUSCO}" ]]; then
        echo "WARNING: Missing protein BUSCO summary for ${CODE}" >&2
        STATUS="FAIL_MISSING_PROTEIN_BUSCO"
    fi

    if [[ ! -s "${GENOME_BUSCO}" ]]; then
        echo "WARNING: Missing genome BUSCO summary for ${CODE}" >&2

        if [[ "${STATUS}" == "PASS" ]]; then
            STATUS="FAIL_MISSING_GENOME_BUSCO"
        else
            STATUS="${STATUS};FAIL_MISSING_GENOME_BUSCO"
        fi
    fi

    {
        echo "===== ${CODE}: protein BUSCO ====="

        if [[ -s "${PROTEIN_BUSCO}" ]]; then
            cat "${PROTEIN_BUSCO}"
        else
            echo "MISSING: ${PROTEIN_BUSCO}"
        fi

        echo
        echo "===== ${CODE}: genome BUSCO ====="

        if [[ -s "${GENOME_BUSCO}" ]]; then
            cat "${GENOME_BUSCO}"
        else
            echo "MISSING: ${GENOME_BUSCO}"
        fi

        echo
        echo "===== ${CODE}: general BUSCO summary ====="

        if [[ -s "${GENERAL_BUSCO}" ]]; then
            cat "${GENERAL_BUSCO}"
        else
            echo "MISSING: ${GENERAL_BUSCO}"
        fi

        echo
    } >> "${RAW_SUMMARIES}"

    FEATURE_ROW=$(
        awk -F'\t' -v code="${CODE}" '
            NR > 1 && $1 == code {
                print
            }
        ' "${FEATURE_FILE}"
    )

    if [[ -z "${FEATURE_ROW}" ]]; then
        echo "ERROR: Missing annotation count row for ${CODE}" >&2
        exit 1
    fi

    GENES=$(awk -F'\t' '{print $2}' <<< "${FEATURE_ROW}")
    TRANSCRIPTS=$(awk -F'\t' '{print $4}' <<< "${FEATURE_ROW}")
    ALL_PROTEINS=$(awk -F'\t' '{print $5}' <<< "${FEATURE_ROW}")
    CDS=$(awk -F'\t' '{print $6}' <<< "${FEATURE_ROW}")
    REPRESENTATIVE=$(awk -F'\t' '{print $7}' <<< "${FEATURE_ROW}")
    AVG=$(awk -F'\t' '{print $8}' <<< "${FEATURE_ROW}")

    IFS=$'\t' read -r \
        PROT_C \
        PROT_S \
        PROT_D \
        PROT_F \
        PROT_M \
        PROT_N \
        <<< "$(parse_busco_summary "${PROTEIN_BUSCO}")"

    IFS=$'\t' read -r \
        GENOME_C \
        GENOME_S \
        GENOME_D \
        GENOME_F \
        GENOME_M \
        GENOME_N \
        <<< "$(parse_busco_summary "${GENOME_BUSCO}")"

    if [[ "${PROT_C}" == "NA" ]]; then
        echo "WARNING: Could not parse protein BUSCO for ${CODE}" >&2

        if [[ "${STATUS}" == "PASS" ]]; then
            STATUS="FAIL_PARSE_PROTEIN_BUSCO"
        else
            STATUS="${STATUS};FAIL_PARSE_PROTEIN_BUSCO"
        fi
    fi

    if [[ "${GENOME_C}" == "NA" ]]; then
        echo "WARNING: Could not parse genome BUSCO for ${CODE}" >&2

        if [[ "${STATUS}" == "PASS" ]]; then
            STATUS="FAIL_PARSE_GENOME_BUSCO"
        else
            STATUS="${STATUS};FAIL_PARSE_GENOME_BUSCO"
        fi
    fi

    printf "%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${PROTEIN_BUSCO}" \
        "${GENOME_BUSCO}" \
        "${STATUS}" \
        >> "${BUSCO_CHECK}"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${SCIENTIFIC_NAME[${CODE}]}" \
        "${EVIDENCE[${CODE}]}" \
        "${GENES}" \
        "${TRANSCRIPTS}" \
        "${ALL_PROTEINS}" \
        "${CDS}" \
        "${REPRESENTATIVE}" \
        "${AVG}" \
        "${PROT_C}" \
        "${PROT_S}" \
        "${PROT_D}" \
        "${PROT_F}" \
        "${PROT_M}" \
        "${PROT_N}" \
        "${GENOME_C}" \
        "${GENOME_S}" \
        "${GENOME_D}" \
        "${GENOME_F}" \
        "${GENOME_M}" \
        "${GENOME_N}" \
        >> "${OUTPUT}"
done

echo
echo "============================================================"
echo "BUSCO file validation"
echo "============================================================"

column -t -s $'\t' "${BUSCO_CHECK}"

echo
echo "============================================================"
echo "Final manuscript annotation statistics"
echo "============================================================"

column -t -s $'\t' "${OUTPUT}"

echo
echo "Final table:"
echo "${OUTPUT}"

echo
echo "BUSCO file check:"
echo "${BUSCO_CHECK}"

echo
echo "Raw BUSCO summaries:"
echo "${RAW_SUMMARIES}"
