#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-11
#SBATCH --job-name=wgdi36B
#SBATCH --output=11_wgdi/logs/step36B_%A_%a.out
#SBATCH --error=11_wgdi/logs/step36B_%A_%a.err

set -euo pipefail

###############################################################################
# PROJECT CONFIGURATION
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"
PEP_DIR="${WGDI_ROOT}/01_inputs/pep"
MANIFEST="${WGDI_ROOT}/00_admin/step36B_homology_comparisons.tsv"

DB_DIR="${WGDI_ROOT}/04_homology/diamond_databases"
RESULT_DIR="${WGDI_ROOT}/04_homology/results"
RAW_DIR="${WGDI_ROOT}/04_homology/raw"
QC_DIR="${WGDI_ROOT}/02_qc/homology"
LOG_DIR="${WGDI_ROOT}/logs"

mkdir -p \
    "${DB_DIR}" \
    "${RESULT_DIR}" \
    "${RAW_DIR}" \
    "${QC_DIR}" \
    "${LOG_DIR}"

###############################################################################
# ENVIRONMENT
###############################################################################

cd "${PROJECT_ROOT}"

module purge
module load hpc-env/13.1 2>/dev/null || true

CONDA_EXE_PATH=""

if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
    CONDA_EXE_PATH="${CONDA_EXE}"
elif command -v conda >/dev/null 2>&1; then
    CONDA_EXE_PATH="$(command -v conda)"
elif [[ -x "${HOME}/miniforge3/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/miniforge3/bin/conda"
elif [[ -x "${HOME}/mambaforge/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/mambaforge/bin/conda"
elif [[ -x "${HOME}/miniconda3/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/miniconda3/bin/conda"
else
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate wgdi_env

for PROGRAM in diamond python
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: Required program unavailable: ${PROGRAM}" >&2
        exit 1
    fi
done

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Comparison manifest is missing:" >&2
    echo "${MANIFEST}" >&2
    exit 1
fi

###############################################################################
# READ ARRAY TASK
###############################################################################

TASK_LINE="$(
    awk -F $'\t' \
        -v task="${SLURM_ARRAY_TASK_ID}" \
        'NR > 1 && $1 == task {
            print
            exit
        }' \
        "${MANIFEST}"
)"

if [[ -z "${TASK_LINE}" ]]; then
    echo "ERROR: No comparison for array task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

IFS=$'\t' read -r \
    TASK_ID \
    COMPARISON_TYPE \
    SPECIES1 \
    SPECIES2 \
    <<< "${TASK_LINE}"

RUN_ID="${SPECIES1}_${SPECIES2}"

PEP1="${PEP_DIR}/${SPECIES1}.wgdi.pep.fa"
PEP2="${PEP_DIR}/${SPECIES2}.wgdi.pep.fa"

for FILE in "${PEP1}" "${PEP2}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required protein FASTA missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# PARAMETERS
###############################################################################

EVALUE="1e-5"
MAX_TARGET_SEQS="100"
SENSITIVE_MODE="--very-sensitive"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Standard BLAST tabular columns required by WGDI:
# qseqid sseqid pident length mismatch gapopen
# qstart qend sstart send evalue bitscore

###############################################################################
# OUTPUT PATHS
###############################################################################

RUN_RAW_DIR="${RAW_DIR}/${RUN_ID}"
RUN_QC_DIR="${QC_DIR}/${RUN_ID}"

mkdir -p \
    "${RUN_RAW_DIR}" \
    "${RUN_QC_DIR}"

FINAL_BLAST="${RESULT_DIR}/${RUN_ID}.blast.tsv"
FINAL_BLAST_GZ="${FINAL_BLAST}.gz"

TASK_SUMMARY="${RUN_QC_DIR}/${RUN_ID}.homology_qc.tsv"
COMMAND_LOG="${RUN_RAW_DIR}/${RUN_ID}.commands.txt"

rm -f \
    "${FINAL_BLAST}" \
    "${FINAL_BLAST_GZ}" \
    "${TASK_SUMMARY}" \
    "${COMMAND_LOG}" \
    "${RUN_RAW_DIR}"/*.tmp.tsv \
    "${RUN_RAW_DIR}"/*.raw.tsv \
    "${RUN_RAW_DIR}"/*.sorted.tsv

###############################################################################
# HELPERS
###############################################################################

count_fasta_ids()
{
    grep -c '^>' "$1"
}

build_database()
{
    local SPECIES="$1"
    local PROTEIN_FASTA="$2"
    local DATABASE="${DB_DIR}/${SPECIES}"

    if [[ ! -s "${DATABASE}.dmnd" ]]; then
        local TEMP_DATABASE="${DB_DIR}/${SPECIES}.tmp.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

        rm -f "${TEMP_DATABASE}.dmnd"

        diamond makedb \
            --in "${PROTEIN_FASTA}" \
            --db "${TEMP_DATABASE}" \
            --threads "${THREADS}"

        if [[ ! -s "${TEMP_DATABASE}.dmnd" ]]; then
            echo "ERROR: DIAMOND database was not created for ${SPECIES}." >&2
            exit 1
        fi

        # A second array task might create the same database concurrently.
        # Keep an already completed database if one appeared meanwhile.
        if [[ -s "${DATABASE}.dmnd" ]]; then
            rm -f "${TEMP_DATABASE}.dmnd"
        else
            mv \
                "${TEMP_DATABASE}.dmnd" \
                "${DATABASE}.dmnd"
        fi
    fi

    if [[ ! -s "${DATABASE}.dmnd" ]]; then
        echo "ERROR: DIAMOND database unavailable for ${SPECIES}." >&2
        exit 1
    fi
}

run_diamond()
{
    local QUERY_FASTA="$1"
    local DATABASE_SPECIES="$2"
    local OUTPUT_FILE="$3"

    diamond blastp \
        --query "${QUERY_FASTA}" \
        --db "${DB_DIR}/${DATABASE_SPECIES}" \
        --out "${OUTPUT_FILE}" \
        --outfmt 6 \
            qseqid \
            sseqid \
            pident \
            length \
            mismatch \
            gapopen \
            qstart \
            qend \
            sstart \
            send \
            evalue \
            bitscore \
        --evalue "${EVALUE}" \
        --max-target-seqs "${MAX_TARGET_SEQS}" \
        ${SENSITIVE_MODE} \
        --threads "${THREADS}" \
        --tmpdir "${RUN_RAW_DIR}"
}

###############################################################################
# RUN HOMOLOGY SEARCH
###############################################################################

echo "============================================================"
echo "Step 36B: protein homology search"
echo "============================================================"
echo "Task: ${TASK_ID}"
echo "Type: ${COMPARISON_TYPE}"
echo "Species 1: ${SPECIES1}"
echo "Species 2: ${SPECIES2}"
echo "Run ID: ${RUN_ID}"
echo "Start: $(date --iso-8601=seconds)"
echo "DIAMOND: $(diamond version)"
echo "Threads: ${THREADS}"

N_PROTEINS1="$(count_fasta_ids "${PEP1}")"
N_PROTEINS2="$(count_fasta_ids "${PEP2}")"

if [[ "${N_PROTEINS1}" -le 0 || "${N_PROTEINS2}" -le 0 ]]; then
    echo "ERROR: Empty protein input." >&2
    exit 1
fi

build_database "${SPECIES1}" "${PEP1}"

if [[ "${SPECIES2}" != "${SPECIES1}" ]]; then
    build_database "${SPECIES2}" "${PEP2}"
fi

if [[ "${COMPARISON_TYPE}" == "self" ]]; then
    RAW_SELF="${RUN_RAW_DIR}/${RUN_ID}.raw.tsv"

    cat > "${COMMAND_LOG}" <<EOF2
diamond blastp --query ${PEP1} --db ${DB_DIR}/${SPECIES1} --outfmt 6 --evalue ${EVALUE} --max-target-seqs ${MAX_TARGET_SEQS} --very-sensitive
EOF2

    run_diamond \
        "${PEP1}" \
        "${SPECIES1}" \
        "${RAW_SELF}"

    if [[ ! -s "${RAW_SELF}" ]]; then
        echo "ERROR: Empty self-comparison result for ${RUN_ID}." >&2
        exit 1
    fi

    # Remove exact duplicate rows only.
    LC_ALL=C sort -u \
        "${RAW_SELF}" \
        > "${FINAL_BLAST}"

elif [[ "${COMPARISON_TYPE}" == "pairwise" ]]; then
    RAW_FORWARD="${RUN_RAW_DIR}/${SPECIES1}_to_${SPECIES2}.raw.tsv"
    RAW_REVERSE="${RUN_RAW_DIR}/${SPECIES2}_to_${SPECIES1}.raw.tsv"
    REVERSED_NORMALIZED="${RUN_RAW_DIR}/${SPECIES2}_to_${SPECIES1}.normalized.tmp.tsv"

    cat > "${COMMAND_LOG}" <<EOF2
diamond blastp --query ${PEP1} --db ${DB_DIR}/${SPECIES2} --outfmt 6 --evalue ${EVALUE} --max-target-seqs ${MAX_TARGET_SEQS} --very-sensitive
diamond blastp --query ${PEP2} --db ${DB_DIR}/${SPECIES1} --outfmt 6 --evalue ${EVALUE} --max-target-seqs ${MAX_TARGET_SEQS} --very-sensitive
EOF2

    run_diamond \
        "${PEP1}" \
        "${SPECIES2}" \
        "${RAW_FORWARD}"

    run_diamond \
        "${PEP2}" \
        "${SPECIES1}" \
        "${RAW_REVERSE}"

    if [[ ! -s "${RAW_FORWARD}" || ! -s "${RAW_REVERSE}" ]]; then
        echo "ERROR: One reciprocal result is empty for ${RUN_ID}." >&2
        exit 1
    fi

    # Normalize reverse hits so that species1 IDs are always column 1
    # and species2 IDs are always column 2.
    #
    # BLAST columns:
    # 1 qseqid, 2 sseqid, 3 pident, 4 length, 5 mismatch,
    # 6 gapopen, 7 qstart, 8 qend, 9 sstart, 10 send,
    # 11 evalue, 12 bitscore
    awk -F $'\t' '
        BEGIN {
            OFS="\t"
        }
        {
            print \
                $2, $1, $3, $4, $5, $6, \
                $9, $10, $7, $8, $11, $12
        }
    ' "${RAW_REVERSE}" \
        > "${REVERSED_NORMALIZED}"

    # Reciprocal searches can produce the same pair twice with slightly
    # different scores. Retain one row per gene pair, prioritizing:
    # lowest E-value, then highest bitscore, then highest identity.
    cat \
        "${RAW_FORWARD}" \
        "${REVERSED_NORMALIZED}" |
    LC_ALL=C sort \
        -t $'\t' \
        -k1,1 \
        -k2,2 \
        -k11,11g \
        -k12,12gr \
        -k3,3gr |
    awk -F $'\t' '
        BEGIN {
            OFS="\t"
        }
        {
            pair=$1 SUBSEP $2

            if (!(pair in seen)) {
                print
                seen[pair]=1
            }
        }
    ' > "${FINAL_BLAST}"

else
    echo "ERROR: Unsupported comparison type: ${COMPARISON_TYPE}" >&2
    exit 1
fi

if [[ ! -s "${FINAL_BLAST}" ]]; then
    echo "ERROR: Final BLAST-format table is empty." >&2
    exit 1
fi

###############################################################################
# STRICT OUTPUT VALIDATION
###############################################################################

python - \
    "${FINAL_BLAST}" \
    "${PEP1}" \
    "${PEP2}" \
    "${SPECIES1}" \
    "${SPECIES2}" \
    "${COMPARISON_TYPE}" \
    "${TASK_SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import math
import sys
from collections import Counter
from pathlib import Path

blast_path = Path(sys.argv[1])
pep1_path = Path(sys.argv[2])
pep2_path = Path(sys.argv[3])
species1 = sys.argv[4]
species2 = sys.argv[5]
comparison_type = sys.argv[6]
summary_path = Path(sys.argv[7])


def fasta_ids(path: Path):
    result = set()

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                identifier = line[1:].split()[0]

                if identifier in result:
                    raise SystemExit(
                        f"ERROR: Duplicate FASTA ID in {path}: "
                        f"{identifier}"
                    )

                result.add(identifier)

    if not result:
        raise SystemExit(f"ERROR: No FASTA IDs in {path}")

    return result


ids1 = fasta_ids(pep1_path)
ids2 = fasta_ids(pep2_path)

row_count = 0
malformed_rows = 0
invalid_id_rows = 0
self_identity_rows = 0
unique_queries = set()
unique_subjects = set()
pair_counter = Counter()

identity_sum = 0.0
bitscore_sum = 0.0
evalue_zero_count = 0

with blast_path.open("r", encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, start=1):
        fields = line.rstrip("\n").split("\t")

        if len(fields) != 12:
            malformed_rows += 1
            continue

        query_id = fields[0]
        subject_id = fields[1]

        try:
            pident = float(fields[2])
            alignment_length = int(fields[3])
            mismatch = int(fields[4])
            gapopen = int(fields[5])
            qstart = int(fields[6])
            qend = int(fields[7])
            sstart = int(fields[8])
            send = int(fields[9])
            evalue = float(fields[10])
            bitscore = float(fields[11])
        except ValueError:
            malformed_rows += 1
            continue

        if comparison_type == "self":
            valid_ids = (
                query_id in ids1
                and subject_id in ids1
            )
        else:
            valid_ids = (
                query_id in ids1
                and subject_id in ids2
            )

        if not valid_ids:
            invalid_id_rows += 1

        if query_id == subject_id:
            self_identity_rows += 1

        if not (0.0 <= pident <= 100.0):
            malformed_rows += 1

        if alignment_length <= 0:
            malformed_rows += 1

        if (
            mismatch < 0
            or gapopen < 0
            or qstart <= 0
            or qend <= 0
            or sstart <= 0
            or send <= 0
            or evalue < 0
            or bitscore < 0
        ):
            malformed_rows += 1

        row_count += 1
        unique_queries.add(query_id)
        unique_subjects.add(subject_id)
        pair_counter[(query_id, subject_id)] += 1

        identity_sum += pident
        bitscore_sum += bitscore

        if evalue == 0:
            evalue_zero_count += 1

duplicate_pairs = sum(
    count - 1
    for count in pair_counter.values()
    if count > 1
)

if row_count == 0:
    raise SystemExit("ERROR: No valid homology rows.")

if malformed_rows != 0:
    raise SystemExit(
        f"ERROR: {malformed_rows} malformed rows detected."
    )

if invalid_id_rows != 0:
    raise SystemExit(
        f"ERROR: {invalid_id_rows} rows contain invalid IDs."
    )

if duplicate_pairs != 0:
    raise SystemExit(
        f"ERROR: {duplicate_pairs} duplicate gene pairs remain."
    )

query_coverage = len(unique_queries) / len(ids1)

if comparison_type == "self":
    subject_denominator = len(ids1)
else:
    subject_denominator = len(ids2)

subject_coverage = (
    len(unique_subjects)
    / subject_denominator
)

status = "PASS"

row = {
    "comparison": f"{species1}_{species2}",
    "comparison_type": comparison_type,
    "species1": species1,
    "species2": species2,
    "species1_protein_count": len(ids1),
    "species2_protein_count": len(ids2),
    "homology_rows": row_count,
    "unique_gene_pairs": len(pair_counter),
    "unique_queries": len(unique_queries),
    "unique_subjects": len(unique_subjects),
    "query_gene_coverage": f"{query_coverage:.8f}",
    "subject_gene_coverage": f"{subject_coverage:.8f}",
    "mean_percent_identity": f"{identity_sum / row_count:.8f}",
    "mean_bitscore": f"{bitscore_sum / row_count:.8f}",
    "zero_evalue_rows": evalue_zero_count,
    "self_identity_rows": self_identity_rows,
    "malformed_rows": malformed_rows,
    "invalid_id_rows": invalid_id_rows,
    "duplicate_gene_pairs": duplicate_pairs,
    "status": status,
}

with summary_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(row),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerow(row)

print(
    f"{species1}_{species2}: "
    f"rows={row_count:,}; "
    f"query coverage={query_coverage:.4f}; "
    f"subject coverage={subject_coverage:.4f}; "
    f"status=PASS"
)
PY

###############################################################################
# COMPRESS FINAL RESULT
###############################################################################

gzip -f "${FINAL_BLAST}"

if [[ ! -s "${FINAL_BLAST_GZ}" ]]; then
    echo "ERROR: Compressed result was not created." >&2
    exit 1
fi

if [[ ! -s "${TASK_SUMMARY}" ]]; then
    echo "ERROR: Task QC summary was not created." >&2
    exit 1
fi

echo
echo "============================================================"
echo "Step 36B task summary"
echo "============================================================"

column -t -s $'\t' \
    "${TASK_SUMMARY}"

echo
echo "Output:"
echo "${FINAL_BLAST_GZ}"
echo
echo "Completed: $(date --iso-8601=seconds)"
echo "Step 36B task ${TASK_ID}: PASS"
