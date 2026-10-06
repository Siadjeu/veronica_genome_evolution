#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-2%3
#SBATCH --job-name=wgdi36X1
#SBATCH --output=11_wgdi/logs/step36X1_%A_%a.out
#SBATCH --error=11_wgdi/logs/step36X1_%A_%a.err

set -euo pipefail

###############################################################################
# PROJECT CONFIGURATION
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

PEP_DIR="${WGDI_ROOT}/01_inputs/pep"

DB_DIR="${WGDI_ROOT}/04_homology/diamond_databases"
RESULT_DIR="${WGDI_ROOT}/04_homology/results"
RAW_DIR="${WGDI_ROOT}/04_homology/raw"
QC_DIR="${WGDI_ROOT}/02_qc/homology"
LOG_DIR="${WGDI_ROOT}/logs"

ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X1"

X0_MANIFEST="${ADMIN_DIR}/step36X_additional_pairwise_comparisons.tsv"

X1_MANIFEST="${ADMIN_DIR}/step36X1_homology_comparisons.tsv"

mkdir -p \
    "${DB_DIR}" \
    "${RESULT_DIR}" \
    "${RAW_DIR}" \
    "${QC_DIR}" \
    "${LOG_DIR}" \
    "${ADMIN_DIR}" \
    "${CHECKPOINT_DIR}"

###############################################################################
# REQUIRE STEP 36X0
###############################################################################

X0_CHECKPOINT="${WGDI_ROOT}/checkpoints/step36X0/STEP36X0_COMPLETE.txt"

if [[ ! -s "${X0_CHECKPOINT}" ]]; then
    echo "ERROR: Missing Step 36X0 checkpoint." >&2
    exit 1
fi

if ! grep -q '^status=PASS$' "${X0_CHECKPOINT}"; then
    echo "ERROR: Step 36X0 checkpoint is not PASS." >&2
    exit 1
fi

###############################################################################
# CREATE EXACT X1 MANIFEST
#
# Same schema as original Step 36B:
# task_id, comparison_type, species1, species2
###############################################################################

{
    printf '%s\n' \
        $'task_id\tcomparison_type\tspecies1\tspecies2' \
        $'0\tpairwise\tVSCU\tVSER' \
        $'1\tpairwise\tVSCU\tVPAN' \
        $'2\tpairwise\tVPAN\tVPER'
} > "${X1_MANIFEST}"

###############################################################################
# ENVIRONMENT
###############################################################################

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

for PROGRAM in diamond python gzip awk sort
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: Required program unavailable: ${PROGRAM}" >&2
        exit 1
    fi
done

###############################################################################
# READ ARRAY TASK
###############################################################################

TASK_LINE="$(
    awk -F $'\t' \
        -v task="${SLURM_ARRAY_TASK_ID}" \
        '
        NR > 1 && $1 == task {
            print
            exit
        }
        ' \
        "${X1_MANIFEST}"
)"

if [[ -z "${TASK_LINE}" ]]; then
    echo "ERROR: No comparison for task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

IFS=$'\t' read -r \
    TASK_ID \
    COMPARISON_TYPE \
    SPECIES1 \
    SPECIES2 \
    <<< "${TASK_LINE}"

RUN_ID="${SPECIES1}_${SPECIES2}"

###############################################################################
# INPUTS
###############################################################################

PEP1="${PEP_DIR}/${SPECIES1}.wgdi.pep.fa"
PEP2="${PEP_DIR}/${SPECIES2}.wgdi.pep.fa"

for FILE in \
    "${PEP1}" \
    "${PEP2}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required protein FASTA missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# EXACT ORIGINAL STEP 36B PARAMETERS
###############################################################################

EVALUE="1e-5"
MAX_TARGET_SEQS="100"
SENSITIVE_MODE="--very-sensitive"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

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

###############################################################################
# TARGETED CLEANUP FOR THIS NEW COMPARISON ONLY
###############################################################################

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

        rm -f \
            "${TEMP_DATABASE}.dmnd"

        diamond makedb \
            --in "${PROTEIN_FASTA}" \
            --db "${TEMP_DATABASE}" \
            --threads "${THREADS}"

        if [[ ! -s "${TEMP_DATABASE}.dmnd" ]]; then
            echo "ERROR: DIAMOND database was not created for ${SPECIES}." >&2
            exit 1
        fi

        #
        # Array-safe database creation.
        #
        # Another task can finish the same species database first.
        #
        if [[ -s "${DATABASE}.dmnd" ]]; then
            rm -f \
                "${TEMP_DATABASE}.dmnd"
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
echo "Step 36X1: additional pairwise protein homology"
echo "============================================================"
echo "Task: ${TASK_ID}"
echo "Type: ${COMPARISON_TYPE}"
echo "Species 1: ${SPECIES1}"
echo "Species 2: ${SPECIES2}"
echo "Run ID: ${RUN_ID}"
echo "Start: $(date --iso-8601=seconds)"
echo "DIAMOND: $(diamond version)"
echo "Threads: ${THREADS}"
echo "E-value: ${EVALUE}"
echo "Max target sequences: ${MAX_TARGET_SEQS}"
echo "Sensitivity: ${SENSITIVE_MODE}"

###############################################################################
# FASTA COUNTS
###############################################################################

N_PROTEINS1="$(
    count_fasta_ids \
        "${PEP1}"
)"

N_PROTEINS2="$(
    count_fasta_ids \
        "${PEP2}"
)"

if [[ "${N_PROTEINS1}" -le 0 || "${N_PROTEINS2}" -le 0 ]]; then
    echo "ERROR: Empty protein input." >&2
    exit 1
fi

echo "Species1 proteins: ${N_PROTEINS1}"
echo "Species2 proteins: ${N_PROTEINS2}"

###############################################################################
# DATABASES
###############################################################################

build_database \
    "${SPECIES1}" \
    "${PEP1}"

build_database \
    "${SPECIES2}" \
    "${PEP2}"

###############################################################################
# PAIRWISE RECIPROCAL DIAMOND
###############################################################################

if [[ "${COMPARISON_TYPE}" != "pairwise" ]]; then
    echo "ERROR: Step 36X1 supports pairwise comparisons only." >&2
    exit 1
fi

RAW_FORWARD="${RUN_RAW_DIR}/${SPECIES1}_to_${SPECIES2}.raw.tsv"

RAW_REVERSE="${RUN_RAW_DIR}/${SPECIES2}_to_${SPECIES1}.raw.tsv"

REVERSED_NORMALIZED="${RUN_RAW_DIR}/${SPECIES2}_to_${SPECIES1}.normalized.tmp.tsv"

###############################################################################
# COMMAND PROVENANCE
###############################################################################

cat > "${COMMAND_LOG}" <<EOF2
diamond blastp --query ${PEP1} --db ${DB_DIR}/${SPECIES2} --outfmt 6 --evalue ${EVALUE} --max-target-seqs ${MAX_TARGET_SEQS} --very-sensitive
diamond blastp --query ${PEP2} --db ${DB_DIR}/${SPECIES1} --outfmt 6 --evalue ${EVALUE} --max-target-seqs ${MAX_TARGET_SEQS} --very-sensitive
EOF2

###############################################################################
# FORWARD
###############################################################################

run_diamond \
    "${PEP1}" \
    "${SPECIES2}" \
    "${RAW_FORWARD}"

###############################################################################
# REVERSE
###############################################################################

run_diamond \
    "${PEP2}" \
    "${SPECIES1}" \
    "${RAW_REVERSE}"

if [[ ! -s "${RAW_FORWARD}" ]]; then
    echo "ERROR: Forward DIAMOND result is empty for ${RUN_ID}." >&2
    exit 1
fi

if [[ ! -s "${RAW_REVERSE}" ]]; then
    echo "ERROR: Reverse DIAMOND result is empty for ${RUN_ID}." >&2
    exit 1
fi

###############################################################################
# NORMALIZE REVERSE DIRECTION
#
# Species1 must always be column 1.
# Species2 must always be column 2.
#
# Original columns:
# 1  qseqid
# 2  sseqid
# 3  pident
# 4  length
# 5  mismatch
# 6  gapopen
# 7  qstart
# 8  qend
# 9  sstart
# 10 send
# 11 evalue
# 12 bitscore
###############################################################################

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

if [[ ! -s "${REVERSED_NORMALIZED}" ]]; then
    echo "ERROR: Reverse-normalized table is empty." >&2
    exit 1
fi

###############################################################################
# MERGE RECIPROCAL RESULTS
#
# Keep exactly one row per species1/species2 gene pair.
#
# Priority:
#   1. lowest E-value
#   2. highest bitscore
#   3. highest percentage identity
#
# This reproduces original Step 36B logic.
###############################################################################

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


###############################################################################
# FASTA IDENTIFIERS
###############################################################################

def fasta_ids(path: Path):
    result = set()

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:

        for line in handle:

            if not line.startswith(">"):
                continue

            identifier = (
                line[1:]
                .split()[0]
            )

            if identifier in result:
                raise SystemExit(
                    f"ERROR: Duplicate FASTA ID in {path}: "
                    f"{identifier}"
                )

            result.add(identifier)

    if not result:
        raise SystemExit(
            f"ERROR: No FASTA IDs in {path}"
        )

    return result


ids1 = fasta_ids(
    pep1_path
)

ids2 = fasta_ids(
    pep2_path
)


###############################################################################
# VALIDATE BLAST TABLE
###############################################################################

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


with blast_path.open(
    "r",
    encoding="utf-8",
) as handle:

    for line_number, line in enumerate(
        handle,
        start=1,
    ):

        line = line.rstrip("\n")

        if not line:
            continue

        fields = line.split("\t")

        if len(fields) != 12:
            malformed_rows += 1
            continue

        query_id = fields[0]
        subject_id = fields[1]

        #######################################################################
        # ID VALIDATION
        #######################################################################

        if (
            query_id not in ids1
            or subject_id not in ids2
        ):
            invalid_id_rows += 1

        if (
            comparison_type == "self"
            and query_id == subject_id
        ):
            self_identity_rows += 1

        #######################################################################
        # NUMERIC VALIDATION
        #######################################################################

        try:

            pident = float(
                fields[2]
            )

            length = int(
                fields[3]
            )

            mismatch = int(
                fields[4]
            )

            gapopen = int(
                fields[5]
            )

            qstart = int(
                fields[6]
            )

            qend = int(
                fields[7]
            )

            sstart = int(
                fields[8]
            )

            send = int(
                fields[9]
            )

            evalue = float(
                fields[10]
            )

            bitscore = float(
                fields[11]
            )

        except ValueError:

            malformed_rows += 1
            continue

        #######################################################################
        # VALUE RANGE VALIDATION
        #######################################################################

        numeric_values = [
            pident,
            float(length),
            float(mismatch),
            float(gapopen),
            float(qstart),
            float(qend),
            float(sstart),
            float(send),
            evalue,
            bitscore,
        ]

        if not all(
            math.isfinite(value)
            for value in numeric_values
        ):
            malformed_rows += 1
            continue

        if (
            pident < 0
            or pident > 100
            or length <= 0
        ):
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

        #######################################################################
        # COUNTS
        #######################################################################

        row_count += 1

        unique_queries.add(
            query_id
        )

        unique_subjects.add(
            subject_id
        )

        pair_counter[
            (
                query_id,
                subject_id,
            )
        ] += 1

        identity_sum += pident
        bitscore_sum += bitscore

        if evalue == 0:
            evalue_zero_count += 1


###############################################################################
# DUPLICATE PAIRS
###############################################################################

duplicate_pairs = sum(
    count - 1
    for count in pair_counter.values()
    if count > 1
)


###############################################################################
# HARD QC
###############################################################################

if row_count == 0:
    raise SystemExit(
        "ERROR: No valid homology rows."
    )

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


###############################################################################
# COVERAGE
###############################################################################

query_coverage = (
    len(unique_queries)
    / len(ids1)
)

if comparison_type == "self":
    subject_denominator = len(ids1)
else:
    subject_denominator = len(ids2)

subject_coverage = (
    len(unique_subjects)
    / subject_denominator
)


###############################################################################
# SUMMARY
###############################################################################

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

gzip -f \
    "${FINAL_BLAST}"

if [[ ! -s "${FINAL_BLAST_GZ}" ]]; then
    echo "ERROR: Compressed result was not created." >&2
    exit 1
fi

if [[ ! -s "${TASK_SUMMARY}" ]]; then
    echo "ERROR: Task QC summary was not created." >&2
    exit 1
fi

###############################################################################
# TASK-LEVEL CHECKSUMS
###############################################################################

TASK_CHECKSUM="${RUN_QC_DIR}/${RUN_ID}.homology.sha256"

sha256sum \
    "${FINAL_BLAST_GZ}" \
    "${TASK_SUMMARY}" \
    "${COMMAND_LOG}" \
    > "${TASK_CHECKSUM}"

###############################################################################
# DISPLAY
###############################################################################

echo
echo "============================================================"
echo "Step 36X1 task summary"
echo "============================================================"

column -t -s $'\t' \
    "${TASK_SUMMARY}"

echo
echo "Output:"
echo "${FINAL_BLAST_GZ}"

echo
echo "Checksum:"
cat \
    "${TASK_CHECKSUM}"

echo
echo "Completed: $(date --iso-8601=seconds)"
echo "Step 36X1 task ${TASK_ID}: PASS"
