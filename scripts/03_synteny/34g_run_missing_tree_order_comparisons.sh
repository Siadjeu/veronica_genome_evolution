#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --array=1-3%2
#SBATCH --job-name=jcvi34g
#SBATCH --output=10_synteny/refined_macrosynteny/missing_comparisons/logs/jcvi34g_%A_%a.out
#SBATCH --error=10_synteny/refined_macrosynteny/missing_comparisons/logs/jcvi34g_%A_%a.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

INPUT_MANIFEST="${SYNTENY_DIR}/manifests/jcvi_input_manifest.tsv"

RUN_ROOT="${SYNTENY_DIR}/refined_macrosynteny/missing_comparisons"
COMPARISON_MANIFEST="${RUN_ROOT}/manifests/step34g_comparisons.tsv"
WORK_ROOT="${RUN_ROOT}/comparisons"
STATUS_DIR="${RUN_ROOT}/status"

mkdir -p \
    "${WORK_ROOT}" \
    "${STATUS_DIR}"

cd "${PROJECT_DIR}"

for FILE in \
    "${INPUT_MANIFEST}" \
    "${COMPARISON_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required file is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Activate jcvi_env
# ============================================================

module purge
module load hpc-env/13.1

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
elif [[ -x "${HOME}/anaconda3/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/anaconda3/bin/conda"
else
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate jcvi_env

export MPLBACKEND=Agg

echo "Python:"
python --version

echo "JCVI:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY

if ! command -v diamond >/dev/null 2>&1; then
    echo "ERROR: DIAMOND is unavailable in jcvi_env." >&2
    exit 1
fi

echo "DIAMOND:"
diamond version

# ============================================================
# Read task metadata
# ============================================================

TASK_ENV="${SLURM_TMPDIR:-/tmp}/step34g_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.env"

python - \
    "${COMPARISON_MANIFEST}" \
    "${INPUT_MANIFEST}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${TASK_ENV}" <<'PY'
from __future__ import annotations

import csv
import shlex
import sys
from pathlib import Path

comparison_manifest = Path(sys.argv[1])
input_manifest = Path(sys.argv[2])
task_id = str(sys.argv[3])
output_env = Path(sys.argv[4])

with comparison_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    comparisons = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

selected = [
    row
    for row in comparisons
    if row["task_id"] == task_id
]

if len(selected) != 1:
    raise SystemExit(
        f"ERROR: Expected one comparison for task {task_id}; "
        f"found {len(selected)}."
    )

comparison = selected[0]

with input_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    inputs = {
        row["species_code"]: row
        for row in csv.DictReader(
            handle,
            delimiter="\t",
        )
    }

query = comparison["query_species"]
subject = comparison["subject_species"]

for species in [query, subject]:
    if species not in inputs:
        raise SystemExit(
            f"ERROR: {species} is absent from the JCVI input manifest."
        )

values = {
    "TASK_ID": comparison["task_id"],
    "COMPARISON_ID": comparison["comparison_id"],
    "QUERY": query,
    "SUBJECT": subject,
    "CSCORE": comparison["cscore"],
    "MINSPAN": comparison["minspan"],
    "QUERY_BED": inputs[query]["bed_file"],
    "QUERY_PEP": inputs[query]["protein_fasta"],
    "SUBJECT_BED": inputs[subject]["bed_file"],
    "SUBJECT_PEP": inputs[subject]["protein_fasta"],
}

with output_env.open(
    "w",
    encoding="utf-8",
) as handle:
    for key, value in values.items():
        handle.write(
            f"{key}={shlex.quote(str(value))}\n"
        )
PY

source "${TASK_ENV}"
rm -f "${TASK_ENV}"

for FILE in \
    "${QUERY_BED}" \
    "${QUERY_PEP}" \
    "${SUBJECT_BED}" \
    "${SUBJECT_PEP}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing or empty JCVI input:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Prepare comparison directory
# ============================================================

WORK_DIR="${WORK_ROOT}/${COMPARISON_ID}"
STATUS_FILE="${STATUS_DIR}/${COMPARISON_ID}.status.tsv"

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

rm -f \
    "${QUERY}.bed" \
    "${QUERY}.pep" \
    "${SUBJECT}.bed" \
    "${SUBJECT}.pep"

ln -s "${QUERY_BED}" "${QUERY}.bed"
ln -s "${QUERY_PEP}" "${QUERY}.pep"
ln -s "${SUBJECT_BED}" "${SUBJECT}.bed"
ln -s "${SUBJECT_PEP}" "${SUBJECT}.pep"

rm -f \
    RUN_COMPLETE.txt \
    RUN_FAILED.txt \
    command.txt \
    simple_screen.stdout.txt \
    simple_screen.stderr.txt \
    *.anchors \
    *.anchors.simple \
    *.last \
    *.last.filtered \
    *.blast \
    *.blast.filtered \
    *.dmnd

START_TIME="$(date --iso-8601=seconds)"

echo "============================================================"
echo "Step 34G missing phylogenetic-order comparison"
echo "============================================================"
echo "Task:        ${TASK_ID}"
echo "Comparison:  ${COMPARISON_ID}"
echo "Query:       ${QUERY}"
echo "Subject:     ${SUBJECT}"
echo "C-score:     ${CSCORE}"
echo "Minspan:     ${MINSPAN}"
echo "Start:       ${START_TIME}"

cat > command.txt <<EOF2
python -m jcvi.compara.catalog ortholog \
${QUERY} \
${SUBJECT} \
--dbtype=prot \
--align_soft=diamond_blastp \
--cscore=${CSCORE} \
--no_strip_names \
--cpus=${SLURM_CPUS_PER_TASK}
EOF2

# ============================================================
# Run JCVI
# ============================================================

set +e

python -m jcvi.compara.catalog ortholog \
    "${QUERY}" \
    "${SUBJECT}" \
    --dbtype=prot \
    --align_soft=diamond_blastp \
    --cscore="${CSCORE}" \
    --no_strip_names \
    --cpus="${SLURM_CPUS_PER_TASK}"

JCVI_EXIT=$?

set -e

if [[ "${JCVI_EXIT}" -ne 0 ]]; then
    cat > RUN_FAILED.txt <<EOF2
task_id=${TASK_ID}
comparison_id=${COMPARISON_ID}
query=${QUERY}
subject=${SUBJECT}
exit_code=${JCVI_EXIT}
status=FAIL
reason=JCVI_ORTHOLOG_COMMAND_FAILED
EOF2

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" \
        "${COMPARISON_ID}" \
        "${QUERY}" \
        "${SUBJECT}" \
        "FAIL" \
        "JCVI_ORTHOLOG_COMMAND_FAILED" \
        > "${STATUS_FILE}"

    exit "${JCVI_EXIT}"
fi

EXPECTED_ANCHOR="${QUERY}.${SUBJECT}.anchors"
REVERSE_ANCHOR="${SUBJECT}.${QUERY}.anchors"

if [[ -s "${EXPECTED_ANCHOR}" ]]; then
    ANCHOR_FILE="${EXPECTED_ANCHOR}"
elif [[ -s "${REVERSE_ANCHOR}" ]]; then
    ANCHOR_FILE="${REVERSE_ANCHOR}"
else
    ANCHOR_FILE=$(
        find . \
            -maxdepth 1 \
            -type f \
            -name '*.anchors' \
            -size +0c \
            -printf '%f\n' \
            | sort \
            | head -n 1
    )
fi

if [[ -z "${ANCHOR_FILE:-}" || ! -s "${ANCHOR_FILE}" ]]; then
    cat > RUN_FAILED.txt <<EOF2
task_id=${TASK_ID}
comparison_id=${COMPARISON_ID}
query=${QUERY}
subject=${SUBJECT}
status=FAIL
reason=ANCHOR_FILE_MISSING_OR_EMPTY
EOF2

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" \
        "${COMPARISON_ID}" \
        "${QUERY}" \
        "${SUBJECT}" \
        "FAIL" \
        "ANCHOR_FILE_MISSING_OR_EMPTY" \
        > "${STATUS_FILE}"

    exit 1
fi

ANCHOR_PAIR_COUNT=$(
    awk '
        NF >= 2 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${ANCHOR_FILE}"
)

ANCHOR_BLOCK_COUNT=$(
    awk '
        /^###/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${ANCHOR_FILE}"
)

if [[ "${ANCHOR_PAIR_COUNT}" -eq 0 ]]; then
    echo "ERROR: Zero anchor pairs for ${COMPARISON_ID}." >&2
    exit 1
fi

# ============================================================
# Generate simple anchor file
# ============================================================

SIMPLE_ANCHOR_FILE="${ANCHOR_FILE%.anchors}.anchors.simple"

set +e

python -m jcvi.compara.synteny screen \
    --minspan="${MINSPAN}" \
    --simple \
    "${ANCHOR_FILE}" \
    "${SIMPLE_ANCHOR_FILE}" \
    > simple_screen.stdout.txt \
    2> simple_screen.stderr.txt

SCREEN_EXIT=$?

set -e

if [[ "${SCREEN_EXIT}" -ne 0 || ! -s "${SIMPLE_ANCHOR_FILE}" ]]; then
    cat > RUN_FAILED.txt <<EOF2
task_id=${TASK_ID}
comparison_id=${COMPARISON_ID}
query=${QUERY}
subject=${SUBJECT}
status=FAIL
reason=SIMPLE_ANCHOR_SCREEN_FAILED
anchor_file=${WORK_DIR}/${ANCHOR_FILE}
EOF2

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" \
        "${COMPARISON_ID}" \
        "${QUERY}" \
        "${SUBJECT}" \
        "FAIL" \
        "SIMPLE_ANCHOR_SCREEN_FAILED" \
        > "${STATUS_FILE}"

    exit 1
fi

SIMPLE_PAIR_COUNT=$(
    awk '
        NF >= 2 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${SIMPLE_ANCHOR_FILE}"
)

END_TIME="$(date --iso-8601=seconds)"

cat > RUN_COMPLETE.txt <<EOF2
task_id=${TASK_ID}
comparison_id=${COMPARISON_ID}
query=${QUERY}
subject=${SUBJECT}
cscore=${CSCORE}
minspan=${MINSPAN}
align_soft=diamond_blastp
anchor_file=${WORK_DIR}/${ANCHOR_FILE}
anchor_pair_count=${ANCHOR_PAIR_COUNT}
anchor_block_count=${ANCHOR_BLOCK_COUNT}
simple_anchor_file=${WORK_DIR}/${SIMPLE_ANCHOR_FILE}
simple_anchor_pair_count=${SIMPLE_PAIR_COUNT}
start_time=${START_TIME}
end_time=${END_TIME}
status=PASS
EOF2

printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${TASK_ID}" \
    "${COMPARISON_ID}" \
    "${QUERY}" \
    "${SUBJECT}" \
    "PASS" \
    "PASS" \
    > "${STATUS_FILE}"

echo
echo "============================================================"
echo "Step 34G comparison completed"
echo "============================================================"
echo "Anchor file:         ${ANCHOR_FILE}"
echo "Anchor pairs:        ${ANCHOR_PAIR_COUNT}"
echo "Anchor blocks:       ${ANCHOR_BLOCK_COUNT}"
echo "Simple anchor file:  ${SIMPLE_ANCHOR_FILE}"
echo "Simple anchor pairs: ${SIMPLE_PAIR_COUNT}"
