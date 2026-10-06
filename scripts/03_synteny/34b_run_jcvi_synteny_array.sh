#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-00:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --array=1-28%2
#SBATCH --job-name=jcvi34
#SBATCH --output=10_synteny/jcvi_step34/logs/jcvi34_%A_%a.out
#SBATCH --error=10_synteny/jcvi_step34/logs/jcvi34_%A_%a.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

INPUT_MANIFEST="${SYNTENY_DIR}/manifests/jcvi_input_manifest.tsv"
COMPARISON_MANIFEST="${SYNTENY_DIR}/jcvi_step34/manifests/jcvi_comparisons.tsv"
SETUP_CHECKPOINT="${SYNTENY_DIR}/checkpoint_jcvi_step34_setup/JCVI_STEP34_SETUP_COMPLETE.txt"

STEP_DIR="${SYNTENY_DIR}/jcvi_step34"
WORK_ROOT="${STEP_DIR}/comparisons"
STATUS_DIR="${STEP_DIR}/status"
LOG_DIR="${STEP_DIR}/logs"

mkdir -p \
    "${WORK_ROOT}" \
    "${STATUS_DIR}" \
    "${LOG_DIR}"

cd "${PROJECT_DIR}"

# ============================================================
# Validate prerequisites
# ============================================================

for FILE in \
    "${INPUT_MANIFEST}" \
    "${COMPARISON_MANIFEST}" \
    "${SETUP_CHECKPOINT}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

if ! grep -q '^status=PASS$' "${SETUP_CHECKPOINT}"; then
    echo "ERROR: Step 34A checkpoint is not PASS." >&2
    cat "${SETUP_CHECKPOINT}" >&2
    exit 1
fi

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

echo "============================================================"
echo "Software environment"
echo "============================================================"

echo "Python:"
python --version

echo "JCVI:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY

if ! command -v diamond >/dev/null 2>&1; then
    echo "ERROR: DIAMOND is unavailable after activating jcvi_env." >&2
    echo "Install it with:" >&2
    echo "conda install -n jcvi_env -c conda-forge -c bioconda diamond" >&2
    exit 1
fi

echo "DIAMOND executable:"
command -v diamond

echo "DIAMOND version:"
diamond version

# ============================================================
# Read comparison metadata
# ============================================================

TASK_ENV="${SLURM_TMPDIR:-/tmp}/jcvi34_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.env"

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
    comparison_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

matches = [
    row
    for row in comparison_rows
    if row["task_id"] == task_id
]

if len(matches) != 1:
    raise SystemExit(
        f"ERROR: Expected exactly one comparison for task "
        f"{task_id}; found {len(matches)}."
    )

comparison = matches[0]

with input_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    input_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

inputs = {
    row["species_code"]: row
    for row in input_rows
}

query = comparison["query_species"]
subject = comparison["subject_species"]

for species in [query, subject]:
    if species not in inputs:
        raise SystemExit(
            f"ERROR: Species {species} is missing from the "
            "JCVI input manifest."
        )

values = {
    "TASK_ID": comparison["task_id"],
    "COMPARISON_ID": comparison["comparison_id"],
    "QUERY": query,
    "SUBJECT": subject,
    "COMPARISON_TYPE": comparison["comparison_type"],
    "CSCORE": comparison["cscore"],
    "MINSPAN": comparison["minspan"],
    "PRIORITY": comparison["priority"],
    "RATIONALE": comparison["rationale"],
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
        echo "ERROR: Required comparison input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Prepare isolated workspace
# ============================================================

WORK_DIR="${WORK_ROOT}/${COMPARISON_ID}"
STATUS_FILE="${STATUS_DIR}/${COMPARISON_ID}.status.tsv"

mkdir -p "${WORK_DIR}"
cd "${WORK_DIR}"

rm -f \
    "${QUERY}.bed" \
    "${QUERY}.pep"

ln -s "${QUERY_BED}" "${QUERY}.bed"
ln -s "${QUERY_PEP}" "${QUERY}.pep"

if [[ "${QUERY}" != "${SUBJECT}" ]]; then
    rm -f \
        "${SUBJECT}.bed" \
        "${SUBJECT}.pep"

    ln -s "${SUBJECT_BED}" "${SUBJECT}.bed"
    ln -s "${SUBJECT_PEP}" "${SUBJECT}.pep"
fi

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

echo
echo "============================================================"
echo "JCVI Step 34 comparison"
echo "============================================================"
echo "Task ID:       ${TASK_ID}"
echo "Comparison:    ${COMPARISON_ID}"
echo "Query:         ${QUERY}"
echo "Subject:       ${SUBJECT}"
echo "Type:          ${COMPARISON_TYPE}"
echo "Priority:      ${PRIORITY}"
echo "C-score:       ${CSCORE}"
echo "Minimum span:  ${MINSPAN}"
echo "CPUs:          ${SLURM_CPUS_PER_TASK}"
echo "Working dir:   ${WORK_DIR}"
echo "Start:         ${START_TIME}"

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
comparison_type=${COMPARISON_TYPE}
exit_code=${JCVI_EXIT}
status=FAIL
reason=JCVI_ORTHOLOG_COMMAND_FAILED
EOF2

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" \
        "${COMPARISON_ID}" \
        "${QUERY}" \
        "${SUBJECT}" \
        "${COMPARISON_TYPE}" \
        "FAIL" \
        "JCVI_ORTHOLOG_COMMAND_FAILED" \
        > "${STATUS_FILE}"

    exit "${JCVI_EXIT}"
fi

# ============================================================
# Locate anchors
# ============================================================

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
comparison_type=${COMPARISON_TYPE}
status=FAIL
reason=ANCHOR_FILE_MISSING_OR_EMPTY
EOF2

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" \
        "${COMPARISON_ID}" \
        "${QUERY}" \
        "${SUBJECT}" \
        "${COMPARISON_TYPE}" \
        "FAIL" \
        "ANCHOR_FILE_MISSING_OR_EMPTY" \
        > "${STATUS_FILE}"

    exit 1
fi

# ============================================================
# Count anchor pairs and blocks
# ============================================================

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

UNIQUE_QUERY_ANCHORS=$(
    awk '
        NF >= 2 && $1 !~ /^#/ {
            seen[$1] = 1
        }
        END {
            print length(seen)
        }
    ' "${ANCHOR_FILE}"
)

UNIQUE_SUBJECT_ANCHORS=$(
    awk '
        NF >= 2 && $1 !~ /^#/ {
            seen[$2] = 1
        }
        END {
            print length(seen)
        }
    ' "${ANCHOR_FILE}"
)

if [[ "${ANCHOR_PAIR_COUNT}" -eq 0 ]]; then
    cat > RUN_FAILED.txt <<EOF2
task_id=${TASK_ID}
comparison_id=${COMPARISON_ID}
query=${QUERY}
subject=${SUBJECT}
comparison_type=${COMPARISON_TYPE}
status=FAIL
reason=ZERO_ANCHOR_PAIRS
anchor_file=${WORK_DIR}/${ANCHOR_FILE}
EOF2

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" \
        "${COMPARISON_ID}" \
        "${QUERY}" \
        "${SUBJECT}" \
        "${COMPARISON_TYPE}" \
        "FAIL" \
        "ZERO_ANCHOR_PAIRS" \
        > "${STATUS_FILE}"

    exit 1
fi

# ============================================================
# Optional simplified anchor file
# ============================================================

SIMPLE_ANCHOR_FILE="${ANCHOR_FILE%.anchors}.anchors.simple"
SIMPLE_STATUS="NOT_RUN"
SIMPLE_PAIR_COUNT=0

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

if [[ "${SCREEN_EXIT}" -eq 0 && -s "${SIMPLE_ANCHOR_FILE}" ]]; then
    SIMPLE_STATUS="PASS"

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
else
    SIMPLE_STATUS="OPTIONAL_SCREEN_FAILED"
    rm -f "${SIMPLE_ANCHOR_FILE}"
fi

END_TIME="$(date --iso-8601=seconds)"

cat > RUN_COMPLETE.txt <<EOF2
task_id=${TASK_ID}
comparison_id=${COMPARISON_ID}
query=${QUERY}
subject=${SUBJECT}
comparison_type=${COMPARISON_TYPE}
priority=${PRIORITY}
cscore=${CSCORE}
minspan=${MINSPAN}
align_soft=diamond_blastp
diamond_version=$(diamond version | head -n 1)
anchor_file=${WORK_DIR}/${ANCHOR_FILE}
anchor_pair_count=${ANCHOR_PAIR_COUNT}
anchor_block_count=${ANCHOR_BLOCK_COUNT}
unique_query_anchor_genes=${UNIQUE_QUERY_ANCHORS}
unique_subject_anchor_genes=${UNIQUE_SUBJECT_ANCHORS}
simple_anchor_file=${WORK_DIR}/${SIMPLE_ANCHOR_FILE}
simple_anchor_pair_count=${SIMPLE_PAIR_COUNT}
simple_status=${SIMPLE_STATUS}
start_time=${START_TIME}
end_time=${END_TIME}
status=PASS
EOF2

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${TASK_ID}" \
    "${COMPARISON_ID}" \
    "${QUERY}" \
    "${SUBJECT}" \
    "${COMPARISON_TYPE}" \
    "PASS" \
    "${SIMPLE_STATUS}" \
    > "${STATUS_FILE}"

echo
echo "============================================================"
echo "Comparison completed"
echo "============================================================"
echo "Anchor file:                 ${ANCHOR_FILE}"
echo "Anchor pairs:                ${ANCHOR_PAIR_COUNT}"
echo "Anchor blocks:               ${ANCHOR_BLOCK_COUNT}"
echo "Unique query anchor genes:   ${UNIQUE_QUERY_ANCHORS}"
echo "Unique subject anchor genes: ${UNIQUE_SUBJECT_ANCHORS}"
echo "Simple status:               ${SIMPLE_STATUS}"
echo "Simple pairs:                ${SIMPLE_PAIR_COUNT}"
echo "End:                         ${END_TIME}"
