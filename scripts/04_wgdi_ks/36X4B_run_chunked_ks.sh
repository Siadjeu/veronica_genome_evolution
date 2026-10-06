#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-22%23
#SBATCH --job-name=wgdi36X4B
#SBATCH --output=11_wgdi/logs/step36X4B_%A_%a.out
#SBATCH --error=11_wgdi/logs/step36X4B_%A_%a.err

set -euo pipefail

###############################################################################
# PROJECT
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"
MANIFEST="${WGDI_ROOT}/00_admin/step36X4_chunk_manifest.tsv"

WORK_ROOT="${WGDI_ROOT}/06_ks/05_chunked_additional/work"
LOG_DIR="${WGDI_ROOT}/logs"

mkdir -p \
    "${WORK_ROOT}" \
    "${LOG_DIR}"

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

for PROGRAM in wgdi mafft yn00 python
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: Required program unavailable: ${PROGRAM}" >&2
        exit 1
    fi
done

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Chunk manifest missing: ${MANIFEST}" >&2
    exit 1
fi

###############################################################################
# READ TASK
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
    echo "ERROR: No manifest row for task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

IFS=$'\t' read -r \
    TASK_ID \
    COMPARISON \
    COMPARISON_TYPE \
    SPECIES1 \
    SPECIES2 \
    CHUNK_NUMBER \
    CHUNK_COUNT \
    START_PAIR \
    END_PAIR \
    EXPECTED_PAIRS \
    PAIRS_REL \
    PEP_REL \
    CDS_REL \
    CONFIG_REL \
    RESULT_REL \
    STDOUT_REL \
    STDERR_REL \
    QC_REL \
    STATUS \
    <<< "${TASK_LINE}"

PAIRS="${PROJECT_ROOT}/${PAIRS_REL}"
PEP="${PROJECT_ROOT}/${PEP_REL}"
CDS="${PROJECT_ROOT}/${CDS_REL}"
SOURCE_CONFIG="${PROJECT_ROOT}/${CONFIG_REL}"

RESULT="${PROJECT_ROOT}/${RESULT_REL}"
WGDI_STDOUT="${PROJECT_ROOT}/${STDOUT_REL}"
WGDI_STDERR="${PROJECT_ROOT}/${STDERR_REL}"
QC="${PROJECT_ROOT}/${QC_REL}"

CHUNK_ID="${COMPARISON}.chunk_$(printf '%03d' "${CHUNK_NUMBER}")"
WORK_DIR="${WORK_ROOT}/${CHUNK_ID}"

for FILE in \
    "${PAIRS}" \
    "${PEP}" \
    "${CDS}" \
    "${SOURCE_CONFIG}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required task input missing: ${FILE}" >&2
        exit 1
    fi
done

mkdir -p \
    "$(dirname "${RESULT}")" \
    "$(dirname "${WGDI_STDOUT}")" \
    "$(dirname "${QC}")"

###############################################################################
# REUSE A PREVIOUSLY VALIDATED RESULT
###############################################################################

if [[ -s "${RESULT}" && -s "${QC}" ]]; then
    EXISTING_STATUS="$(
        awk -F $'\t' '
            NR == 1 {
                for (i = 1; i <= NF; i++) {
                    if ($i == "status") {
                        status_col = i
                    }
                }
                next
            }
            NR == 2 && status_col > 0 {
                print $status_col
            }
        ' "${QC}"
    )"

    if [[ "${EXISTING_STATUS}" == "PASS" ]]; then
        echo "Task ${TASK_ID}: existing validated result reused."
        column -t -s $'\t' "${QC}"
        exit 0
    fi
fi

###############################################################################
# ISOLATED WORKING DIRECTORY
###############################################################################

rm -rf "${WORK_DIR}"

mkdir -p "${WORK_DIR}"

LOCAL_CONFIG="${WORK_DIR}/${CHUNK_ID}.ks.conf"
LOCAL_RESULT="${WORK_DIR}/${CHUNK_ID}.ks.tsv"

cat > "${LOCAL_CONFIG}" <<EOF2
[ks]
cds_file = ${CDS}
pep_file = ${PEP}
align_software = mafft
pairs_file = ${PAIRS}
ks_file = ${LOCAL_RESULT}
EOF2

rm -f \
    "${RESULT}" \
    "${WGDI_STDOUT}" \
    "${WGDI_STDERR}" \
    "${QC}"

###############################################################################
# RUN WGDI INSIDE THE ISOLATED DIRECTORY
###############################################################################

echo "============================================================"
echo "Step 36X4B: isolated WGDI Ks chunk"
echo "============================================================"
echo "Task: ${TASK_ID}"
echo "Chunk: ${CHUNK_ID}"
echo "Comparison: ${COMPARISON}"
echo "Expected pairs: ${EXPECTED_PAIRS}"
echo "Work directory: ${WORK_DIR}"
echo "Start: $(date --iso-8601=seconds)"

START_EPOCH="$(date +%s)"

set +e

(
    cd "${WORK_DIR}"

    wgdi -ks "${LOCAL_CONFIG}"
) \
    > "${WGDI_STDOUT}" \
    2> "${WGDI_STDERR}"

WGDI_EXIT_CODE=$?

set -e

END_EPOCH="$(date +%s)"
ELAPSED_SECONDS="$((END_EPOCH - START_EPOCH))"

if [[ "${WGDI_EXIT_CODE}" -ne 0 ]]; then
    echo "ERROR: WGDI -ks exited with code ${WGDI_EXIT_CODE}." >&2

    echo "===== WGDI stderr tail =====" >&2
    tail -n 120 "${WGDI_STDERR}" >&2 || true

    echo "Work directory retained for diagnosis:" >&2
    echo "${WORK_DIR}" >&2

    exit "${WGDI_EXIT_CODE}"
fi

if [[ ! -s "${LOCAL_RESULT}" ]]; then
    echo "ERROR: WGDI produced no result in ${WORK_DIR}." >&2
    exit 1
fi

cp -f \
    "${LOCAL_RESULT}" \
    "${RESULT}"

###############################################################################
# OUTPUT VALIDATION
###############################################################################

python - \
    "${PAIRS}" \
    "${RESULT}" \
    "${EXPECTED_PAIRS}" \
    "${COMPARISON}" \
    "${CHUNK_NUMBER}" \
    "${ELAPSED_SECONDS}" \
    "${QC}" <<'PY'
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

pairs_path = Path(sys.argv[1])
result_path = Path(sys.argv[2])
expected_pairs = int(sys.argv[3])
comparison = sys.argv[4]
chunk_number = int(sys.argv[5])
elapsed_seconds = int(sys.argv[6])
qc_path = Path(sys.argv[7])

expected_columns = [
    "id1",
    "id2",
    "ka_NG86",
    "ks_NG86",
    "ka_YN00",
    "ks_YN00",
]

value_fields = [
    "ka_NG86",
    "ks_NG86",
    "ka_YN00",
    "ks_YN00",
]


def classify_value(raw: str):
    text = "" if raw is None else raw.strip()

    if text == "":
        return "missing", None

    try:
        value = float(text)
    except ValueError:
        return "nonnumeric", None

    if math.isnan(value):
        return "nan", value

    if math.isinf(value):
        return "infinite", value

    return "finite", value


with pairs_path.open(
    encoding="utf-8",
) as handle:
    input_pairs = [
        tuple(line.rstrip("\n").split("\t")[:2])
        for line in handle
        if line.strip()
    ]

if len(input_pairs) != expected_pairs:
    raise SystemExit(
        f"ERROR: Input pair count is {len(input_pairs)}, "
        f"expected {expected_pairs}."
    )

if len(set(input_pairs)) != len(input_pairs):
    raise SystemExit(
        "ERROR: Duplicate input pairs detected."
    )

with result_path.open(
    newline="",
    encoding="utf-8-sig",
    errors="replace",
) as handle:
    reader = csv.DictReader(
        handle,
        delimiter="\t",
    )

    if reader.fieldnames != expected_columns:
        raise SystemExit(
            "ERROR: Unexpected Ks columns: "
            f"{reader.fieldnames}"
        )

    result_rows = list(reader)

result_pairs = [
    (row["id1"], row["id2"])
    for row in result_rows
]

input_pair_set = set(input_pairs)
result_pair_set = set(result_pairs)

duplicate_result_pairs = (
    len(result_pairs)
    - len(result_pair_set)
)

unexpected_pairs = (
    result_pair_set
    - input_pair_set
)

missing_pairs = (
    input_pair_set
    - result_pair_set
)

if duplicate_result_pairs != 0:
    raise SystemExit(
        f"ERROR: {duplicate_result_pairs} duplicate result pairs."
    )

if unexpected_pairs:
    raise SystemExit(
        f"ERROR: {len(unexpected_pairs)} unexpected result pairs."
    )

if missing_pairs:
    raise SystemExit(
        f"ERROR: {len(missing_pairs)} input pairs lack output rows."
    )

if len(result_rows) != expected_pairs:
    raise SystemExit(
        f"ERROR: Result row count is {len(result_rows)}, "
        f"expected {expected_pairs}."
    )

counts = {
    field: {
        "finite": 0,
        "missing": 0,
        "nonnumeric": 0,
        "nan": 0,
        "infinite": 0,
        "negative": 0,
    }
    for field in value_fields
}

rows_with_all_four_finite = 0
rows_with_finite_ks_ng86 = 0
rows_with_finite_ks_yn00 = 0
rows_with_both_ks_finite = 0

for row in result_rows:
    row_classes = {}

    for field in value_fields:
        value_class, value = classify_value(row[field])
        row_classes[field] = value_class
        counts[field][value_class] += 1

        if (
            value_class == "finite"
            and value is not None
            and field.startswith("ks_")
            and value < -1e-12
        ):
            counts[field]["negative"] += 1

    if all(
        row_classes[field] == "finite"
        for field in value_fields
    ):
        rows_with_all_four_finite += 1

    if row_classes["ks_NG86"] == "finite":
        rows_with_finite_ks_ng86 += 1

    if row_classes["ks_YN00"] == "finite":
        rows_with_finite_ks_yn00 += 1

    if (
        row_classes["ks_NG86"] == "finite"
        and row_classes["ks_YN00"] == "finite"
    ):
        rows_with_both_ks_finite += 1

if rows_with_both_ks_finite == 0:
    raise SystemExit(
        "ERROR: No pair has finite NG86 and YN00 Ks estimates."
    )

qc_row = {
    "comparison": comparison,
    "chunk_number": chunk_number,
    "input_pairs": expected_pairs,
    "result_rows": len(result_rows),
    "completion_fraction": (
        f"{len(result_rows) / expected_pairs:.8f}"
    ),
    "duplicate_result_pairs": duplicate_result_pairs,
    "unexpected_result_pairs": len(unexpected_pairs),
    "missing_result_pairs": len(missing_pairs),
    "rows_with_all_four_finite": rows_with_all_four_finite,
    "rows_with_finite_ks_NG86": rows_with_finite_ks_ng86,
    "rows_with_finite_ks_YN00": rows_with_finite_ks_yn00,
    "rows_with_both_ks_finite": rows_with_both_ks_finite,
    "finite_ks_NG86_fraction": (
        f"{rows_with_finite_ks_ng86 / expected_pairs:.8f}"
    ),
    "finite_ks_YN00_fraction": (
        f"{rows_with_finite_ks_yn00 / expected_pairs:.8f}"
    ),
    "both_ks_finite_fraction": (
        f"{rows_with_both_ks_finite / expected_pairs:.8f}"
    ),
    "nonnumeric_ka_NG86": counts["ka_NG86"]["nonnumeric"],
    "nonnumeric_ks_NG86": counts["ks_NG86"]["nonnumeric"],
    "nonnumeric_ka_YN00": counts["ka_YN00"]["nonnumeric"],
    "nonnumeric_ks_YN00": counts["ks_YN00"]["nonnumeric"],
    "nan_ks_NG86": counts["ks_NG86"]["nan"],
    "nan_ks_YN00": counts["ks_YN00"]["nan"],
    "infinite_ks_NG86": counts["ks_NG86"]["infinite"],
    "infinite_ks_YN00": counts["ks_YN00"]["infinite"],
    "negative_finite_ks_NG86": counts["ks_NG86"]["negative"],
    "negative_finite_ks_YN00": counts["ks_YN00"]["negative"],
    "elapsed_seconds": elapsed_seconds,
    "seconds_per_pair": (
        f"{elapsed_seconds / expected_pairs:.8f}"
    ),
    "status": "PASS",
}

with qc_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(qc_row),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerow(qc_row)

print(
    f"{comparison} chunk {chunk_number}: "
    f"{len(result_rows)}/{expected_pairs} rows; "
    f"finite YN00 Ks={rows_with_finite_ks_yn00}; "
    f"finite both Ks={rows_with_both_ks_finite}; "
    f"elapsed={elapsed_seconds}s; PASS"
)
PY

###############################################################################
# SUCCESSFUL CLEANUP
###############################################################################

rm -rf "${WORK_DIR}"

echo
echo "============================================================"
echo "Chunk QC"
echo "============================================================"

column -t -s $'\t' "${QC}"

echo
echo "Completed: $(date --iso-8601=seconds)"
