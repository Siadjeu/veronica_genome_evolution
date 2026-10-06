#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36E2
#SBATCH --output=11_wgdi/logs/step36E2_%A_%a.out
#SBATCH --error=11_wgdi/logs/step36E2_%A_%a.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"
MANIFEST="${WGDI_ROOT}/00_admin/step36E_blockinfo_manifest.tsv"

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

if ! command -v wgdi >/dev/null 2>&1; then
    echo "ERROR: WGDI is unavailable." >&2
    exit 1
fi

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Missing blockinfo manifest: ${MANIFEST}" >&2
    exit 1
fi

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
    echo "ERROR: No row for array task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

IFS=$'\t' read -r \
    TASK_ID \
    COMPARISON \
    COMPARISON_TYPE \
    SPECIES1 \
    SPECIES2 \
    BLAST_REL \
    GFF1_REL \
    GFF2_REL \
    LENS1_REL \
    LENS2_REL \
    COLLINEARITY_REL \
    KS_REL \
    KS_COL \
    CONFIG_REL \
    RESULT_REL \
    STDOUT_REL \
    STDERR_REL \
    QC_REL \
    WORK_REL \
    STATUS \
    <<< "${TASK_LINE}"

CONFIG="${PROJECT_ROOT}/${CONFIG_REL}"
RESULT="${PROJECT_ROOT}/${RESULT_REL}"
WGDI_STDOUT="${PROJECT_ROOT}/${STDOUT_REL}"
WGDI_STDERR="${PROJECT_ROOT}/${STDERR_REL}"
QC="${PROJECT_ROOT}/${QC_REL}"
WORK_DIR="${PROJECT_ROOT}/${WORK_REL}"

if [[ ! -s "${CONFIG}" ]]; then
    echo "ERROR: Missing configuration: ${CONFIG}" >&2
    exit 1
fi

mkdir -p \
    "$(dirname "${RESULT}")" \
    "$(dirname "${QC}")"

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
        echo "${COMPARISON}: existing PASS result reused."
        column -t -s $'\t' "${QC}"
        exit 0
    fi
fi

rm -rf "${WORK_DIR}"
mkdir -p "${WORK_DIR}"

rm -f \
    "${RESULT}" \
    "${WGDI_STDOUT}" \
    "${WGDI_STDERR}" \
    "${QC}"

START_EPOCH="$(date +%s)"

set +e

(
    cd "${WORK_DIR}"
    wgdi -bi "${CONFIG}"
) \
    > "${WGDI_STDOUT}" \
    2> "${WGDI_STDERR}"

WGDI_EXIT_CODE=$?

set -e

END_EPOCH="$(date +%s)"
ELAPSED_SECONDS="$((END_EPOCH - START_EPOCH))"

if [[ "${WGDI_EXIT_CODE}" -ne 0 ]]; then
    echo "ERROR: WGDI -bi exited ${WGDI_EXIT_CODE} for ${COMPARISON}." >&2

    echo "===== WGDI stderr tail =====" >&2
    tail -n 120 "${WGDI_STDERR}" >&2 || true

    echo "===== WGDI stdout tail =====" >&2
    tail -n 120 "${WGDI_STDOUT}" >&2 || true

    echo "Work directory retained: ${WORK_DIR}" >&2
    exit "${WGDI_EXIT_CODE}"
fi

if [[ ! -s "${RESULT}" ]]; then
    echo "ERROR: WGDI produced no blockinfo CSV for ${COMPARISON}." >&2
    exit 1
fi

python - \
    "${RESULT}" \
    "${COMPARISON}" \
    "${COMPARISON_TYPE}" \
    "${SPECIES1}" \
    "${SPECIES2}" \
    "${KS_COL}" \
    "${ELAPSED_SECONDS}" \
    "${QC}" <<'PY'
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

result_path = Path(sys.argv[1])
comparison = sys.argv[2]
comparison_type = sys.argv[3]
species1 = sys.argv[4]
species2 = sys.argv[5]
ks_column_requested = sys.argv[6]
elapsed_seconds = int(sys.argv[7])
qc_path = Path(sys.argv[8])

with result_path.open(
    newline="",
    encoding="utf-8-sig",
    errors="replace",
) as handle:
    reader = csv.DictReader(handle)
    fieldnames = reader.fieldnames
    rows = list(reader)

if not fieldnames:
    raise SystemExit(
        f"ERROR: Blockinfo output lacks a CSV header: {result_path}"
    )

if not rows:
    raise SystemExit(
        f"ERROR: Blockinfo output contains no data rows: {result_path}"
    )

normalized_fields = {
    field.strip(): field
    for field in fieldnames
    if field is not None
}

ks_candidates = [
    field
    for field in fieldnames
    if field is not None
    and "ks" in field.lower()
]

if not ks_candidates:
    raise SystemExit(
        "ERROR: No Ks-related column detected in blockinfo output. "
        f"Columns: {fieldnames}"
    )

completely_empty_rows = 0
rows_with_any_ks_value = 0
finite_scalar_ks_values = 0
nonfinite_scalar_ks_values = 0
nonnumeric_scalar_ks_values = 0

for row in rows:
    if not any(
        value is not None and value.strip()
        for value in row.values()
    ):
        completely_empty_rows += 1
        continue

    row_has_ks = False

    for column in ks_candidates:
        raw = row.get(column)

        if raw is None or not raw.strip():
            continue

        row_has_ks = True
        text = raw.strip()

        try:
            value = float(text)
        except ValueError:
            # WGDI may store lists of Ks values in some columns.
            nonnumeric_scalar_ks_values += 1
            continue

        if math.isfinite(value):
            finite_scalar_ks_values += 1
        else:
            nonfinite_scalar_ks_values += 1

    if row_has_ks:
        rows_with_any_ks_value += 1

if completely_empty_rows != 0:
    raise SystemExit(
        f"ERROR: {completely_empty_rows} completely empty blockinfo rows."
    )

qc_row = {
    "comparison": comparison,
    "comparison_type": comparison_type,
    "species1": species1,
    "species2": species2,
    "ks_column_requested": ks_column_requested,
    "csv_columns": len(fieldnames),
    "blockinfo_rows": len(rows),
    "ks_related_columns": ",".join(ks_candidates),
    "rows_with_any_ks_value": rows_with_any_ks_value,
    "finite_scalar_ks_values": finite_scalar_ks_values,
    "nonfinite_scalar_ks_values": nonfinite_scalar_ks_values,
    "nonnumeric_or_list_ks_cells": nonnumeric_scalar_ks_values,
    "elapsed_seconds": elapsed_seconds,
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
    f"{comparison}: {len(rows)} blockinfo rows, "
    f"{len(fieldnames)} columns, PASS"
)
PY

rm -rf "${WORK_DIR}"

echo
column -t -s $'\t' "${QC}"
