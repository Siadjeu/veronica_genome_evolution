#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=wgdi36X2A
#SBATCH --output=11_wgdi/logs/step36X2A_%j.out
#SBATCH --error=11_wgdi/logs/step36X2A_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

SOURCE_RUNNER="05_scripts/36C_run_wgdi_improved_collinearity.sh"
NEW_RUNNER="05_scripts/36X2B_run_additional_collinearity.sh"

MANIFEST="${WGDI_ROOT}/00_admin/step36X2_collinearity_manifest.tsv"

QC_DIR="${WGDI_ROOT}/02_qc/additional_pairwise_comparisons"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X2A"

mkdir -p \
    "${WGDI_ROOT}/00_admin" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

###############################################################################
# REQUIRE STEP 36X1 RESULTS
###############################################################################

COMPARISONS=(
    VSCU_VSER
    VSCU_VPAN
    VPAN_VPER
)

for C in "${COMPARISONS[@]}"
do
    FILE="${WGDI_ROOT}/04_homology/results/${C}.blast.tsv.gz"

    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing DIAMOND input: ${FILE}" >&2
        exit 1
    fi

    if ! gzip -t "${FILE}"; then
        echo "ERROR: Invalid gzip file: ${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# VERIFY ORIGINAL SUCCESSFUL RUNNER EXISTS
###############################################################################

if [[ ! -s "${SOURCE_RUNNER}" ]]; then
    echo "ERROR: Missing original Step 36C runner:" >&2
    echo "${SOURCE_RUNNER}" >&2
    exit 1
fi

###############################################################################
# CREATE EXACT SIX-TASK MANIFEST
#
# Columns match original Step 36C:
#
# task_id
# comparison_type
# species1
# species2
# comparison
# parameter_set
# mg
###############################################################################

{
    printf '%s\n' \
        $'task_id\tcomparison_type\tspecies1\tspecies2\tcomparison\tparameter_set\tmg' \
        $'0\tpairwise\tVSCU\tVSER\tVSCU_VSER\tstrict\t25,25' \
        $'1\tpairwise\tVSCU\tVSER\tVSCU_VSER\tsensitive\t50,50' \
        $'2\tpairwise\tVSCU\tVPAN\tVSCU_VPAN\tstrict\t25,25' \
        $'3\tpairwise\tVSCU\tVPAN\tVSCU_VPAN\tsensitive\t50,50' \
        $'4\tpairwise\tVPAN\tVPER\tVPAN_VPER\tstrict\t25,25' \
        $'5\tpairwise\tVPAN\tVPER\tVPAN_VPER\tsensitive\t50,50'
} > "${MANIFEST}"

###############################################################################
# MANIFEST QC
###############################################################################

python - \
    "${MANIFEST}" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])

expected = [
    ["0", "pairwise", "VSCU", "VSER", "VSCU_VSER", "strict", "25,25"],
    ["1", "pairwise", "VSCU", "VSER", "VSCU_VSER", "sensitive", "50,50"],
    ["2", "pairwise", "VSCU", "VPAN", "VSCU_VPAN", "strict", "25,25"],
    ["3", "pairwise", "VSCU", "VPAN", "VSCU_VPAN", "sensitive", "50,50"],
    ["4", "pairwise", "VPAN", "VPER", "VPAN_VPER", "strict", "25,25"],
    ["5", "pairwise", "VPAN", "VPER", "VPAN_VPER", "sensitive", "50,50"],
]

with path.open(newline="", encoding="utf-8") as handle:
    reader = csv.reader(handle, delimiter="\t")
    rows = list(reader)

expected_header = [
    "task_id",
    "comparison_type",
    "species1",
    "species2",
    "comparison",
    "parameter_set",
    "mg",
]

if not rows:
    raise SystemExit("ERROR: Manifest is empty.")

if rows[0] != expected_header:
    raise SystemExit(
        f"ERROR: Incorrect manifest header: {rows[0]}"
    )

if rows[1:] != expected:
    raise SystemExit(
        "ERROR: Manifest rows do not match expected six tasks."
    )

for row_number, row in enumerate(rows, start=1):
    if len(row) != 7:
        raise SystemExit(
            f"ERROR: Manifest row {row_number} has "
            f"{len(row)} columns, expected 7."
        )

print("PASS: 6 tasks, 7 columns, exact comparison definitions.")
PY

###############################################################################
# CLONE ORIGINAL VALIDATED RUNNER
###############################################################################

cp -f \
    "${SOURCE_RUNNER}" \
    "${NEW_RUNNER}"

###############################################################################
# PATCH ONLY WHAT MUST DIFFER
###############################################################################

python - \
    "${NEW_RUNNER}" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])

text = path.read_text()

replacements = [
    (
        "#SBATCH --array=0-23",
        "#SBATCH --array=0-5%6",
    ),
    (
        "#SBATCH --job-name=wgdi36C",
        "#SBATCH --job-name=wgdi36X2B",
    ),
    (
        "#SBATCH --output=11_wgdi/logs/step36C_%A_%a.out",
        "#SBATCH --output=11_wgdi/logs/step36X2B_%A_%a.out",
    ),
    (
        "#SBATCH --error=11_wgdi/logs/step36C_%A_%a.err",
        "#SBATCH --error=11_wgdi/logs/step36X2B_%A_%a.err",
    ),
    (
        'MANIFEST="${WGDI_ROOT}/00_admin/step36C_collinearity_manifest.tsv"',
        'MANIFEST="${WGDI_ROOT}/00_admin/step36X2_collinearity_manifest.tsv"',
    ),
    (
        'echo "Step 36C: WGDI improved collinearity"',
        'echo "Step 36X2B: additional WGDI improved collinearity"',
    ),
    (
        'echo "Step 36C task summary"',
        'echo "Step 36X2B task summary"',
    ),
    (
        'echo "Step 36C task ${TASK_ID}: PASS"',
        'echo "Step 36X2B task ${TASK_ID}: PASS"',
    ),
]

for old, new in replacements:

    count = text.count(old)

    if count != 1:
        raise SystemExit(
            f"ERROR: Expected exactly one occurrence of:\n"
            f"{old}\n"
            f"Found: {count}"
        )

    text = text.replace(
        old,
        new,
        1,
    )

path.write_text(text)

print("PASS: runner cloned and patched.")
PY

chmod +x \
    "${NEW_RUNNER}"

###############################################################################
# SYNTAX CHECK
###############################################################################

bash -n \
    "${NEW_RUNNER}"

###############################################################################
# VERIFY CORE WGDI PARAMETERS DID NOT CHANGE
###############################################################################

python - \
    "${SOURCE_RUNNER}" \
    "${NEW_RUNNER}" <<'PY'
from pathlib import Path
import re
import sys

source = Path(sys.argv[1]).read_text()
new = Path(sys.argv[2]).read_text()

patterns = {
    "comparison": r"comparison = genomes",
    "multiple": r"multiple = 1",
    "evalue": r"evalue = 1e-5",
    "score": r"score = 100",
    "grading": r"grading = 50,30,25",
    "mg": r"mg = \$\{MG\}",
    "pvalue": r"pvalue = 1",
    "repeat_number": r"repeat_number = 20",
    "position": r"positon = order",
    "wgdi_command": r"wgdi -icl",
}

for name, pattern in patterns.items():

    source_match = bool(
        re.search(pattern, source)
    )

    new_match = bool(
        re.search(pattern, new)
    )

    if not source_match:
        raise SystemExit(
            f"ERROR: Parameter {name} missing from original runner."
        )

    if not new_match:
        raise SystemExit(
            f"ERROR: Parameter {name} missing from new runner."
        )

print("PASS: all core WGDI parameters preserved.")
PY

###############################################################################
# RECORD SHA256 OF ORIGINAL AND CLONED RUNNER
###############################################################################

SOURCE_SHA="$(
    sha256sum "${SOURCE_RUNNER}" |
    awk '{print $1}'
)"

NEW_SHA="$(
    sha256sum "${NEW_RUNNER}" |
    awk '{print $1}'
)"

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36X2A_COMPLETE.txt" <<EOF2
checkpoint=step36X2A_prepare_additional_collinearity
date=$(date --iso-8601=seconds)
source_runner=05_scripts/36C_run_wgdi_improved_collinearity.sh
new_runner=05_scripts/36X2B_run_additional_collinearity.sh
source_runner_sha256=${SOURCE_SHA}
new_runner_sha256=${NEW_SHA}
manifest=11_wgdi/00_admin/step36X2_collinearity_manifest.tsv
tasks=6
comparisons=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
parameter_sets=strict,sensitive
strict_mg=25,25
sensitive_mg=50,50
evalue=1e-5
score=100
grading=50,30,25
pvalue=1
repeat_number=20
position_parameter_key=positon
position=order
original_runner_modified=false
status=PASS
next_step=step36X2B_run_additional_collinearity
EOF2

cp -f \
    "${MANIFEST}" \
    "${NEW_RUNNER}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# DISPLAY
###############################################################################

echo
echo "============================================================"
echo "STEP 36X2 MANIFEST"
echo "============================================================"

column -t -s $'\t' \
    "${MANIFEST}"

echo
echo "============================================================"
echo "CORE CONFIGURATION IN NEW RUNNER"
echo "============================================================"

grep -nE \
    'comparison = genomes|multiple = 1|evalue =|score =|grading =|mg =|pvalue =|repeat_number =|positon =|wgdi -icl' \
    "${NEW_RUNNER}"

echo
echo "============================================================"
echo "CHECKPOINT"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X2A_COMPLETE.txt"

echo
echo "Step 36X2A completed successfully."
