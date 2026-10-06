#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36X3
#SBATCH --output=11_wgdi/logs/step36X3_%j.out
#SBATCH --error=11_wgdi/logs/step36X3_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

SOURCE_SCRIPT="05_scripts/36D1_prepare_wgdi_ks_inputs.sh"
CLONED_SCRIPT="05_scripts/36X3_run_additional_ks_preparation.sh"

NEW_COLLINEARITY_MANIFEST="${WGDI_ROOT}/00_admin/step36X2_wgdi_collinearity_manifest.tsv"

QC_DIR="${WGDI_ROOT}/02_qc/additional_ks_preparation"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X3"

mkdir -p \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

###############################################################################
# REQUIRE PREVIOUS STEP
###############################################################################

PREVIOUS_CHECKPOINT="${WGDI_ROOT}/checkpoints/step36X2C/STEP36X2C_COMPLETE.txt"

if [[ ! -s "${PREVIOUS_CHECKPOINT}" ]]; then
    echo "ERROR: Missing Step 36X2C checkpoint." >&2
    exit 1
fi

if ! grep -q '^status=PASS$' "${PREVIOUS_CHECKPOINT}"; then
    echo "ERROR: Step 36X2C checkpoint is not PASS." >&2
    exit 1
fi

###############################################################################
# VERIFY SOURCE FILES
###############################################################################

if [[ ! -s "${SOURCE_SCRIPT}" ]]; then
    echo "ERROR: Missing original validated Step 36D1 script:" >&2
    echo "${SOURCE_SCRIPT}" >&2
    exit 1
fi

if [[ ! -s "${NEW_COLLINEARITY_MANIFEST}" ]]; then
    echo "ERROR: Missing additional collinearity manifest:" >&2
    echo "${NEW_COLLINEARITY_MANIFEST}" >&2
    exit 1
fi

###############################################################################
# VALIDATE EXTENSION COLLINEARITY MANIFEST
###############################################################################

python - \
    "${NEW_COLLINEARITY_MANIFEST}" <<'PY'
import csv
import sys
from pathlib import Path

path = Path(sys.argv[1])

with path.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

if len(rows) != 6:
    raise SystemExit(
        f"ERROR: Expected 6 collinearity rows; found {len(rows)}."
    )

strict = [
    row
    for row in rows
    if row["parameter_set"] == "strict"
]

if len(strict) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 strict runs; found {len(strict)}."
    )

expected = {
    ("VSCU", "VSER", "VSCU_VSER"),
    ("VSCU", "VPAN", "VSCU_VPAN"),
    ("VPAN", "VPER", "VPAN_VPER"),
}

observed = {
    (
        row["species1"],
        row["species2"],
        row["comparison"],
    )
    for row in strict
}

if observed != expected:
    raise SystemExit(
        f"ERROR: Unexpected strict comparison set: {observed}"
    )

for row in strict:

    if row["mg"] != "25,25":
        raise SystemExit(
            f"ERROR: Strict mg is not 25,25 for "
            f"{row['comparison']}."
        )

    if row["status"] != "PASS":
        raise SystemExit(
            f"ERROR: Status is not PASS for "
            f"{row['comparison']}."
        )

print(
    "PASS: 3 exact strict extension comparisons validated."
)
PY

###############################################################################
# TARGETED CLEANUP OF FAILED PREPARATION ATTEMPT
#
# Only generated X3 preparation files are removed.
# No scientific outputs from previous completed steps are touched.
###############################################################################

rm -f \
    "${CLONED_SCRIPT}" \
    "${CHECKPOINT_DIR}/STEP36X3_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# CLONE ORIGINAL VALIDATED STEP 36D1
###############################################################################

cp \
    "${SOURCE_SCRIPT}" \
    "${CLONED_SCRIPT}"

###############################################################################
# PATCH CLONE
#
# Important:
# We patch exact contexts rather than assuming that checkpoint filenames
# occur only once in the source script.
###############################################################################

python - \
    "${CLONED_SCRIPT}" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])

text = path.read_text(
    encoding="utf-8"
)


def replace_once(old: str, new: str, label: str):
    global text

    count = text.count(old)

    if count != 1:
        raise SystemExit(
            f"ERROR: Patch '{label}' expected exactly one occurrence.\n"
            f"Found: {count}\n"
            f"Target:\n{old}"
        )

    text = text.replace(
        old,
        new,
        1,
    )


###############################################################################
# SLURM METADATA
###############################################################################

replace_once(
    "#SBATCH --job-name=wgdi36D1",
    "#SBATCH --job-name=wgdi36X3R",
    "job name",
)

replace_once(
    "#SBATCH --output=11_wgdi/logs/step36D1_%j.out",
    "#SBATCH --output=11_wgdi/logs/step36X3R_%j.out",
    "stdout log",
)

replace_once(
    "#SBATCH --error=11_wgdi/logs/step36D1_%j.err",
    "#SBATCH --error=11_wgdi/logs/step36X3R_%j.err",
    "stderr log",
)


###############################################################################
# INPUT MANIFEST
###############################################################################

replace_once(
    'COLLINEARITY_MANIFEST="${WGDI_ROOT}/00_admin/wgdi_collinearity_manifest.tsv"',
    'COLLINEARITY_MANIFEST="${WGDI_ROOT}/00_admin/step36X2_wgdi_collinearity_manifest.tsv"',
    "collinearity manifest",
)


###############################################################################
# EXTENSION-SPECIFIC QC / CHECKPOINT / MANIFEST NAMES
###############################################################################

replace_once(
    'QC_DIR="${WGDI_ROOT}/02_qc/ks_preparation"',
    'QC_DIR="${WGDI_ROOT}/02_qc/additional_ks_preparation"',
    "QC directory",
)

replace_once(
    'CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36D1"',
    'CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X3R"',
    "checkpoint directory",
)

replace_once(
    'KS_MANIFEST="${ADMIN_DIR}/step36D_ks_manifest.tsv"',
    'KS_MANIFEST="${ADMIN_DIR}/step36X3_ks_manifest.tsv"',
    "Ks manifest",
)

replace_once(
    'SUMMARY="${QC_DIR}/step36D1_ks_input_summary.tsv"',
    'SUMMARY="${QC_DIR}/step36X3_ks_input_summary.tsv"',
    "Ks preparation summary",
)


###############################################################################
# CLEANUP CHECKPOINT REFERENCE
#
# This is one of the three legitimate STEP36D1_COMPLETE occurrences.
###############################################################################

replace_once(
    '''    "${CHECKPOINT_DIR}/STEP36D1_COMPLETE.txt" \\
    "${CHECKPOINT_DIR}/sha256_checksums.txt"''',
    '''    "${CHECKPOINT_DIR}/STEP36X3R_COMPLETE.txt" \\
    "${CHECKPOINT_DIR}/sha256_checksums.txt"''',
    "checkpoint cleanup reference",
)


###############################################################################
# EXPECT THREE STRICT COMPARISONS, NOT TWELVE
###############################################################################

replace_once(
    "if len(strict_rows) != 12:",
    "if len(strict_rows) != 3:",
    "strict row count condition",
)

replace_once(
    'f"ERROR: Expected 12 strict collinearity runs; "',
    'f"ERROR: Expected 3 strict collinearity runs; "',
    "strict row count error",
)


###############################################################################
# FINAL READY COUNT
###############################################################################

replace_once(
    'if [[ "${READY_COUNT}" -ne 12 ]]; then',
    'if [[ "${READY_COUNT}" -ne 3 ]]; then',
    "READY_COUNT condition",
)

replace_once(
    'echo "ERROR: Only ${READY_COUNT}/12 Ks datasets are ready." >&2',
    'echo "ERROR: Only ${READY_COUNT}/3 Ks datasets are ready." >&2',
    "READY_COUNT error",
)


###############################################################################
# CHECKPOINT CREATION
#
# Second legitimate STEP36D1_COMPLETE occurrence.
###############################################################################

replace_once(
    'cat > "${CHECKPOINT_DIR}/STEP36D1_COMPLETE.txt" <<EOF2',
    'cat > "${CHECKPOINT_DIR}/STEP36X3R_COMPLETE.txt" <<EOF2',
    "checkpoint creation",
)

replace_once(
    "checkpoint=step36D1_prepare_WGDI_Ks_inputs",
    "checkpoint=step36X3R_prepare_additional_WGDI_Ks_inputs",
    "checkpoint name",
)

replace_once(
    "comparisons_expected=12",
    "comparisons_expected=3",
    "checkpoint expected comparisons",
)

replace_once(
    "manifest=11_wgdi/00_admin/step36D_ks_manifest.tsv",
    "manifest=11_wgdi/00_admin/step36X3_ks_manifest.tsv",
    "checkpoint manifest path",
)

replace_once(
    "summary=11_wgdi/02_qc/ks_preparation/step36D1_ks_input_summary.tsv",
    "summary=11_wgdi/02_qc/additional_ks_preparation/step36X3_ks_input_summary.tsv",
    "checkpoint summary path",
)

replace_once(
    "next_step=step36D2_run_WGDI_Ks",
    "next_step=step36X4_prepare_chunked_Ks",
    "checkpoint next step",
)


###############################################################################
# DISPLAY LABELS
###############################################################################

replace_once(
    'echo "Step 36D1 Ks input summary"',
    'echo "Step 36X3R additional Ks input summary"',
    "summary display title",
)

replace_once(
    'echo "Step 36D Ks manifest"',
    'echo "Step 36X3 Ks manifest"',
    "manifest display title",
)

replace_once(
    'echo "Step 36D1 checkpoint"',
    'echo "Step 36X3R checkpoint"',
    "checkpoint display title",
)


###############################################################################
# CHECKPOINT DISPLAY
#
# Third legitimate STEP36D1_COMPLETE occurrence.
###############################################################################

replace_once(
    'cat "${CHECKPOINT_DIR}/STEP36D1_COMPLETE.txt"',
    'cat "${CHECKPOINT_DIR}/STEP36X3R_COMPLETE.txt"',
    "checkpoint display path",
)


###############################################################################
# WRITE PATCHED SCRIPT
###############################################################################

path.write_text(
    text,
    encoding="utf-8",
)

print(
    "PASS: cloned Step 36D1 and applied all "
    "extension-specific patches."
)
PY

###############################################################################
# PERMISSIONS / SYNTAX
###############################################################################

chmod +x \
    "${CLONED_SCRIPT}"

bash -n \
    "${CLONED_SCRIPT}"

###############################################################################
# VERIFY IMPORTANT SCIENTIFIC LOGIC IS IDENTICAL
###############################################################################

python - \
    "${SOURCE_SCRIPT}" \
    "${CLONED_SCRIPT}" <<'PY'
from pathlib import Path
import sys

source = Path(
    sys.argv[1]
).read_text(
    encoding="utf-8"
)

clone = Path(
    sys.argv[2]
).read_text(
    encoding="utf-8"
)

required_fragments = [
    'gene1 = fields[0]',
    'gene2 = fields[2]',
    'if gene1 == gene2:',
    'canonical = tuple(',
    'sorted((gene1, gene2))',
    'if pair in seen:',
    'pairs.sort()',
    'pep_order != cds_order',
    'if gene not in pep_records:',
    'if gene not in cds_records:',
    '"align_software = mafft"',
    '"parameter_set": "strict"',
    '"status": "READY"',
]

for fragment in required_fragments:

    if fragment not in source:
        raise SystemExit(
            f"ERROR: Required logic absent from source: "
            f"{fragment}"
        )

    if fragment not in clone:
        raise SystemExit(
            f"ERROR: Required logic lost from clone: "
            f"{fragment}"
        )

print(
    "PASS: core Ks-input scientific logic preserved."
)
PY

###############################################################################
# VERIFY ALL ORIGINAL D1 OUTPUT REFERENCES WERE REMOVED FROM CLONE
###############################################################################

BAD_REFERENCES="$(
    grep -nE \
        '/00_admin/wgdi_collinearity_manifest\.tsv"|/00_admin/step36D_ks_manifest\.tsv"|/02_qc/ks_preparation/step36D1_ks_input_summary\.tsv"|STEP36D1_COMPLETE\.txt|checkpoints/step36D1"' \
        "${CLONED_SCRIPT}" \
        || true
)"

if [[ -n "${BAD_REFERENCES}" ]]; then

    echo "ERROR: Original Step 36D1 output references remain:" >&2
    echo "${BAD_REFERENCES}" >&2
    exit 1

fi

echo "PASS: no original Step 36D1 output/checkpoint paths remain."

###############################################################################
# VERIFY NEW REFERENCES
###############################################################################

for REQUIRED_TEXT in \
    'step36X2_wgdi_collinearity_manifest.tsv' \
    'step36X3_ks_manifest.tsv' \
    'step36X3_ks_input_summary.tsv' \
    'checkpoints/step36X3R' \
    'STEP36X3R_COMPLETE.txt'
do

    if ! grep -q \
        "${REQUIRED_TEXT}" \
        "${CLONED_SCRIPT}"
    then
        echo "ERROR: Missing new reference: ${REQUIRED_TEXT}" >&2
        exit 1
    fi

done

echo "PASS: all extension-specific paths are present."

###############################################################################
# PREPARATION CHECKPOINT
###############################################################################

SOURCE_SHA="$(
    sha256sum \
        "${SOURCE_SCRIPT}" |
    awk '{print $1}'
)"

CLONE_SHA="$(
    sha256sum \
        "${CLONED_SCRIPT}" |
    awk '{print $1}'
)"

cat > "${CHECKPOINT_DIR}/STEP36X3_COMPLETE.txt" <<EOF2
checkpoint=step36X3_prepare_additional_Ks_runner
date=$(date --iso-8601=seconds)
source_script=05_scripts/36D1_prepare_wgdi_ks_inputs.sh
cloned_script=05_scripts/36X3_run_additional_ks_preparation.sh
source_sha256=${SOURCE_SHA}
clone_sha256=${CLONE_SHA}
collinearity_manifest=11_wgdi/00_admin/step36X2_wgdi_collinearity_manifest.tsv
comparisons_expected=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
parameter_set=strict
mg=25,25
gene_columns=1,3
self_pairs_removed=true
self_reverse_pairs_canonicalized=true
overlapping_pair_assignments_deduplicated=true
protein_cds_presence_required=true
alignment_software=mafft
original_step36D1_modified=false
status=PASS
next_step=run_36X3_run_additional_ks_preparation
EOF2

cp -f \
    "${CLONED_SCRIPT}" \
    "${CHECKPOINT_DIR}/"

find \
    "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# DISPLAY
###############################################################################

echo
echo "============================================================"
echo "STEP 36X3 PREPARATION CHECKPOINT"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X3_COMPLETE.txt"

echo
echo "============================================================"
echo "PATCHED SCRIPT KEY SETTINGS"
echo "============================================================"

grep -nE \
    'COLLINEARITY_MANIFEST=|KS_MANIFEST=|SUMMARY=|strict_rows|READY_COUNT|comparisons_expected=|align_software|STEP36X3R_COMPLETE' \
    "${CLONED_SCRIPT}"

echo
echo "============================================================"
echo "FINAL STATUS"
echo "============================================================"

echo "Step 36X3 preparation: PASS"
