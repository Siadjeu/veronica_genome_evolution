#!/usr/bin/env bash
#SBATCH --job-name=wgdi36X10
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --output=11_wgdi/logs/step36X10_%j.out
#SBATCH --error=11_wgdi/logs/step36X10_%j.err

set -euo pipefail

###############################################################################
# ENVIRONMENT
###############################################################################

if [[ -z "${CONDA_EXE:-}" ]]; then
    if command -v conda >/dev/null 2>&1; then
        CONDA_BASE="$(conda info --base)"
    else
        echo "ERROR: conda not found in PATH." >&2
        exit 1
    fi
else
    CONDA_BASE="$(dirname "$(dirname "${CONDA_EXE}")")"
fi

# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate wgdi_env

###############################################################################
# PATHS
###############################################################################

WGDI_ROOT="11_wgdi"
LOG_DIR="${WGDI_ROOT}/logs"

BASE_R1="05_scripts/36G5R1_prepare_multipeak_wgdi_ksfit.sh"
DERIVED_R1="05_scripts/36X10R1_prepare_additional_multipeak_wgdi_ksfit.sh"

X7_VALIDATED="${WGDI_ROOT}/08_block_ks/09_additional_corrected_peak_validation/tables/step36X7_validated_peaks.tsv"
STANDARDIZED_DIR="${WGDI_ROOT}/11_ksfigure_additional/00_standardized_inputs"
STANDARDIZED_PEAKS="${STANDARDIZED_DIR}/step36X7_validated_peaks.standardized.tsv"

KSFIG_ROOT="${WGDI_ROOT}/11_ksfigure_additional"
KSFIT_DIR="${KSFIG_ROOT}/06_multipeak_ksfit"
CONFIG_DIR="${KSFIG_ROOT}/07_multipeak_configs"
PLOT_DIR="${KSFIG_ROOT}/08_final_plots"

QC_DIR="${WGDI_ROOT}/02_qc/additional_ksfigure_multipeak"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X10"
MANIFEST="${WGDI_ROOT}/00_admin/step36X10_ksfigure_manifest.tsv"

mkdir -p \
    "${LOG_DIR}" \
    "${STANDARDIZED_DIR}" \
    "${KSFIT_DIR}" \
    "${CONFIG_DIR}" \
    "${PLOT_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/00_admin"

rm -f \
    "${STANDARDIZED_PEAKS}" \
    "${MANIFEST}" \
    "${CONFIG_DIR}"/*.conf \
    "${PLOT_DIR}"/* \
    "${CHECKPOINT_DIR}/STEP36X10_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# STEP 1: STANDARDIZE X7 VALIDATED PEAKS
###############################################################################

python - "${X7_VALIDATED}" "${STANDARDIZED_PEAKS}" "${MANIFEST}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

x7 = Path(sys.argv[1])
out = Path(sys.argv[2])
manifest = Path(sys.argv[3])

if not x7.is_file():
    raise SystemExit(f"ERROR: Missing Step 36X7 validated peaks file: {x7}")

rows = []
comparisons = set()

with x7.open() as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    required = {
        "comparison", "comparison_type", "species1", "species2",
        "peak_rank", "primary_peak_ks", "validation_class"
    }
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise SystemExit(
            f"ERROR: Step 36X7 validated peaks header missing required columns. "
            f"Header: {reader.fieldnames}"
        )

    for rec in reader:
        comparison = rec["comparison"]
        comparisons.add(comparison)
        component_id = f"{comparison}_P{int(rec['peak_rank']):02d}"
        rows.append(
            {
                "comparison": comparison,
                "comparison_type": rec["comparison_type"],
                "species1": rec["species1"],
                "species2": rec["species2"],
                "peak_rank": rec["peak_rank"],
                "component_id": component_id,
                "peak_ks": rec["primary_peak_ks"],
                "validation_class": rec["validation_class"],
            }
        )

if len(comparisons) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 additional comparisons, found {len(comparisons)}."
    )

with out.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=[
            "comparison",
            "comparison_type",
            "species1",
            "species2",
            "peak_rank",
            "component_id",
            "peak_ks",
            "validation_class",
        ],
    )
    writer.writeheader()
    writer.writerows(rows)

with manifest.open("w", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        delimiter="\t",
        fieldnames=["comparison", "comparison_type", "species1", "species2", "status"],
    )
    writer.writeheader()
    for comparison in sorted(comparisons):
        first = next(x for x in rows if x["comparison"] == comparison)
        writer.writerow(
            {
                "comparison": comparison,
                "comparison_type": first["comparison_type"],
                "species1": first["species1"],
                "species2": first["species2"],
                "status": "READY",
            }
        )

print(f"Standardized {len(rows)} validated peaks across {len(comparisons)} comparisons.")
PY

###############################################################################
# STEP 2: DERIVE A WORKING ADDITIONAL VERSION OF 36G5R1
###############################################################################

if [[ ! -s "${BASE_R1}" ]]; then
    echo "ERROR: Missing verified base script: ${BASE_R1}" >&2
    exit 1
fi

python - "${BASE_R1}" "${DERIVED_R1}" "${STANDARDIZED_PEAKS}" <<'PY'
from __future__ import annotations

import sys
from pathlib import Path

base = Path(sys.argv[1])
derived = Path(sys.argv[2])
standardized_peaks = sys.argv[3]

text = base.read_text()

replacements = [
    (
        'VALIDATED_PEAKS="${WGDI_ROOT}/11_ksfigure/00_standardized_inputs/step36F4R_validated_peaks.standardized.tsv"',
        f'VALIDATED_PEAKS="{standardized_peaks}"',
    ),
    (
        '11_wgdi/11_ksfigure/06_multipeak_ksfit',
        '11_wgdi/11_ksfigure_additional/06_multipeak_ksfit',
    ),
    (
        '11_wgdi/11_ksfigure/07_multipeak_configs',
        '11_wgdi/11_ksfigure_additional/07_multipeak_configs',
    ),
    (
        '11_wgdi/11_ksfigure/09_multipeak_diagnostics',
        '11_wgdi/11_ksfigure_additional/09_multipeak_diagnostics',
    ),
    (
        '11_wgdi/02_qc/ksfigure_multipeak',
        '11_wgdi/02_qc/additional_ksfigure_multipeak',
    ),
    (
        '11_wgdi/checkpoints/step36G5R1',
        '11_wgdi/checkpoints/step36X10R1',
    ),
    (
        '11_wgdi/00_admin/step36G5R_multipeak_manifest.tsv',
        '11_wgdi/00_admin/step36X10R_multipeak_manifest.tsv',
    ),
    (
        "step36G5R1_prepare_multipeak_WGDI_KsFigure",
        "step36X10R1_prepare_additional_multipeak_WGDI_KsFigure",
    ),
    (
        "step36G5R1",
        "step36X10R1",
    ),
    (
        "if len(manifest_rows) != 12:",
        "if len(manifest_rows) != 3:",
    ),
    (
        'f"ERROR: Expected 12 comparisons; found {len(manifest_rows)}."',
        'f"ERROR: Expected 3 comparisons; found {len(manifest_rows)}."',
    ),
]

for old, new in replacements:
    if old not in text:
        raise SystemExit(
            "ERROR: required replacement target not found in base script:\n"
            f"{old}"
        )
    text = text.replace(old, new)

derived.write_text(text)
print(f"Derived script written to {derived}")
PY

chmod +x "${DERIVED_R1}"
bash -n "${DERIVED_R1}"

###############################################################################
# STEP 3: RUN DERIVED MULTIPEAK PREPARATION
###############################################################################

bash "${DERIVED_R1}"

###############################################################################
# STEP 4: PREPARE AUTHENTIC WGDI -kf CONFIGS
###############################################################################

KSFIT_FILE=""
if [[ -s "${KSFIT_DIR}/pairwise_comparisons.multipeak.ksfit.csv" ]]; then
    KSFIT_FILE="${KSFIT_DIR}/pairwise_comparisons.multipeak.ksfit.csv"
elif [[ -s "${KSFIT_DIR}/all_comparisons.multipeak.ksfit.csv" ]]; then
    KSFIT_FILE="${KSFIT_DIR}/all_comparisons.multipeak.ksfit.csv"
else
    echo "ERROR: No ksfit CSV found in ${KSFIT_DIR}" >&2
    exit 1
fi

SUMMARY="${QC_DIR}/step36X10_ksfigure_summary.tsv"
printf '%s\n' \
    "figure_set	format	config	ksfit	savefig	exists	size_bytes	status" \
    > "${SUMMARY}"

for FMT in pdf svg png
do
    CONFIG="${CONFIG_DIR}/additional_pairwise.${FMT}.total.conf"
    SAVEFIG="${PLOT_DIR}/additional_pairwise_ksfigure.${FMT}"

    cat > "${CONFIG}" <<EOF2
[ksfigure]
ksfit = ${KSFIT_FILE}
labelfontsize = 15
legendfontsize = 15
xlabel = Ks
ylabel = kernel density of syntenic blocks
title = Additional pairwise Veronica comparisons
area = 0,3
figsize = 10,6.18
shadow = true
savefig = ${SAVEFIG}
EOF2

    wgdi -kf "${CONFIG}" \
        > "${PLOT_DIR}/additional_pairwise_ksfigure.${FMT}.stdout.txt" \
        2> "${PLOT_DIR}/additional_pairwise_ksfigure.${FMT}.stderr.txt"

    if [[ -s "${SAVEFIG}" ]]; then
        EXISTS="YES"
        SIZE="$(stat -c '%s' "${SAVEFIG}")"
        STATUS="PASS"
    else
        EXISTS="NO"
        SIZE="0"
        STATUS="FAIL"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "additional_pairwise" \
        "${FMT}" \
        "${CONFIG}" \
        "${KSFIT_FILE}" \
        "${SAVEFIG}" \
        "${EXISTS}" \
        "${SIZE}" \
        "${STATUS}" \
        >> "${SUMMARY}"
done

###############################################################################
# STEP 5: CHECKPOINT
###############################################################################

PASS_COUNT="$(
    awk -F $'\t' 'NR>1 && $8=="PASS"{n++} END{print n+0}' "${SUMMARY}"
)"

if [[ "${PASS_COUNT}" -ne 3 ]]; then
    echo "ERROR: ${PASS_COUNT}/3 KsFigure outputs passed." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36X10_COMPLETE.txt" <<EOF2
checkpoint=step36X10_additional_WGDI_ksfigure
date=$(date --iso-8601=seconds)
comparisons=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
validated_peak_source=step36X7
standardized_peak_input=11_wgdi/11_ksfigure_additional/00_standardized_inputs/step36X7_validated_peaks.standardized.tsv
ksfit_directory=11_wgdi/11_ksfigure_additional/06_multipeak_ksfit
config_directory=11_wgdi/11_ksfigure_additional/07_multipeak_configs
plot_directory=11_wgdi/11_ksfigure_additional/08_final_plots
wgdi_command=wgdi_-kf
summary=11_wgdi/02_qc/additional_ksfigure_multipeak/step36X10_ksfigure_summary.tsv
status=PASS
next_step=interpret_additional_pairwise_KsFigure_and_integrate_with_multisynteny
EOF2

cp -f \
    "${STANDARDIZED_PEAKS}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "============================================================"
echo "STEP 36X10 SUMMARY"
echo "============================================================"
column -t -s $'\t' "${SUMMARY}"

echo
echo "============================================================"
echo "CHECKPOINT"
echo "============================================================"
cat "${CHECKPOINT_DIR}/STEP36X10_COMPLETE.txt"

echo
echo "Step 36X10 completed successfully."
