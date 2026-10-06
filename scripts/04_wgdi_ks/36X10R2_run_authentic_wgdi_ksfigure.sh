#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36X10R2
#SBATCH --output=11_wgdi/logs/step36X10R2_%j.out
#SBATCH --error=11_wgdi/logs/step36X10R2_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

MANIFEST="11_wgdi/00_admin/step36X10_multipeak_manifest.tsv"

CHECKPOINT_DIR="11_wgdi/checkpoints/step36X10"
QC_DIR="11_wgdi/02_qc/additional_ksfigure_multipeak"

SUMMARY="${QC_DIR}/step36X10_authentic_ksfigure_summary.tsv"

R1_CHECKPOINT="11_wgdi/checkpoints/step36X10R1/STEP36X10R1_COMPLETE.txt"

###############################################################################
# REQUIRE X10R1
###############################################################################

if [[ ! -s "${R1_CHECKPOINT}" ]]; then
    echo "ERROR: Missing X10R1 checkpoint." >&2
    exit 1
fi

if ! grep -q '^status=PASS$' "${R1_CHECKPOINT}"; then
    echo "ERROR: X10R1 did not PASS." >&2
    exit 1
fi

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Missing X10 manifest: ${MANIFEST}" >&2
    exit 1
fi

mkdir -p \
    "${CHECKPOINT_DIR}" \
    "${QC_DIR}" \
    "11_wgdi/logs"

rm -f \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36X10_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# Read the single prepared dataset robustly with Python
# Avoid the READY\r issue from the failed first X9 implementation.
###############################################################################

mapfile -t VALUES < <(
python - "${MANIFEST}" <<'PY'
import csv
import sys

path = sys.argv[1]

with open(
    path,
    newline="",
    encoding="utf-8-sig",
) as handle:

    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

if len(rows) != 1:
    raise SystemExit(
        f"ERROR: Expected exactly one X10 dataset; "
        f"found {len(rows)}."
    )

row = rows[0]

if row["status"].strip() != "READY":
    raise SystemExit(
        f"ERROR: X10 dataset status is "
        f"{row['status']!r}, not READY."
    )

for key in [
    "dataset",
    "comparison_count",
    "ksfit",
    "config",
    "pdf_output",
    "svg_output",
    "png_output",
]:
    print(row[key].strip())
PY
)

if [[ "${#VALUES[@]}" -ne 7 ]]; then
    echo "ERROR: Failed to parse X10 manifest." >&2
    exit 1
fi

DATASET="${VALUES[0]}"
COMPARISON_COUNT="${VALUES[1]}"
KSFIT="${VALUES[2]}"
CONFIG="${VALUES[3]}"
PDF_OUTPUT="${VALUES[4]}"
SVG_OUTPUT="${VALUES[5]}"
PNG_OUTPUT="${VALUES[6]}"

if [[ "${DATASET}" != "additional_pairwise_comparisons" ]]; then
    echo "ERROR: Unexpected dataset: ${DATASET}" >&2
    exit 1
fi

if [[ "${COMPARISON_COUNT}" -ne 3 ]]; then
    echo "ERROR: Comparison count is ${COMPARISON_COUNT}; expected 3." >&2
    exit 1
fi

for FILE in \
    "${KSFIT}" \
    "${CONFIG}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing X10R2 input: ${FILE}" >&2
        exit 1
    fi
done

WORK_DIR="11_wgdi/11_ksfigure_additional/10_multipeak_work/${DATASET}"

mkdir -p \
    "${WORK_DIR}" \
    "$(dirname "${PDF_OUTPUT}")"

rm -f \
    "${PDF_OUTPUT}" \
    "${SVG_OUTPUT}" \
    "${PNG_OUTPUT}"

###############################################################################
# ENVIRONMENT — identical successful pattern
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

export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

if ! command -v wgdi >/dev/null 2>&1; then
    echo "ERROR: wgdi unavailable in wgdi_env." >&2
    exit 1
fi

echo "============================================================"
echo "STEP 36X10R2: AUTHENTIC WGDI KSFigure"
echo "============================================================"
echo "Dataset:          ${DATASET}"
echo "Comparisons:      ${COMPARISON_COUNT}"
echo "Ksfit:            ${KSFIT}"
echo "Config:           ${CONFIG}"
echo "WGDI command:     wgdi -kf"
echo

###############################################################################
# AUTHENTIC WGDI RUN — PDF/SVG/PNG
###############################################################################

for FORMAT in pdf svg png
do
    case "${FORMAT}" in
        pdf)
            OUTPUT="${PDF_OUTPUT}"
            ;;
        svg)
            OUTPUT="${SVG_OUTPUT}"
            ;;
        png)
            OUTPUT="${PNG_OUTPUT}"
            ;;
        *)
            echo "ERROR: unsupported format ${FORMAT}" >&2
            exit 1
            ;;
    esac

    RUN_CONFIG="${WORK_DIR}/total.${FORMAT}.conf"

    sed \
        -e "s|^ksfit = .*|ksfit = ${PROJECT_ROOT}/${KSFIT}|" \
        -e "s|^savefig = .*|savefig = ${PROJECT_ROOT}/${OUTPUT}|" \
        "${CONFIG}" \
        > "${RUN_CONFIG}"

    if [[ ! -s "${RUN_CONFIG}" ]]; then
        echo "ERROR: Failed to create ${RUN_CONFIG}" >&2
        exit 1
    fi

    echo
    echo "Running authentic WGDI -kf: ${FORMAT}"

    (
        cd "${WORK_DIR}"

        wgdi -kf \
            "total.${FORMAT}.conf"
    )

    if [[ ! -s "${PROJECT_ROOT}/${OUTPUT}" ]]; then
        echo "ERROR: Missing WGDI output: ${OUTPUT}" >&2
        exit 1
    fi
done

###############################################################################
# FILE SIGNATURE VALIDATION
###############################################################################

if [[ "$(head -c 5 "${PDF_OUTPUT}")" != "%PDF-" ]]; then
    echo "ERROR: Invalid PDF signature: ${PDF_OUTPUT}" >&2
    exit 1
fi

if ! grep -q '<svg' "${SVG_OUTPUT}"; then
    echo "ERROR: SVG markup absent: ${SVG_OUTPUT}" >&2
    exit 1
fi

PNG_SIGNATURE="$(
    od -An -t x1 -N 8 "${PNG_OUTPUT}" \
        | tr -d ' \n'
)"

if [[ "${PNG_SIGNATURE}" != "89504e470d0a1a0a" ]]; then
    echo "ERROR: Invalid PNG signature: ${PNG_OUTPUT}" >&2
    exit 1
fi

###############################################################################
# SUMMARY
###############################################################################

printf \
    'dataset\tcomparison_count\tpdf_bytes\tsvg_bytes\tpng_bytes\tstatus\n' \
    > "${SUMMARY}"

printf \
    '%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "${DATASET}" \
    "${COMPARISON_COUNT}" \
    "$(stat -c '%s' "${PDF_OUTPUT}")" \
    "$(stat -c '%s' "${SVG_OUTPUT}")" \
    "$(stat -c '%s' "${PNG_OUTPUT}")" \
    >> "${SUMMARY}"

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36X10_COMPLETE.txt" <<EOF2
checkpoint=step36X10_authentic_additional_WGDI_KsFigure
date=$(date --iso-8601=seconds)
comparisons=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
validated_components=6
formats=PDF,SVG,PNG
wgdi_command=wgdi_-kf
ksfit=11_wgdi/11_ksfigure_additional/06_multipeak_ksfit/additional_pairwise_comparisons.multipeak.ksfit.csv
ksfit_format=WGDI_blank_Gaussian_parameter_headers
curve_normalization=false
y_axis=Kernel_density_of_syntenic_blocks
area=0,3
shadow=true
figsize=11,7
pdf=${PDF_OUTPUT}
svg=${SVG_OUTPUT}
png=${PNG_OUTPUT}
summary=${SUMMARY}
status=PASS
next_step=integrate_BlockKs_KsFigure_and_multisynteny
EOF2

cp -f \
    "${SUMMARY}" \
    "${KSFIT}" \
    "${CONFIG}" \
    "${MANIFEST}" \
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
echo "STEP 36X10 CHECKPOINT"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36X10_COMPLETE.txt"

echo
echo "Step 36X10 completed successfully."
