#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36G5R2
#SBATCH --output=11_wgdi/logs/step36G5R2_%A_%a.out
#SBATCH --error=11_wgdi/logs/step36G5R2_%A_%a.err
#SBATCH --array=0-2%3

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

MANIFEST="11_wgdi/00_admin/step36G5R_multipeak_manifest.tsv"
TASK_ID="${SLURM_ARRAY_TASK_ID}"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Missing manifest: ${MANIFEST}" >&2
    exit 1
fi

ROW="$(
    awk -F $'\t' \
        -v task="${TASK_ID}" \
        'NR > 1 && $1 == task {
            print
            exit
        }' \
        "${MANIFEST}"
)"

if [[ -z "${ROW}" ]]; then
    echo "ERROR: No manifest row for task ${TASK_ID}." >&2
    exit 1
fi

IFS=$'\t' read -r \
    TASK_FROM_MANIFEST \
    DATASET \
    COMPARISON_COUNT \
    KSFIT \
    CONFIG \
    PDF_OUTPUT \
    SVG_OUTPUT \
    PNG_OUTPUT \
    STATUS \
    <<< "${ROW}"

if [[ "${STATUS}" != "READY" ]]; then
    echo "ERROR: Dataset ${DATASET} is not READY." >&2
    exit 1
fi

for FILE in "${KSFIT}" "${CONFIG}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing input: ${FILE}" >&2
        exit 1
    fi
done

WORK_DIR="11_wgdi/11_ksfigure/10_multipeak_work/${DATASET}"

mkdir -p \
    "${WORK_DIR}" \
    "$(dirname "${PDF_OUTPUT}")" \
    "11_wgdi/logs"

rm -f \
    "${PDF_OUTPUT}" \
    "${SVG_OUTPUT}" \
    "${PNG_OUTPUT}"

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
    echo "ERROR: wgdi is unavailable in wgdi_env." >&2
    exit 1
fi

echo "============================================================"
echo "Step 36G5R2: authentic WGDI ksfigure"
echo "============================================================"
echo "Task ID:          ${TASK_ID}"
echo "Dataset:          ${DATASET}"
echo "Comparison count: ${COMPARISON_COUNT}"
echo "Ks-fit table:     ${KSFIT}"
echo "Base config:      ${CONFIG}"
echo "Command:          wgdi -kf total.conf"

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
            echo "ERROR: Unsupported format: ${FORMAT}" >&2
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
        echo "ERROR: Failed to create run config: ${RUN_CONFIG}" >&2
        exit 1
    fi

    (
        cd "${WORK_DIR}"
        wgdi -kf "total.${FORMAT}.conf"
    )

    if [[ ! -s "${PROJECT_ROOT}/${OUTPUT}" ]]; then
        echo "ERROR: Missing output: ${OUTPUT}" >&2
        exit 1
    fi
done

if [[ "$(head -c 5 "${PDF_OUTPUT}")" != "%PDF-" ]]; then
    echo "ERROR: Invalid PDF signature: ${PDF_OUTPUT}" >&2
    exit 1
fi

if ! grep -q '<svg' "${SVG_OUTPUT}"; then
    echo "ERROR: SVG markup not detected: ${SVG_OUTPUT}" >&2
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

STATUS_FILE="${WORK_DIR}/${DATASET}.status.tsv"

printf \
    'task_id\tdataset\tcomparison_count\tpdf_bytes\tsvg_bytes\tpng_bytes\tstatus\n' \
    > "${STATUS_FILE}"

printf \
    '%s\t%s\t%s\t%s\t%s\t%s\tPASS\n' \
    "${TASK_ID}" \
    "${DATASET}" \
    "${COMPARISON_COUNT}" \
    "$(stat -c '%s' "${PDF_OUTPUT}")" \
    "$(stat -c '%s' "${SVG_OUTPUT}")" \
    "$(stat -c '%s' "${PNG_OUTPUT}")" \
    >> "${STATUS_FILE}"

echo
echo "PASS: ${DATASET}"
echo "PDF: ${PDF_OUTPUT}"
echo "SVG: ${SVG_OUTPUT}"
echo "PNG: ${PNG_OUTPUT}"
