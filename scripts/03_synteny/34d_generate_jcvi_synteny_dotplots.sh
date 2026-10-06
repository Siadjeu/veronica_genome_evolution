#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=jcvi34_plot
#SBATCH --output=10_synteny/logs/jcvi34_plot_%j.out
#SBATCH --error=10_synteny/logs/jcvi34_plot_%j.err

set -euo pipefail

# ============================================================
# Project paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

STEP34_CHECKPOINT="${SYNTENY_DIR}/checkpoint_jcvi_step34/JCVI_STEP34_COMPLETE.txt"
STEP34_RESULTS="${SYNTENY_DIR}/jcvi_step34/tables/jcvi_step34_results.tsv"
JCVI_MANIFEST="${SYNTENY_DIR}/manifests/jcvi_input_manifest.tsv"

PLOT_ROOT="${SYNTENY_DIR}/jcvi_step34/plots"
PDF_DIR="${PLOT_ROOT}/pdf"
PNG_DIR="${PLOT_ROOT}/png"
LOG_DIR="${PLOT_ROOT}/logs"
TABLE_DIR="${PLOT_ROOT}/tables"
WORK_DIR="${PLOT_ROOT}/work"

CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_jcvi_step34_plots"

mkdir -p \
    "${PDF_DIR}" \
    "${PNG_DIR}" \
    "${LOG_DIR}" \
    "${TABLE_DIR}" \
    "${WORK_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${SYNTENY_DIR}/logs"

cd "${PROJECT_DIR}"

# ============================================================
# Validate inputs
# ============================================================

for FILE in \
    "${STEP34_CHECKPOINT}" \
    "${STEP34_RESULTS}" \
    "${JCVI_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

if ! grep -q '^status=PASS$' "${STEP34_CHECKPOINT}"; then
    echo "ERROR: Step 34 checkpoint is not PASS." >&2
    cat "${STEP34_CHECKPOINT}" >&2
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

export MPLBACKEND=Agg

echo "Python:"
python --version

echo "JCVI:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY

# ============================================================
# Determine supported JCVI plotting options
# ============================================================

HELP_FILE="${WORK_DIR}/jcvi_dotplot_help.txt"

python -m jcvi.graphics.dotplot --help \
    > "${HELP_FILE}" \
    2>&1 || true

if ! grep -q -- '--qbed' "${HELP_FILE}"; then
    echo "ERROR: JCVI dotplot help does not contain --qbed." >&2
    cat "${HELP_FILE}" >&2
    exit 1
fi

EXTRA_OPTIONS=()

if grep -q -- '--notex' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--notex")
fi

if grep -q -- '--skipempty' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--skipempty")
fi

if grep -q -- '--nostdpf' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--nostdpf")
fi

if grep -q -- '--nochpf' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--nochpf")
fi

# ============================================================
# Clean old plotting products
# ============================================================

rm -f "${PDF_DIR}"/*
rm -f "${PNG_DIR}"/*
rm -f "${LOG_DIR}"/*
rm -f "${TABLE_DIR}"/*
rm -f "${CHECKPOINT_DIR}"/*

# ============================================================
# Create the selected plot manifest safely
# ============================================================

PLOT_MANIFEST="${TABLE_DIR}/selected_dotplots.tsv"

python - \
    "${STEP34_RESULTS}" \
    "${JCVI_MANIFEST}" \
    "${PLOT_MANIFEST}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

step34_results_file = Path(sys.argv[1])
jcvi_manifest_file = Path(sys.argv[2])
output_file = Path(sys.argv[3])

with step34_results_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    result_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

with jcvi_manifest_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    manifest_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

results = {
    row["comparison_id"]: row
    for row in result_rows
}

species = {
    row["species_code"]: row
    for row in manifest_rows
}

selected = [
    # All self-synteny plots
    ("VPAN__VPAN", "self", "VPAN self-synteny"),
    ("VSCU__VSCU", "self", "VSCU self-synteny"),
    ("VANA__VANA", "self", "VANA self-synteny"),
    ("VARV__VARV", "self", "VARV self-synteny"),
    ("VPER__VPER", "self", "VPER self-synteny"),
    ("VSER__VSER", "self", "VSER self-synteny"),
    ("VTRI__VTRI", "self", "VTRI self-synteny"),
    ("VVER__VVER", "self", "VVER self-synteny"),
    ("PMAJ__PMAJ", "self", "PMAJ self-synteny"),

    # Priority interspecies plots
    ("VANA__VPER", "pairwise", "VANA versus VPER"),
    ("VANA__VPAN", "pairwise", "VANA versus VPAN"),
    ("VANA__VSCU", "pairwise", "VANA versus VSCU"),
    ("VPER__VSER", "pairwise", "VPER versus VSER"),
    ("VPER__VVER", "pairwise", "VPER versus VVER"),
    ("VPAN__VSCU", "pairwise", "VPAN versus VSCU"),
    ("PMAJ__VANA", "outgroup", "PMAJ versus VANA"),
    ("PMAJ__VPER", "outgroup", "PMAJ versus VPER"),
]

output_rows = []

for plot_number, (
    comparison_id,
    plot_group,
    title,
) in enumerate(selected, start=1):

    if comparison_id not in results:
        raise SystemExit(
            f"ERROR: Comparison missing from Step 34: "
            f"{comparison_id}"
        )

    result = results[comparison_id]

    if result["status"] != "PASS":
        raise SystemExit(
            f"ERROR: Comparison is not PASS: "
            f"{comparison_id}"
        )

    query = result["query_species"]
    subject = result["subject_species"]

    if query not in species or subject not in species:
        raise SystemExit(
            f"ERROR: Missing JCVI input metadata for "
            f"{comparison_id}"
        )

    anchor_file = Path(result["anchor_file"])
    query_bed = Path(species[query]["bed_file"])
    subject_bed = Path(species[subject]["bed_file"])

    for path in [
        anchor_file,
        query_bed,
        subject_bed,
    ]:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: Missing or empty plotting input: "
                f"{path}"
            )

    output_rows.append(
        {
            "plot_number": plot_number,
            "comparison_id": comparison_id,
            "plot_group": plot_group,
            "query_species": query,
            "subject_species": subject,
            "anchor_file": str(anchor_file),
            "query_bed": str(query_bed),
            "subject_bed": str(subject_bed),
            "anchor_pair_count": result[
                "anchor_pair_count"
            ],
            "title": title,
        }
    )

with output_file.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "plot_number",
            "comparison_id",
            "plot_group",
            "query_species",
            "subject_species",
            "anchor_file",
            "query_bed",
            "subject_bed",
            "anchor_pair_count",
            "title",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(output_rows)

print("Selected dotplots:", len(output_rows))
PY

# ============================================================
# Generate the plots
# ============================================================

QC_TABLE="${TABLE_DIR}/dotplot_qc.tsv"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "comparison_id" \
    "plot_group" \
    "anchor_pair_count" \
    "pdf_file" \
    "png_file" \
    "pdf_size_bytes" \
    "png_size_bytes" \
    "status" \
    > "${QC_TABLE}"

tail -n +2 "${PLOT_MANIFEST}" |
while IFS=$'\t' read -r \
    PLOT_NUMBER \
    COMPARISON_ID \
    PLOT_GROUP \
    QUERY \
    SUBJECT \
    ANCHOR_FILE \
    QUERY_BED \
    SUBJECT_BED \
    ANCHOR_PAIR_COUNT \
    TITLE
do
    echo
    echo "============================================================"
    echo "Plotting ${COMPARISON_ID}"
    echo "============================================================"

    PDF_FILE="${PDF_DIR}/${COMPARISON_ID}.synteny_dotplot.pdf"
    PNG_FILE="${PNG_DIR}/${COMPARISON_ID}.synteny_dotplot.png"

    PDF_LOG="${LOG_DIR}/${COMPARISON_ID}.pdf.log"
    PNG_LOG="${LOG_DIR}/${COMPARISON_ID}.png.log"

    rm -f \
        "${PDF_FILE}" \
        "${PNG_FILE}" \
        "${PDF_LOG}" \
        "${PNG_LOG}"

    COMMON_OPTIONS=(
        "${ANCHOR_FILE}"
        "--qbed=${QUERY_BED}"
        "--sbed=${SUBJECT_BED}"
        "--figsize=12x12"
        "--dpi=300"
        "--nmax=100000"
        "--minfont=5"
        "--colororientation"
        "--title=${TITLE}"
    )

    set +e

    python -m jcvi.graphics.dotplot \
        "${COMMON_OPTIONS[@]}" \
        "${EXTRA_OPTIONS[@]}" \
        --format=pdf \
        --outfile="${PDF_FILE}" \
        > "${PDF_LOG}" \
        2>&1

    PDF_EXIT=$?

    python -m jcvi.graphics.dotplot \
        "${COMMON_OPTIONS[@]}" \
        "${EXTRA_OPTIONS[@]}" \
        --format=png \
        --outfile="${PNG_FILE}" \
        > "${PNG_LOG}" \
        2>&1

    PNG_EXIT=$?

    set -e

    if (
        [[ "${PDF_EXIT}" -eq 0 ]] &&
        [[ "${PNG_EXIT}" -eq 0 ]] &&
        [[ -s "${PDF_FILE}" ]] &&
        [[ -s "${PNG_FILE}" ]]
    ); then
        STATUS="PASS"
    else
        STATUS="FAIL"

        echo "ERROR: Plot failed for ${COMPARISON_ID}" >&2
        echo "PDF exit: ${PDF_EXIT}" >&2
        echo "PNG exit: ${PNG_EXIT}" >&2
        echo "PDF log:" >&2
        cat "${PDF_LOG}" >&2
        echo "PNG log:" >&2
        cat "${PNG_LOG}" >&2
    fi

    PDF_SIZE=$(
        stat -c '%s' "${PDF_FILE}" \
        2>/dev/null || echo 0
    )

    PNG_SIZE=$(
        stat -c '%s' "${PNG_FILE}" \
        2>/dev/null || echo 0
    )

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${COMPARISON_ID}" \
        "${PLOT_GROUP}" \
        "${ANCHOR_PAIR_COUNT}" \
        "${PDF_FILE}" \
        "${PNG_FILE}" \
        "${PDF_SIZE}" \
        "${PNG_SIZE}" \
        "${STATUS}" \
        >> "${QC_TABLE}"

    if [[ "${STATUS}" != "PASS" ]]; then
        exit 1
    fi
done

# ============================================================
# Final QC
# ============================================================

EXPECTED_PLOTS=$(
    awk -F'\t' '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
    ' "${PLOT_MANIFEST}"
)

PASS_PLOTS=$(
    awk -F'\t' '
        NR > 1 && $8 == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${QC_TABLE}"
)

PDF_COUNT=$(
    find "${PDF_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.pdf' \
        -size +0c \
        | wc -l
)

PNG_COUNT=$(
    find "${PNG_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.png' \
        -size +0c \
        | wc -l
)

if [[ "${EXPECTED_PLOTS}" -ne 17 ]]; then
    echo "ERROR: Expected 17 selected plots; found ${EXPECTED_PLOTS}." >&2
    exit 1
fi

if [[ "${PASS_PLOTS}" -ne "${EXPECTED_PLOTS}" ]]; then
    echo "ERROR: Only ${PASS_PLOTS}/${EXPECTED_PLOTS} plots passed." >&2
    exit 1
fi

if [[ "${PDF_COUNT}" -ne "${EXPECTED_PLOTS}" ]]; then
    echo "ERROR: Expected ${EXPECTED_PLOTS} PDFs; found ${PDF_COUNT}." >&2
    exit 1
fi

if [[ "${PNG_COUNT}" -ne "${EXPECTED_PLOTS}" ]]; then
    echo "ERROR: Expected ${EXPECTED_PLOTS} PNGs; found ${PNG_COUNT}." >&2
    exit 1
fi

# ============================================================
# Summary and checkpoint
# ============================================================

SUMMARY_TABLE="${TABLE_DIR}/dotplot_summary.tsv"

cat > "${SUMMARY_TABLE}" <<EOF2
metricvalue
plots_selected${EXPECTED_PLOTS}
plots_pass${PASS_PLOTS}
pdf_files${PDF_COUNT}
png_files${PNG_COUNT}
statusPASS
EOF2

echo
echo "============================================================"
echo "Dotplot summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY_TABLE}"

echo
echo "============================================================"
echo "Dotplot QC"
echo "============================================================"

column -t -s $'\t' \
    "${QC_TABLE}"

cp -f \
    "${SUMMARY_TABLE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${QC_TABLE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${PLOT_MANIFEST}" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/JCVI_SYNTENY_DOTPLOTS_COMPLETE.txt" <<EOF2
checkpoint=jcvi_synteny_dotplots
date=$(date --iso-8601=seconds)
plots_selected=${EXPECTED_PLOTS}
plots_pass=${PASS_PLOTS}
pdf_files=${PDF_COUNT}
png_files=${PNG_COUNT}
status=PASS
next_step=inspect_priority_plots_and_generate_karyotype_ribbon_figures
EOF2

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "============================================================"
echo "JCVI synteny dotplots completed"
echo "============================================================"
echo "PDF directory:"
echo "${PDF_DIR}"
echo
echo "PNG directory:"
echo "${PNG_DIR}"
