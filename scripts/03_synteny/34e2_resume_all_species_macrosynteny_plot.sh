#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=macro_all9
#SBATCH --output=10_synteny/logs/macrosynteny_all9_resume_%j.out
#SBATCH --error=10_synteny/logs/macrosynteny_all9_resume_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

OUTDIR="${SYNTENY_DIR}/jcvi_step34/macrosynteny_all_species"
BEDDIR="${OUTDIR}/numeric_beds"
ANCHORDIR="${OUTDIR}/anchors"
TABLEDIR="${OUTDIR}/tables"
CONFIGDIR="${OUTDIR}/config"
PLOTDIR="${OUTDIR}/plots"
LOGDIR="${OUTDIR}/logs"

CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_macrosynteny_all_species"

SEQIDS_FILE="${CONFIGDIR}/all_species.numeric.seqids"
LAYOUT_FILE="${CONFIGDIR}/all_species.numeric.layout"

PDF_FILE="${PLOTDIR}/all_9_species.macrosynteny.numeric_chromosomes.pdf"
PNG_FILE="${PLOTDIR}/all_9_species.macrosynteny.numeric_chromosomes.png"
SVG_FILE="${PLOTDIR}/all_9_species.macrosynteny.numeric_chromosomes.svg"

mkdir -p \
    "${PLOTDIR}" \
    "${LOGDIR}" \
    "${CHECKPOINT_DIR}" \
    "${SYNTENY_DIR}/logs"

cd "${PROJECT_DIR}"

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

echo "============================================================"
echo "Software environment"
echo "============================================================"

python --version

python - <<'PY'
import ete4
import jcvi
import jcvi.graphics.karyotype

print("ETE4:", getattr(ete4, "__version__", "unknown"))
print("JCVI:", getattr(jcvi, "__version__", "unknown"))
print("JCVI karyotype import: PASS")
PY

# ============================================================
# Validate all prepared inputs
# ============================================================

for FILE in \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    "${TABLEDIR}/numeric_chromosome_mapping.tsv" \
    "${TABLEDIR}/numeric_bed_qc.tsv" \
    "${TABLEDIR}/macrosynteny_species_order.tsv" \
    "${TABLEDIR}/macrosynteny_edges.tsv"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required prepared file is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

for CODE in \
    PMAJ \
    VPAN \
    VSCU \
    VANA \
    VARV \
    VVER \
    VPER \
    VSER \
    VTRI
do
    BED_FILE="${BEDDIR}/${CODE}.numeric.bed"

    if [[ ! -s "${BED_FILE}" ]]; then
        echo "ERROR: Numeric BED missing or empty:" >&2
        echo "${BED_FILE}" >&2
        exit 1
    fi

    NON_NUMERIC=$(
        awk -F'\t' '
            $1 !~ /^[0-9]+$/ {
                count++
            }
            END {
                print count + 0
            }
        ' "${BED_FILE}"
    )

    if [[ "${NON_NUMERIC}" -ne 0 ]]; then
        echo "ERROR: Non-numeric chromosome IDs remain in:" >&2
        echo "${BED_FILE}" >&2
        exit 1
    fi
done

EXPECTED_EDGES=8

OBSERVED_EDGES=$(
    awk -F'\t' '
        NR > 1 && $6 == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${TABLEDIR}/macrosynteny_edges.tsv"
)

if [[ "${OBSERVED_EDGES}" -ne "${EXPECTED_EDGES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_EDGES} prepared edges;" >&2
    echo "observed ${OBSERVED_EDGES}." >&2
    exit 1
fi

while IFS=$'\t' read -r \
    UPPER \
    LOWER \
    COMPARISON \
    SOURCE_ANCHOR \
    PLOT_ANCHOR \
    STATUS
do
    if [[ "${UPPER}" == "upper_species" ]]; then
        continue
    fi

    if [[ ! -s "${PLOT_ANCHOR}" ]]; then
        echo "ERROR: Prepared simple-anchor file missing:" >&2
        echo "${PLOT_ANCHOR}" >&2
        exit 1
    fi
done < "${TABLEDIR}/macrosynteny_edges.tsv"

# ============================================================
# Detect supported JCVI options
# ============================================================

HELP_FILE="${LOGDIR}/jcvi_karyotype_help.txt"

python -m jcvi.graphics.karyotype --help \
    > "${HELP_FILE}" \
    2>&1

EXTRA_OPTIONS=()

if grep -q -- '--notex' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--notex")
fi

if grep -q -- '--shadestyle' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--shadestyle=curve")
fi

echo
echo "============================================================"
echo "Plotting options"
echo "============================================================"
printf ' %q' "${EXTRA_OPTIONS[@]}"
echo

# ============================================================
# Remove only failed plot products
# ============================================================

rm -f \
    "${PDF_FILE}" \
    "${PNG_FILE}" \
    "${SVG_FILE}" \
    "${LOGDIR}/macrosynteny_pdf.log" \
    "${LOGDIR}/macrosynteny_png.log" \
    "${LOGDIR}/macrosynteny_svg.log"

# ============================================================
# PDF plot
# ============================================================

set +e

python -m jcvi.graphics.karyotype \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    --format=pdf \
    --outfile="${PDF_FILE}" \
    --figsize=16x14 \
    --dpi=300 \
    "${EXTRA_OPTIONS[@]}" \
    > "${LOGDIR}/macrosynteny_pdf.log" \
    2>&1

PDF_EXIT=$?

set -e

if [[ "${PDF_EXIT}" -ne 0 || ! -s "${PDF_FILE}" ]]; then
    echo "ERROR: PDF macrosynteny plot failed." >&2
    echo "Exit code: ${PDF_EXIT}" >&2
    cat "${LOGDIR}/macrosynteny_pdf.log" >&2
    exit 1
fi

# ============================================================
# PNG plot
# ============================================================

set +e

python -m jcvi.graphics.karyotype \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    --format=png \
    --outfile="${PNG_FILE}" \
    --figsize=16x14 \
    --dpi=400 \
    "${EXTRA_OPTIONS[@]}" \
    > "${LOGDIR}/macrosynteny_png.log" \
    2>&1

PNG_EXIT=$?

set -e

if [[ "${PNG_EXIT}" -ne 0 || ! -s "${PNG_FILE}" ]]; then
    echo "ERROR: PNG macrosynteny plot failed." >&2
    echo "Exit code: ${PNG_EXIT}" >&2
    cat "${LOGDIR}/macrosynteny_png.log" >&2
    exit 1
fi

# ============================================================
# Optional SVG plot
# ============================================================

set +e

python -m jcvi.graphics.karyotype \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    --format=svg \
    --outfile="${SVG_FILE}" \
    --figsize=16x14 \
    "${EXTRA_OPTIONS[@]}" \
    > "${LOGDIR}/macrosynteny_svg.log" \
    2>&1

SVG_EXIT=$?

set -e

if [[ "${SVG_EXIT}" -eq 0 && -s "${SVG_FILE}" ]]; then
    SVG_STATUS="PASS"
else
    SVG_STATUS="OPTIONAL_SVG_FAILED"
    rm -f "${SVG_FILE}"
fi

# ============================================================
# Dynamic chromosome-count validation
# ============================================================

MAPPING_ROWS=$(
    awk -F'\t' '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
    ' "${TABLEDIR}/numeric_chromosome_mapping.tsv"
)

EXPECTED_MAPPING_ROWS=$(
    awk -F'\t' '
        NR > 1 {
            total += $2
        }
        END {
            print total + 0
        }
    ' "${TABLEDIR}/numeric_bed_qc.tsv"
)

if [[ "${MAPPING_ROWS}" -ne "${EXPECTED_MAPPING_ROWS}" ]]; then
    echo "ERROR: Chromosome mapping count mismatch." >&2
    echo "Mapping rows: ${MAPPING_ROWS}" >&2
    echo "Expected from BED QC: ${EXPECTED_MAPPING_ROWS}" >&2
    exit 1
fi

PDF_SIZE=$(
    stat -c '%s' "${PDF_FILE}"
)

PNG_SIZE=$(
    stat -c '%s' "${PNG_FILE}"
)

if [[ -s "${SVG_FILE}" ]]; then
    SVG_SIZE=$(
        stat -c '%s' "${SVG_FILE}"
    )
else
    SVG_SIZE=0
fi

SUMMARY="${TABLEDIR}/macrosynteny_plot_summary.tsv"

cat > "${SUMMARY}" <<EOF2
metricvalue
species_tracks9
adjacent_synteny_edges8
numeric_chromosome_mappings${MAPPING_ROWS}
pdf_file${PDF_FILE}
pdf_size_bytes${PDF_SIZE}
png_file${PNG_FILE}
png_size_bytes${PNG_SIZE}
svg_file${SVG_FILE}
svg_size_bytes${SVG_SIZE}
svg_status${SVG_STATUS}
statusPASS
EOF2

echo
echo "============================================================"
echo "All-species macrosynteny summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

# ============================================================
# Checkpoint
# ============================================================

rm -f "${CHECKPOINT_DIR}"/*

cp -f \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/numeric_chromosome_mapping.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/numeric_bed_qc.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/macrosynteny_species_order.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/macrosynteny_edges.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${SEQIDS_FILE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${LAYOUT_FILE}" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/ALL_SPECIES_MACROSYNTENY_COMPLETE.txt" <<EOF2
checkpoint=all_species_macrosynteny_numeric_chromosomes
date=$(date --iso-8601=seconds)
species_tracks=9
adjacent_synteny_edges=8
numeric_chromosome_mappings=${MAPPING_ROWS}
pdf_status=PASS
png_status=PASS
svg_status=${SVG_STATUS}
status=PASS
next_step=inspect_and_refine_macrosynteny_figure
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
echo "All-species macrosynteny plot completed"
echo "============================================================"
echo "PDF:"
echo "${PDF_FILE}"
echo
echo "PNG:"
echo "${PNG_FILE}"

if [[ -s "${SVG_FILE}" ]]; then
    echo
    echo "SVG:"
    echo "${SVG_FILE}"
fi
