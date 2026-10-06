#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=VSCU_VSER
#SBATCH --output=10_synteny/synteny_optimized_pairwise/additional_prerequisite_comparisons/logs/VSCU_VSER_%j.out
#SBATCH --error=10_synteny/synteny_optimized_pairwise/additional_prerequisite_comparisons/logs/VSCU_VSER_%j.err

set -euo pipefail

###############################################################################
# Step 35Y1
#
# Generate the one missing JCVI full-anchor comparison required for the
# additional synteny-optimized pairwise macrosynteny figures:
#
#     VSCU -> VSER
#
# This deliberately reproduces the successful Step 35Y0 PMAJ-VSCU method:
#
#   database type    = protein
#   alignment        = diamond_blastp
#   cscore           = 0.70
#   JCVI ortholog    = jcvi.compara.catalog ortholog
#
# Existing VSCU-VPAN and VPAN-VPER Step-34 results are NOT modified.
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

SYNTENY_DIR="10_synteny"

INPUT_BED_DIR="${SYNTENY_DIR}/jcvi_inputs/bed"
INPUT_PEP_DIR="${SYNTENY_DIR}/jcvi_inputs/pep"

OUTPUT_ROOT="${SYNTENY_DIR}/synteny_optimized_pairwise/additional_prerequisite_comparisons"

WORK_DIR="${OUTPUT_ROOT}/VSCU_VSER"
TABLE_DIR="${OUTPUT_ROOT}/tables"
LOG_DIR="${OUTPUT_ROOT}/logs"

CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_step35Y1"

RESULT_FILE="${TABLE_DIR}/step35Y1_additional_prerequisite_results.tsv"
SUMMARY_FILE="${TABLE_DIR}/VSCU_VSER_summary.tsv"

mkdir -p \
    "${WORK_DIR}" \
    "${TABLE_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

###############################################################################
# Activate JCVI environment
###############################################################################

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
echo "Step 35Y1: VSCU-VSER protein-based JCVI comparison"
echo "============================================================"

echo
echo "Date:"
date --iso-8601=seconds

echo
echo "Python:"
python --version

echo
echo "DIAMOND:"
diamond version

echo
echo "JCVI:"
python - <<'PY'
import jcvi

print(
    getattr(
        jcvi,
        "__version__",
        "unknown",
    )
)
PY

###############################################################################
# Required input validation
###############################################################################

REQUIRED_INPUTS=(
    "${INPUT_BED_DIR}/VSCU.bed"
    "${INPUT_BED_DIR}/VSER.bed"
    "${INPUT_PEP_DIR}/VSCU.pep"
    "${INPUT_PEP_DIR}/VSER.pep"
)

for FILE in "${REQUIRED_INPUTS[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required JCVI input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# BED / protein count QC
###############################################################################

VSCU_BED_GENES="$(
    awk '
        NF >= 4 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${INPUT_BED_DIR}/VSCU.bed"
)"

VSER_BED_GENES="$(
    awk '
        NF >= 4 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${INPUT_BED_DIR}/VSER.bed"
)"

VSCU_PEP_GENES="$(
    grep -c '^>' \
        "${INPUT_PEP_DIR}/VSCU.pep"
)"

VSER_PEP_GENES="$(
    grep -c '^>' \
        "${INPUT_PEP_DIR}/VSER.pep"
)"

if [[ "${VSCU_BED_GENES}" -ne "${VSCU_PEP_GENES}" ]]; then
    echo "ERROR: VSCU BED/protein counts differ." >&2
    echo "BED genes: ${VSCU_BED_GENES}" >&2
    echo "Protein IDs: ${VSCU_PEP_GENES}" >&2
    exit 1
fi

if [[ "${VSER_BED_GENES}" -ne "${VSER_PEP_GENES}" ]]; then
    echo "ERROR: VSER BED/protein counts differ." >&2
    echo "BED genes: ${VSER_BED_GENES}" >&2
    echo "Protein IDs: ${VSER_PEP_GENES}" >&2
    exit 1
fi

echo
echo "============================================================"
echo "Input counts"
echo "============================================================"

echo "VSCU BED genes: ${VSCU_BED_GENES}"
echo "VSCU proteins:  ${VSCU_PEP_GENES}"
echo "VSER BED genes: ${VSER_BED_GENES}"
echo "VSER proteins:  ${VSER_PEP_GENES}"

###############################################################################
# Prepare isolated comparison directory
###############################################################################

cd "${PROJECT_ROOT}/${WORK_DIR}"

#
# Targeted cleanup only.
#
# Nothing from the original Step 34 or Step 35Y workflow is removed.
#

rm -f \
    VSCU.bed \
    VSCU.pep \
    VSER.bed \
    VSER.pep \
    VSCU.VSER.last \
    VSCU.VSER.last.filtered \
    VSCU.VSER.anchors \
    VSCU.VSER.anchors.new \
    VSCU.VSER.anchors.simple \
    VSCU.VSER.lifted.anchors \
    VSCU.VSER.pdf \
    VSCU.VSER.png \
    VSCU.VSER.svg \
    VSCU.dmnd \
    VSER.dmnd \
    VSCU.VSER.blast \
    VSCU.VSER.blast.filtered

rm -f \
    "${PROJECT_ROOT}/${RESULT_FILE}" \
    "${PROJECT_ROOT}/${SUMMARY_FILE}"

rm -f \
    "${PROJECT_ROOT}/${CHECKPOINT_DIR}/STEP35Y1_COMPLETE.txt" \
    "${PROJECT_ROOT}/${CHECKPOINT_DIR}/step35Y1_additional_prerequisite_results.tsv" \
    "${PROJECT_ROOT}/${CHECKPOINT_DIR}/VSCU_VSER_summary.tsv" \
    "${PROJECT_ROOT}/${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# Link validated JCVI inputs
###############################################################################

ln -s \
    "${PROJECT_ROOT}/${INPUT_BED_DIR}/VSCU.bed" \
    VSCU.bed

ln -s \
    "${PROJECT_ROOT}/${INPUT_PEP_DIR}/VSCU.pep" \
    VSCU.pep

ln -s \
    "${PROJECT_ROOT}/${INPUT_BED_DIR}/VSER.bed" \
    VSER.bed

ln -s \
    "${PROJECT_ROOT}/${INPUT_PEP_DIR}/VSER.pep" \
    VSER.pep

for FILE in \
    VSCU.bed \
    VSCU.pep \
    VSER.bed \
    VSER.pep
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Working-directory link is missing or broken:" >&2
        echo "${WORK_DIR}/${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# Run the same protein-based JCVI ortholog comparison used in Step 35Y0
###############################################################################

echo
echo "============================================================"
echo "Running JCVI VSCU-VSER protein ortholog comparison"
echo "============================================================"

python -m jcvi.compara.catalog ortholog \
    VSCU \
    VSER \
    --dbtype=prot \
    --align_soft=diamond_blastp \
    --cscore=.70 \
    --no_strip_names \
    --cpus="${SLURM_CPUS_PER_TASK}"

###############################################################################
# Locate full primary anchor file
###############################################################################

ANCHOR_FILE="${PROJECT_ROOT}/${WORK_DIR}/VSCU.VSER.anchors"

if [[ ! -s "${ANCHOR_FILE}" ]]; then
    echo "ERROR: Expected VSCU.VSER.anchors was not generated." >&2

    echo
    echo "Files present in work directory:" >&2

    find "${PROJECT_ROOT}/${WORK_DIR}" \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s bytes\n' \
        | sort \
        >&2

    exit 1
fi

###############################################################################
# Anchor counts
###############################################################################

FULL_ANCHOR_PAIRS="$(
    awk '
        NF >= 2 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${ANCHOR_FILE}"
)"

ANCHOR_BLOCKS="$(
    grep -c '^###' \
        "${ANCHOR_FILE}" \
        || true
)"

if [[ "${FULL_ANCHOR_PAIRS}" -le 0 ]]; then
    echo "ERROR: VSCU-VSER anchor file has no anchor pairs." >&2
    exit 1
fi

if [[ "${ANCHOR_BLOCKS}" -le 0 ]]; then
    echo "ERROR: VSCU-VSER anchor file has no syntenic blocks." >&2
    exit 1
fi

echo
echo "Full anchor pairs: ${FULL_ANCHOR_PAIRS}"
echo "Anchor blocks:     ${ANCHOR_BLOCKS}"

###############################################################################
# Verify anchor orientation and mapping back to BED IDs
###############################################################################

python - \
    "${ANCHOR_FILE}" \
    "${PROJECT_ROOT}/${INPUT_BED_DIR}/VSCU.bed" \
    "${PROJECT_ROOT}/${INPUT_BED_DIR}/VSER.bed" <<'PY'

from __future__ import annotations

import sys
from pathlib import Path


anchor_file = Path(sys.argv[1])
vscu_bed = Path(sys.argv[2])
vser_bed = Path(sys.argv[3])


def read_bed_ids(path: Path):
    ids = set()

    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) >= 4:
                ids.add(fields[3])

    return ids


vscu_ids = read_bed_ids(vscu_bed)
vser_ids = read_bed_ids(vser_bed)

total = 0
mapped = 0

unmapped_vscu = set()
unmapped_vser = set()

with anchor_file.open(
    "r",
    encoding="utf-8",
    errors="replace",
) as handle:

    for line in handle:

        stripped = line.strip()

        if (
            not stripped
            or stripped.startswith("#")
        ):
            continue

        fields = stripped.split()

        if len(fields) < 2:
            continue

        vscu_gene = fields[0]
        vser_gene = fields[1]

        total += 1

        vscu_present = vscu_gene in vscu_ids
        vser_present = vser_gene in vser_ids

        if not vscu_present:
            unmapped_vscu.add(vscu_gene)

        if not vser_present:
            unmapped_vser.add(vser_gene)

        if (
            vscu_present
            and vser_present
        ):
            mapped += 1


if total == 0:
    raise SystemExit(
        "ERROR: No anchor pairs available for mapping QC."
    )


mapping_fraction = mapped / total

print(f"Anchor pairs total:  {total}")
print(f"Anchor pairs mapped: {mapped}")
print(f"Mapping fraction:    {mapping_fraction:.6f}")

print(
    "Unmapped VSCU identifiers: "
    f"{len(unmapped_vscu)}"
)

print(
    "Unmapped VSER identifiers: "
    f"{len(unmapped_vser)}"
)

if mapping_fraction < 0.80:
    raise SystemExit(
        "ERROR: VSCU-VSER anchor mapping fraction "
        f"is below 0.80: {mapping_fraction:.6f}"
    )

PY

###############################################################################
# Write isolated prerequisite result table
###############################################################################

cd "${PROJECT_ROOT}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "comparison_id" \
    "query_species" \
    "subject_species" \
    "comparison_type" \
    "cscore" \
    "minspan" \
    "align_soft" \
    "anchor_file" \
    "anchor_pair_count" \
    "status" \
    > "${RESULT_FILE}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "VSCU__VSER" \
    "VSCU" \
    "VSER" \
    "step35Y1_additional_prerequisite" \
    "0.70" \
    "5" \
    "diamond_blastp" \
    "${WORK_DIR}/VSCU.VSER.anchors" \
    "${FULL_ANCHOR_PAIRS}" \
    "PASS" \
    >> "${RESULT_FILE}"

###############################################################################
# Write detailed summary
###############################################################################

printf '%s\t%s\n' \
    "metric" \
    "value" \
    > "${SUMMARY_FILE}"

{
    printf 'comparison\tVSCU__VSER\n'
    printf 'query_species\tVSCU\n'
    printf 'subject_species\tVSER\n'

    printf 'database_type\tprotein\n'
    printf 'alignment_backend\tdiamond_blastp\n'
    printf 'cscore\t0.70\n'
    printf 'minspan\t5\n'

    printf 'VSCU_bed_genes\t%s\n' "${VSCU_BED_GENES}"
    printf 'VSCU_proteins\t%s\n' "${VSCU_PEP_GENES}"

    printf 'VSER_bed_genes\t%s\n' "${VSER_BED_GENES}"
    printf 'VSER_proteins\t%s\n' "${VSER_PEP_GENES}"

    printf 'full_anchor_pairs\t%s\n' "${FULL_ANCHOR_PAIRS}"
    printf 'anchor_blocks\t%s\n' "${ANCHOR_BLOCKS}"

    printf 'anchor_file\t%s\n' \
        "${WORK_DIR}/VSCU.VSER.anchors"

    printf 'status\tPASS\n'

} >> "${SUMMARY_FILE}"

###############################################################################
# Validate TSV structure
###############################################################################

python - \
    "${RESULT_FILE}" \
    "${SUMMARY_FILE}" <<'PY'

from __future__ import annotations

import csv
import sys
from pathlib import Path


result_file = Path(sys.argv[1])
summary_file = Path(sys.argv[2])


with result_file.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:

    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )


if len(rows) != 1:
    raise SystemExit(
        f"ERROR: Expected one result row; found {len(rows)}."
    )


row = rows[0]

required = {
    "comparison_id",
    "query_species",
    "subject_species",
    "anchor_file",
    "anchor_pair_count",
    "status",
}

missing = required.difference(row)

if missing:
    raise SystemExit(
        "ERROR: Result table missing columns: "
        + ",".join(sorted(missing))
    )


if row["comparison_id"] != "VSCU__VSER":
    raise SystemExit(
        "ERROR: Unexpected comparison ID."
    )


if row["status"] != "PASS":
    raise SystemExit(
        "ERROR: Result status is not PASS."
    )


if int(row["anchor_pair_count"]) <= 0:
    raise SystemExit(
        "ERROR: anchor_pair_count is nonpositive."
    )


if (
    not summary_file.is_file()
    or summary_file.stat().st_size == 0
):
    raise SystemExit(
        "ERROR: Summary file missing or empty."
    )


print(
    "PASS: Step 35Y1 prerequisite tables validated."
)

PY

###############################################################################
# Create checkpoint
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP35Y1_COMPLETE.txt" <<EOF2
checkpoint=step35Y1_VSCU_VSER_comparison
date=$(date --iso-8601=seconds)
comparison=VSCU__VSER
database_type=protein
alignment_backend=diamond_blastp
cscore=0.70
minspan=5
full_anchor_pairs=${FULL_ANCHOR_PAIRS}
anchor_blocks=${ANCHOR_BLOCKS}
anchor_file=${WORK_DIR}/VSCU.VSER.anchors
existing_VSCU_VPAN_reused=true
existing_VPAN_VPER_reused=true
existing_step34_outputs_modified=false
existing_step35Y_outputs_modified=false
status=PASS
next_step=additional_synteny_optimized_pairwise_plots
EOF2

cp -f \
    "${RESULT_FILE}" \
    "${SUMMARY_FILE}" \
    "${CHECKPOINT_DIR}/"

###############################################################################
# Checksum freeze
###############################################################################

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# Display results
###############################################################################

echo
echo "============================================================"
echo "VSCU-VSER comparison summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY_FILE}"

echo
echo "============================================================"
echo "Step 35Y1 result table"
echo "============================================================"

column -t -s $'\t' \
    "${RESULT_FILE}"

echo
echo "============================================================"
echo "Step 35Y1 checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP35Y1_COMPLETE.txt"

echo
echo "============================================================"
echo "Step 35Y1 completed successfully"
echo "============================================================"
