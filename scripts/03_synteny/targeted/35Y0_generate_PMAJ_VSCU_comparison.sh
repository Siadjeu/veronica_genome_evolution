#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=PMAJ_VSCU
#SBATCH --output=10_synteny/synteny_optimized_pairwise/prerequisite_comparisons/logs/PMAJ_VSCU_%j.out
#SBATCH --error=10_synteny/synteny_optimized_pairwise/prerequisite_comparisons/logs/PMAJ_VSCU_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

INPUT_BED_DIR="${SYNTENY_DIR}/jcvi_inputs/bed"
INPUT_PEP_DIR="${SYNTENY_DIR}/jcvi_inputs/pep"

OUTPUT_ROOT="${SYNTENY_DIR}/synteny_optimized_pairwise/prerequisite_comparisons"
WORK_DIR="${OUTPUT_ROOT}/PMAJ_VSCU"
TABLE_DIR="${OUTPUT_ROOT}/tables"
LOG_DIR="${OUTPUT_ROOT}/logs"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_step35Y0"

mkdir -p \
    "${WORK_DIR}" \
    "${TABLE_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

cd "${PROJECT_DIR}"

# ============================================================
# Activate JCVI environment
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

echo "============================================================"
echo "Step 35Y0: PMAJ-VSCU protein-based JCVI comparison"
echo "============================================================"

echo "Date:"
date --iso-8601=seconds

echo
echo "Hostname:"
hostname

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

# ============================================================
# Validate required inputs
# ============================================================

REQUIRED_INPUTS=(
    "${INPUT_BED_DIR}/PMAJ.bed"
    "${INPUT_BED_DIR}/VSCU.bed"
    "${INPUT_PEP_DIR}/PMAJ.pep"
    "${INPUT_PEP_DIR}/VSCU.pep"
)

for FILE in "${REQUIRED_INPUTS[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required JCVI input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

PMAJ_BED_GENES=$(
    awk '
        NF >= 4 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${INPUT_BED_DIR}/PMAJ.bed"
)

VSCU_BED_GENES=$(
    awk '
        NF >= 4 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${INPUT_BED_DIR}/VSCU.bed"
)

PMAJ_PEP_GENES=$(
    grep -c '^>' \
        "${INPUT_PEP_DIR}/PMAJ.pep"
)

VSCU_PEP_GENES=$(
    grep -c '^>' \
        "${INPUT_PEP_DIR}/VSCU.pep"
)

if [[ "${PMAJ_BED_GENES}" -ne "${PMAJ_PEP_GENES}" ]]; then
    echo "ERROR: PMAJ BED/protein counts differ." >&2
    echo "BED genes: ${PMAJ_BED_GENES}" >&2
    echo "Protein IDs: ${PMAJ_PEP_GENES}" >&2
    exit 1
fi

if [[ "${VSCU_BED_GENES}" -ne "${VSCU_PEP_GENES}" ]]; then
    echo "ERROR: VSCU BED/protein counts differ." >&2
    echo "BED genes: ${VSCU_BED_GENES}" >&2
    echo "Protein IDs: ${VSCU_PEP_GENES}" >&2
    exit 1
fi

echo
echo "Input counts:"
echo "PMAJ BED genes: ${PMAJ_BED_GENES}"
echo "PMAJ proteins:  ${PMAJ_PEP_GENES}"
echo "VSCU BED genes: ${VSCU_BED_GENES}"
echo "VSCU proteins:  ${VSCU_PEP_GENES}"

# ============================================================
# Prepare isolated working directory
# ============================================================

cd "${WORK_DIR}"

# Remove only files generated for this PMAJ-VSCU comparison.
rm -f \
    PMAJ.bed \
    PMAJ.pep \
    VSCU.bed \
    VSCU.pep \
    PMAJ.VSCU.last \
    PMAJ.VSCU.last.filtered \
    PMAJ.VSCU.anchors \
    PMAJ.VSCU.anchors.new \
    PMAJ.VSCU.anchors.simple \
    PMAJ.VSCU.pdf \
    PMAJ.VSCU.png \
    PMAJ.VSCU.svg \
    PMAJ.dmnd \
    VSCU.dmnd \
    PMAJ.VSCU.blast \
    PMAJ.VSCU.blast.filtered

rm -f \
    "${TABLE_DIR}/step35Y_prerequisite_results.tsv" \
    "${TABLE_DIR}/PMAJ_VSCU_summary.tsv"

rm -f \
    "${CHECKPOINT_DIR}/STEP35Y0_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/step35Y_prerequisite_results.tsv" \
    "${CHECKPOINT_DIR}/PMAJ_VSCU_summary.tsv" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

ln -s \
    "${INPUT_BED_DIR}/PMAJ.bed" \
    PMAJ.bed

ln -s \
    "${INPUT_PEP_DIR}/PMAJ.pep" \
    PMAJ.pep

ln -s \
    "${INPUT_BED_DIR}/VSCU.bed" \
    VSCU.bed

ln -s \
    "${INPUT_PEP_DIR}/VSCU.pep" \
    VSCU.pep

for FILE in \
    PMAJ.bed \
    PMAJ.pep \
    VSCU.bed \
    VSCU.pep
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Working-directory link is missing or broken:" >&2
        echo "${WORK_DIR}/${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Run protein-based JCVI comparison
#
# --dbtype=prot is essential. Without it, JCVI searches for
# PMAJ.cds and VSCU.cds instead of using the protein FASTAs.
# ============================================================

echo
echo "Running JCVI protein ortholog comparison..."

python -m jcvi.compara.catalog ortholog \
    PMAJ \
    VSCU \
    --dbtype=prot \
    --align_soft=diamond_blastp \
    --cscore=.70 \
    --no_strip_names \
    --cpus="${SLURM_CPUS_PER_TASK}"

# ============================================================
# Locate and validate the full anchor file
# ============================================================

ANCHOR_FILE="${WORK_DIR}/PMAJ.VSCU.anchors"

if [[ ! -s "${ANCHOR_FILE}" ]]; then
    echo "ERROR: Expected PMAJ.VSCU.anchors was not generated." >&2

    echo
    echo "Files present in work directory:" >&2

    find "${WORK_DIR}" \
        -maxdepth 1 \
        -type f \
        -printf '%f\t%s bytes\n' \
        | sort \
        >&2

    exit 1
fi

FULL_ANCHOR_PAIRS=$(
    awk '
        NF >= 2 && $1 !~ /^#/ {
            count++
        }
        END {
            print count + 0
        }
    ' "${ANCHOR_FILE}"
)

ANCHOR_BLOCKS=$(
    grep -c '^###' \
        "${ANCHOR_FILE}" \
        || true
)

if [[ "${FULL_ANCHOR_PAIRS}" -le 0 ]]; then
    echo "ERROR: PMAJ-VSCU anchor file has no anchor pairs." >&2
    exit 1
fi

if [[ "${ANCHOR_BLOCKS}" -le 0 ]]; then
    echo "ERROR: PMAJ-VSCU anchor file has no syntenic blocks." >&2
    exit 1
fi

# ============================================================
# Confirm anchor gene identifiers map to the BED files
# ============================================================

python - \
    "${ANCHOR_FILE}" \
    "${WORK_DIR}/PMAJ.bed" \
    "${WORK_DIR}/VSCU.bed" <<'PY'
from __future__ import annotations

import sys
from pathlib import Path

anchor_file = Path(sys.argv[1])
pmaj_bed = Path(sys.argv[2])
vscu_bed = Path(sys.argv[3])


def read_bed_ids(path):
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

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) >= 4:
                ids.add(fields[3])

    return ids


pmaj_ids = read_bed_ids(
    pmaj_bed
)

vscu_ids = read_bed_ids(
    vscu_bed
)

total = 0
mapped = 0
unmapped_pmaj = set()
unmapped_vscu = set()

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

        pmaj_gene = fields[0]
        vscu_gene = fields[1]

        total += 1

        pmaj_present = (
            pmaj_gene in pmaj_ids
        )

        vscu_present = (
            vscu_gene in vscu_ids
        )

        if not pmaj_present:
            unmapped_pmaj.add(
                pmaj_gene
            )

        if not vscu_present:
            unmapped_vscu.add(
                vscu_gene
            )

        if (
            pmaj_present
            and vscu_present
        ):
            mapped += 1

if total == 0:
    raise SystemExit(
        "ERROR: No anchor pairs available for mapping QC."
    )

mapping_fraction = (
    mapped / total
)

print(
    f"Anchor pairs total: {total}"
)

print(
    f"Anchor pairs mapped: {mapped}"
)

print(
    f"Mapping fraction: {mapping_fraction:.6f}"
)

print(
    f"Unmapped PMAJ identifiers: "
    f"{len(unmapped_pmaj)}"
)

print(
    f"Unmapped VSCU identifiers: "
    f"{len(unmapped_vscu)}"
)

if mapping_fraction < 0.80:
    raise SystemExit(
        "ERROR: PMAJ-VSCU anchor mapping fraction "
        f"is below 0.80: {mapping_fraction:.6f}"
    )
PY

# ============================================================
# Write prerequisite result table
# ============================================================

RESULT_FILE="${TABLE_DIR}/step35Y_prerequisite_results.tsv"
SUMMARY_FILE="${TABLE_DIR}/PMAJ_VSCU_summary.tsv"

cat > "${RESULT_FILE}" <<EOF2
comparison_idquery_speciessubject_speciescomparison_typeanchor_fileanchor_pair_countanchor_block_countstatus
PMAJ__VSCUPMAJVSCUstep35Y_prerequisite${ANCHOR_FILE}${FULL_ANCHOR_PAIRS}${ANCHOR_BLOCKS}PASS
EOF2

cat > "${SUMMARY_FILE}" <<EOF2
metricvalue
comparisonPMAJ__VSCU
query_speciesPMAJ
subject_speciesVSCU
database_typeprotein
alignment_backenddiamond_blastp
cscore0.70
PMAJ_bed_genes${PMAJ_BED_GENES}
PMAJ_proteins${PMAJ_PEP_GENES}
VSCU_bed_genes${VSCU_BED_GENES}
VSCU_proteins${VSCU_PEP_GENES}
full_anchor_pairs${FULL_ANCHOR_PAIRS}
anchor_blocks${ANCHOR_BLOCKS}
anchor_file${ANCHOR_FILE}
statusPASS
EOF2

# ============================================================
# Create checkpoint
# ============================================================

cat > "${CHECKPOINT_DIR}/STEP35Y0_COMPLETE.txt" <<EOF2
checkpoint=step35Y0_PMAJ_VSCU_comparison
date=$(date --iso-8601=seconds)
comparison=PMAJ__VSCU
database_type=protein
alignment_backend=diamond_blastp
cscore=0.70
full_anchor_pairs=${FULL_ANCHOR_PAIRS}
anchor_blocks=${ANCHOR_BLOCKS}
anchor_file=${ANCHOR_FILE}
status=PASS
next_step=run_step35Y_synteny_optimized_pairwise_plots
EOF2

cp -f \
    "${RESULT_FILE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${SUMMARY_FILE}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

# ============================================================
# Display final results
# ============================================================

echo
echo "============================================================"
echo "PMAJ-VSCU comparison summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY_FILE}"

echo
echo "============================================================"
echo "Step 35Y0 checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP35Y0_COMPLETE.txt"

echo
echo "Step 35Y0 completed successfully."
