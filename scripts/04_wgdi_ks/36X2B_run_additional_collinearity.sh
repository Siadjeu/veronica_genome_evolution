#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --array=0-5%6
#SBATCH --job-name=wgdi36X2B
#SBATCH --output=11_wgdi/logs/step36X2B_%A_%a.out
#SBATCH --error=11_wgdi/logs/step36X2B_%A_%a.err

set -euo pipefail

###############################################################################
# PROJECT CONFIGURATION
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

MANIFEST="${WGDI_ROOT}/00_admin/step36X2_collinearity_manifest.tsv"

GFF_DIR="${WGDI_ROOT}/01_inputs/gff"
LENS_DIR="${WGDI_ROOT}/01_inputs/lens"
HOMOLOGY_DIR="${WGDI_ROOT}/04_homology/results"

CONFIG_ROOT="${WGDI_ROOT}/05_collinearity/configs"
RESULT_ROOT="${WGDI_ROOT}/05_collinearity/results"
QC_ROOT="${WGDI_ROOT}/02_qc/collinearity"
LOG_DIR="${WGDI_ROOT}/logs"

mkdir -p \
    "${CONFIG_ROOT}" \
    "${RESULT_ROOT}" \
    "${QC_ROOT}" \
    "${LOG_DIR}"

###############################################################################
# ENVIRONMENT
###############################################################################

cd "${PROJECT_ROOT}"

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

for PROGRAM in wgdi python gzip
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: Required program unavailable: ${PROGRAM}" >&2
        exit 1
    fi
done

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Step 36C manifest missing:" >&2
    echo "${MANIFEST}" >&2
    exit 1
fi

###############################################################################
# READ ARRAY TASK
###############################################################################

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
    echo "ERROR: No manifest entry for task ${SLURM_ARRAY_TASK_ID}." >&2
    exit 1
fi

IFS=$'\t' read -r \
    TASK_ID \
    COMPARISON_TYPE \
    SPECIES1 \
    SPECIES2 \
    COMPARISON \
    PARAMETER_SET \
    MG \
    <<< "${TASK_LINE}"

RUN_ID="${COMPARISON}_${PARAMETER_SET}"

###############################################################################
# INPUTS
###############################################################################

GFF1="${GFF_DIR}/${SPECIES1}.wgdi.gff"
GFF2="${GFF_DIR}/${SPECIES2}.wgdi.gff"

LENS1="${LENS_DIR}/${SPECIES1}.wgdi.lens"
LENS2="${LENS_DIR}/${SPECIES2}.wgdi.lens"

BLAST_GZ="${HOMOLOGY_DIR}/${COMPARISON}.blast.tsv.gz"

for FILE in \
    "${GFF1}" \
    "${GFF2}" \
    "${LENS1}" \
    "${LENS2}" \
    "${BLAST_GZ}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required input missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# OUTPUTS
###############################################################################

CONFIG_DIR="${CONFIG_ROOT}/${RUN_ID}"
RESULT_DIR="${RESULT_ROOT}/${RUN_ID}"
QC_DIR="${QC_ROOT}/${RUN_ID}"

mkdir -p \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${QC_DIR}"

CONFIG_FILE="${CONFIG_DIR}/${RUN_ID}.icl.conf"

BLAST_FILE="${RESULT_DIR}/${COMPARISON}.blast.tsv"
COLLINEARITY_FILE="${RESULT_DIR}/${RUN_ID}.collinearity"

STDOUT_CAPTURE="${RESULT_DIR}/${RUN_ID}.wgdi.stdout.txt"
STDERR_CAPTURE="${RESULT_DIR}/${RUN_ID}.wgdi.stderr.txt"

QC_FILE="${QC_DIR}/${RUN_ID}.collinearity_qc.tsv"

###############################################################################
# ALWAYS REMOVE LARGE TEMPORARY BLAST FILE
###############################################################################

cleanup()
{
    rm -f "${BLAST_FILE}"
}

trap cleanup EXIT

rm -f "${QC_FILE}"

###############################################################################
# CONFIGURATION
###############################################################################

cat > "${CONFIG_FILE}" <<EOF2
[collinearity]
gff1 = ${GFF1}
gff2 = ${GFF2}
lens1 = ${LENS1}
lens2 = ${LENS2}
blast = ${BLAST_FILE}
blast_reverse = false
comparison = genomes
multiple = 1
process = ${SLURM_CPUS_PER_TASK:-8}
evalue = 1e-5
score = 100
grading = 50,30,25
mg = ${MG}
pvalue = 1
repeat_number = 20
positon = order
savefile = ${COLLINEARITY_FILE}
EOF2

###############################################################################
# RUN WGDI ONLY WHEN A USABLE RESULT IS NOT ALREADY PRESENT
###############################################################################

echo "============================================================"
echo "Step 36X2B: additional WGDI improved collinearity"
echo "============================================================"
echo "Task: ${TASK_ID}"
echo "Comparison: ${COMPARISON}"
echo "Type: ${COMPARISON_TYPE}"
echo "Parameter set: ${PARAMETER_SET}"
echo "mg: ${MG}"
echo "Start: $(date --iso-8601=seconds)"
echo "WGDI executable: $(command -v wgdi)"

REUSE_EXISTING="NO"

if [[ -s "${COLLINEARITY_FILE}" ]]; then
    ALIGNMENT_COUNT="$(
        grep -c '^# Alignment ' "${COLLINEARITY_FILE}" || true
    )"

    DATA_COUNT="$(
        awk '
            NF > 0 && $1 !~ /^#/ {
                count++
            }
            END {
                print count + 0
            }
        ' "${COLLINEARITY_FILE}"
    )"

    if [[ "${ALIGNMENT_COUNT}" -gt 0 && "${DATA_COUNT}" -gt 0 ]]; then
        REUSE_EXISTING="YES"
    fi
fi

if [[ "${REUSE_EXISTING}" == "YES" ]]; then
    echo "Existing nonempty WGDI collinearity output detected."
    echo "WGDI computation will be reused; QC will be rerun."

else
    echo "No reusable result detected. Running WGDI."

    rm -f \
        "${COLLINEARITY_FILE}" \
        "${STDOUT_CAPTURE}" \
        "${STDERR_CAPTURE}"

    gzip -cd \
        "${BLAST_GZ}" \
        > "${BLAST_FILE}"

    if [[ ! -s "${BLAST_FILE}" ]]; then
        echo "ERROR: Failed to decompress homology input." >&2
        exit 1
    fi

    set +e

    wgdi -icl \
        "${CONFIG_FILE}" \
        > "${STDOUT_CAPTURE}" \
        2> "${STDERR_CAPTURE}"

    WGDI_EXIT_CODE=$?

    set -e

    if [[ "${WGDI_EXIT_CODE}" -ne 0 ]]; then
        echo "ERROR: WGDI exited with code ${WGDI_EXIT_CODE}." >&2

        echo "WGDI stderr:" >&2
        tail -n 100 "${STDERR_CAPTURE}" >&2 || true

        exit "${WGDI_EXIT_CODE}"
    fi
fi

if [[ ! -s "${COLLINEARITY_FILE}" ]]; then
    echo "ERROR: WGDI produced no collinearity file." >&2

    echo "WGDI stdout:" >&2
    tail -n 100 "${STDOUT_CAPTURE}" >&2 || true

    echo "WGDI stderr:" >&2
    tail -n 100 "${STDERR_CAPTURE}" >&2 || true

    exit 1
fi

###############################################################################
# STRICT OUTPUT QC FOR WGDI 0.75 FORMAT
#
# Header:
# # Alignment 1: score=1248 pvalue=0.1029 N=39 chr1&chr2 plus
#
# Data:
# gene1 order1 gene2 order2 score_class
###############################################################################

python - \
    "${COLLINEARITY_FILE}" \
    "${GFF1}" \
    "${GFF2}" \
    "${SPECIES1}" \
    "${SPECIES2}" \
    "${COMPARISON_TYPE}" \
    "${PARAMETER_SET}" \
    "${MG}" \
    "${REUSE_EXISTING}" \
    "${QC_FILE}" <<'PY'
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

collinearity_path = Path(sys.argv[1])
gff1_path = Path(sys.argv[2])
gff2_path = Path(sys.argv[3])

species1 = sys.argv[4]
species2 = sys.argv[5]
comparison_type = sys.argv[6]
parameter_set = sys.argv[7]
mg = sys.argv[8]
reused_existing = sys.argv[9]

qc_path = Path(sys.argv[10])


def read_gff(path: Path):
    identifiers = set()
    chromosome_by_gene = {}
    order_by_gene = {}
    chromosomes = set()

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")

            if len(fields) != 7:
                raise SystemExit(
                    f"ERROR: Invalid seven-column GFF: "
                    f"{path}:{line_number}"
                )

            chromosome = fields[0]
            gene_id = fields[1]

            try:
                order = int(fields[5])
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: Invalid gene order in "
                    f"{path}:{line_number}"
                ) from error

            if gene_id in identifiers:
                raise SystemExit(
                    f"ERROR: Duplicate GFF gene ID: {gene_id}"
                )

            identifiers.add(gene_id)
            chromosomes.add(chromosome)
            chromosome_by_gene[gene_id] = chromosome
            order_by_gene[gene_id] = order

    if not identifiers:
        raise SystemExit(f"ERROR: Empty GFF: {path}")

    return (
        identifiers,
        chromosomes,
        chromosome_by_gene,
        order_by_gene,
    )


(
    ids1,
    chromosomes1,
    chromosome1_by_gene,
    order1_by_gene,
) = read_gff(gff1_path)

(
    ids2,
    chromosomes2,
    chromosome2_by_gene,
    order2_by_gene,
) = read_gff(gff2_path)

header_pattern = re.compile(
    r"^# Alignment\s+"
    r"(?P<number>\d+):\s+"
    r"score=(?P<score>[^\s]+)\s+"
    r"pvalue=(?P<pvalue>[^\s]+)\s+"
    r"N=(?P<n>\d+)\s+"
    r"(?P<chromosome_pair>[^\s]+)\s+"
    r"(?P<orientation>[^\s]+)"
)

block_count = 0
gene_pair_count = 0

invalid_gene_pairs = 0
invalid_order_values = 0
order_mismatches = 0
malformed_data_lines = 0
malformed_header_lines = 0

same_gene_pairs = 0
duplicate_gene_pairs = 0

genes1_in_blocks = set()
genes2_in_blocks = set()

gene_pair_counter = Counter()
block_sizes = []
reported_block_sizes = []
block_scores = []
block_pvalues = []
orientations = Counter()

current_observed_size = None
current_reported_size = None

comment_lines = 0
data_lines = 0

with collinearity_path.open(
    "r",
    encoding="utf-8",
    errors="replace",
) as handle:
    for line_number, raw_line in enumerate(handle, start=1):
        line = raw_line.strip()

        if not line:
            continue

        if line.startswith("#"):
            comment_lines += 1

            if line.startswith("# Alignment"):
                if current_observed_size is not None:
                    block_sizes.append(current_observed_size)
                    reported_block_sizes.append(
                        current_reported_size
                    )

                match = header_pattern.match(line)

                if match is None:
                    malformed_header_lines += 1
                    current_observed_size = None
                    current_reported_size = None
                    continue

                block_count += 1

                current_observed_size = 0
                current_reported_size = int(
                    match.group("n")
                )

                try:
                    block_scores.append(
                        float(match.group("score"))
                    )
                    block_pvalues.append(
                        float(match.group("pvalue"))
                    )
                except ValueError:
                    malformed_header_lines += 1

                orientations[
                    match.group("orientation")
                ] += 1

            continue

        data_lines += 1
        fields = line.split()

        # WGDI 0.75 expected structure:
        # gene1 order1 gene2 order2 score_class
        if len(fields) < 5:
            malformed_data_lines += 1
            continue

        gene1 = fields[0]
        gene2 = fields[2]

        try:
            reported_order1 = int(fields[1])
            reported_order2 = int(fields[3])
            float(fields[4])
        except ValueError:
            malformed_data_lines += 1
            continue

        if current_observed_size is None:
            malformed_data_lines += 1
            continue

        if gene1 not in ids1 or gene2 not in ids2:
            invalid_gene_pairs += 1
            continue

        expected_order1 = order1_by_gene[gene1]
        expected_order2 = order2_by_gene[gene2]

        if (
            reported_order1 != expected_order1
            or reported_order2 != expected_order2
        ):
            order_mismatches += 1
            continue

        gene_pair_count += 1
        current_observed_size += 1

        genes1_in_blocks.add(gene1)
        genes2_in_blocks.add(gene2)

        gene_pair_counter[(gene1, gene2)] += 1

        if gene1 == gene2:
            same_gene_pairs += 1

if current_observed_size is not None:
    block_sizes.append(current_observed_size)
    reported_block_sizes.append(current_reported_size)

duplicate_gene_pairs = sum(
    count - 1
    for count in gene_pair_counter.values()
    if count > 1
)

if block_count == 0:
    raise SystemExit(
        "ERROR: No WGDI '# Alignment' blocks were detected."
    )

if gene_pair_count == 0:
    raise SystemExit(
        "ERROR: No valid WGDI five-column gene pairs were parsed."
    )

if malformed_header_lines != 0:
    raise SystemExit(
        f"ERROR: {malformed_header_lines} malformed alignment "
        "headers were detected."
    )

if malformed_data_lines != 0:
    raise SystemExit(
        f"ERROR: {malformed_data_lines} malformed data lines "
        "were detected."
    )

if invalid_gene_pairs != 0:
    raise SystemExit(
        f"ERROR: {invalid_gene_pairs} data rows contain "
        "invalid gene IDs."
    )

if invalid_order_values != 0:
    raise SystemExit(
        f"ERROR: {invalid_order_values} invalid order values "
        "were detected."
    )

if order_mismatches != 0:
    raise SystemExit(
        f"ERROR: {order_mismatches} WGDI/GFF order mismatches "
        "were detected."
    )

if len(block_sizes) != block_count:
    raise SystemExit(
        f"ERROR: Parsed {len(block_sizes)} block sizes but "
        f"found {block_count} alignment headers."
    )

block_size_mismatches = sum(
    observed != reported
    for observed, reported in zip(
        block_sizes,
        reported_block_sizes,
    )
)

if block_size_mismatches != 0:
    raise SystemExit(
        f"ERROR: {block_size_mismatches} blocks differ between "
        "reported N and observed gene-pair count."
    )

minimum_block_size = min(block_sizes)
maximum_block_size = max(block_sizes)
mean_block_size = sum(block_sizes) / len(block_sizes)
median_block_size = sorted(block_sizes)[len(block_sizes) // 2]

mean_block_score = (
    sum(block_scores) / len(block_scores)
    if block_scores
    else math.nan
)

mean_block_pvalue = (
    sum(block_pvalues) / len(block_pvalues)
    if block_pvalues
    else math.nan
)

row = {
    "comparison": f"{species1}_{species2}",
    "comparison_type": comparison_type,
    "parameter_set": parameter_set,
    "mg": mg,
    "reused_existing_wgdi_output": reused_existing,
    "species1_gene_count": len(ids1),
    "species2_gene_count": len(ids2),
    "species1_chromosome_count": len(chromosomes1),
    "species2_chromosome_count": len(chromosomes2),
    "block_count": block_count,
    "collinear_gene_pairs": gene_pair_count,
    "unique_collinear_gene_pairs": len(gene_pair_counter),
    "duplicate_gene_pairs": duplicate_gene_pairs,
    "species1_genes_in_blocks": len(genes1_in_blocks),
    "species2_genes_in_blocks": len(genes2_in_blocks),
    "species1_block_gene_coverage": (
        f"{len(genes1_in_blocks) / len(ids1):.8f}"
    ),
    "species2_block_gene_coverage": (
        f"{len(genes2_in_blocks) / len(ids2):.8f}"
    ),
    "minimum_block_size": minimum_block_size,
    "median_block_size": median_block_size,
    "maximum_block_size": maximum_block_size,
    "mean_block_size": f"{mean_block_size:.8f}",
    "mean_block_score": f"{mean_block_score:.8f}",
    "mean_block_pvalue": f"{mean_block_pvalue:.8f}",
    "plus_blocks": orientations.get("plus", 0),
    "minus_blocks": orientations.get("minus", 0),
    "other_orientation_blocks": sum(
        count
        for orientation, count in orientations.items()
        if orientation not in {"plus", "minus"}
    ),
    "same_gene_pairs": same_gene_pairs,
    "comment_lines": comment_lines,
    "data_lines": data_lines,
    "malformed_header_lines": malformed_header_lines,
    "malformed_data_lines": malformed_data_lines,
    "invalid_gene_pairs": invalid_gene_pairs,
    "order_mismatches": order_mismatches,
    "block_size_mismatches": block_size_mismatches,
    "status": "PASS",
}

with qc_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(row),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerow(row)

print(
    f"{row['comparison']} {parameter_set}: "
    f"blocks={block_count:,}; "
    f"pairs={gene_pair_count:,}; "
    f"coverage1={row['species1_block_gene_coverage']}; "
    f"coverage2={row['species2_block_gene_coverage']}; "
    "status=PASS"
)
PY

###############################################################################
# REPORT
###############################################################################

echo
echo "============================================================"
echo "Step 36X2B task summary"
echo "============================================================"

column -t -s $'\t' \
    "${QC_FILE}"

echo
echo "Collinearity output:"
echo "${COLLINEARITY_FILE}"

echo
echo "Existing WGDI output reused: ${REUSE_EXISTING}"
echo "Completed: $(date --iso-8601=seconds)"
echo "Step 36X2B task ${TASK_ID}: PASS"
