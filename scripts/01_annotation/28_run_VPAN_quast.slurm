#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=VPAN_QUAST
#SBATCH --output=logs/VPAN_QUAST_%j.out
#SBATCH --error=logs/VPAN_QUAST_%j.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "${PROJECT_ROOT}"

mkdir -p \
    logs \
    results/annotation/quast

ASSEMBLY="${ASSEMBLY:?Set ASSEMBLY to the VPAN chromosome-scale FASTA}"

# Optional:
# export QUAST_ENV=quast_env
QUAST_ENV="${QUAST_ENV:-}"

if [[ -n "${QUAST_ENV}" ]]; then
    QUAST_CMD=(conda run -n "${QUAST_ENV}" quast.py)
else
    QUAST_CMD=(quast.py)
fi

echo "Assembly:"
echo "  ${ASSEMBLY}"
echo

"${QUAST_CMD[@]}" \
    "${ASSEMBLY}" \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --min-contig 1 \
    --output-dir results/annotation/quast/VPAN

echo
echo "QUAST completed."

ls -lh \
    results/annotation/quast/VPAN/report.tsv \
    results/annotation/quast/VPAN/report.txt
