#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=VPAN_SPROT
#SBATCH --output=logs/VPAN_SPROT_%j.out
#SBATCH --error=logs/VPAN_SPROT_%j.err

set -euo pipefail

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "${PROJECT_ROOT}"

mkdir -p \
    logs \
    results/annotation/functional

QUERY="${QUERY:-/dss/work/mixu4440/Veronica_genome_evolution/08_comparative_inputs/proteomes/VPAN.longest.aa.fa}"

DB="${DB:-databases/swissprot/uniprot_sprot}"

OUT="results/annotation/functional/VPAN_swissprot.diamond.tsv"


###############################################################################
# Validate inputs
###############################################################################

if [[ ! -s "${QUERY}" ]]
then
    echo "ERROR: VPAN representative protein FASTA missing:"
    echo "${QUERY}"
    exit 1
fi

if [[ ! -s "${DB}.dmnd" ]]
then
    echo "ERROR: Swiss-Prot DIAMOND database missing:"
    echo "${DB}.dmnd"
    exit 1
fi

if ! command -v diamond >/dev/null 2>&1
then
    echo "ERROR: diamond not found in PATH."
    exit 1
fi


###############################################################################
# QC
###############################################################################

echo "============================================================"
echo "VPAN Swiss-Prot functional annotation"
echo "============================================================"

echo "DIAMOND:"
diamond version

echo
echo "Query:"
echo "${QUERY}"

echo
echo "Number of representative proteins:"
grep -c '^>' "${QUERY}"

echo
echo "Database:"
echo "${DB}.dmnd"

echo
echo "Started:"
date --iso-8601=seconds


###############################################################################
# DIAMOND BLASTP
###############################################################################

diamond blastp \
    --query "${QUERY}" \
    --db "${DB}" \
    --out "${OUT}" \
    --outfmt 6 \
        qseqid \
        sseqid \
        pident \
        length \
        qlen \
        slen \
        evalue \
        bitscore \
        qcovhsp \
        scovhsp \
        stitle \
    --evalue 1e-5 \
    --max-target-seqs 1 \
    --sensitive \
    --threads "${SLURM_CPUS_PER_TASK}"


###############################################################################
# QC output
###############################################################################

echo
echo "============================================================"
echo "DIAMOND results"
echo "============================================================"

echo "Raw best-hit rows:"
wc -l "${OUT}"

echo
echo "Unique proteins with a Swiss-Prot hit:"
cut -f1 "${OUT}" \
  | sort -u \
  | wc -l

echo
echo "Output:"
ls -lh "${OUT}"

echo
echo "Finished:"
date --iso-8601=seconds
