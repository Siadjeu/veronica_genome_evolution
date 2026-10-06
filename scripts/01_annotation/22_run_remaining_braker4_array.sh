#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-0:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=braker4_remaining
#SBATCH --array=1-8%3
#SBATCH --output=07_annotation/logs/braker4_remaining_%A_%a.out
#SBATCH --error=07_annotation/logs/braker4_remaining_%A_%a.err

set -euo pipefail

# ============================================================
# Project configuration
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SINGULARITY_ARGS="${SINGULARITY_ARGS:--B ${PROJECT_DIR}}" 
BRAKER4_DIR="${BRAKER4_DIR:?Set BRAKER4_DIR to the BRAKER4 installation directory}"
SNAKEFILE="${BRAKER4_DIR}/Snakefile"

PROTEINS="${PROTEINS:?Set PROTEINS to the Viridiplantae protein FASTA}"

SNAKEMAKE_ACTIVATE="${SNAKEMAKE_ACTIVATE:?Set SNAKEMAKE_ACTIVATE to the Snakemake environment activate script}"

CACHE_DIR="${PROJECT_DIR}/07_annotation/.singularity_cache"
LOG_DIR="${PROJECT_DIR}/07_annotation/logs"
TMP_ROOT="${PROJECT_DIR}/07_annotation/tmp"
RUN_ROOT="${PROJECT_DIR}/07_annotation/braker4_runs"
GENOME_ROOT="${PROJECT_DIR}/07_annotation/inputs/genomes"

THREADS="${SLURM_CPUS_PER_TASK:-16}"
TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not defined}"

cd "${PROJECT_DIR}"

mkdir -p \
    "${CACHE_DIR}" \
    "${LOG_DIR}" \
    "${TMP_ROOT}" \
    "${RUN_ROOT}"

# ============================================================
# Species table
#
# VPAN is task 0 and has already been completed.
# Remaining species use protein evidence: EP.
# ============================================================

declare -A CODE
declare -A SPECIES
declare -A MODE

CODE[1]="VSCU"
SPECIES[1]="Veronica_scutellata"
MODE[1]="EP"

CODE[2]="VANA"
SPECIES[2]="Veronica_anagallis_aquatica"
MODE[2]="EP"

CODE[3]="VARV"
SPECIES[3]="Veronica_arvensis"
MODE[3]="EP"

CODE[4]="VPER"
SPECIES[4]="Veronica_persica"
MODE[4]="EP"

CODE[5]="VSER"
SPECIES[5]="Veronica_serpyllifolia"
MODE[5]="EP"

CODE[6]="VTRI"
SPECIES[6]="Veronica_triloba"
MODE[6]="EP"

CODE[7]="VVER"
SPECIES[7]="Veronica_verna"
MODE[7]="EP"

CODE[8]="PMAJ"
SPECIES[8]="Plantago_major"
MODE[8]="EP"

if [[ -z "${CODE[${TASK_ID}]:-}" ]]; then
    echo "ERROR: Invalid array task ID: ${TASK_ID}" >&2
    exit 1
fi

SP_CODE="${CODE[${TASK_ID}]}"
SP_NAME="${SPECIES[${TASK_ID}]}"
EVIDENCE_MODE="${MODE[${TASK_ID}]}"

# ============================================================
# Species-specific paths
# ============================================================

GENOME="${GENOME_ROOT}/${SP_CODE}.nuclear.fa"

RUN_DIR="${RUN_ROOT}/${SP_CODE}"
SAMPLES_FILE="${RUN_DIR}/samples.csv"
CONFIG_FILE="${RUN_DIR}/config.ini"
AUGUSTUS_CONFIG="${RUN_DIR}/augustus_config"

TASK_TMP="${TMP_ROOT}/${SP_CODE}"

STATUS_FILE="${RUN_DIR}/run_status.tsv"
VERSION_FILE="${RUN_DIR}/software_versions.txt"

RESULTS_DIR="${RUN_DIR}/output/${SP_CODE}/results"

mkdir -p \
    "${RUN_DIR}" \
    "${AUGUSTUS_CONFIG}" \
    "${TASK_TMP}" \
    "${RESULTS_DIR}"

# ============================================================
# Load HPC environment
# ============================================================

module purge
module load hpc-env/13.1
module load Python/3.11.3-GCCcore-13.1.0
module load SAMtools/1.18-GCC-13.1.0

if [[ ! -s "${SNAKEMAKE_ACTIVATE}" ]]; then
    echo "ERROR: Snakemake activation file is missing:" >&2
    echo "${SNAKEMAKE_ACTIVATE}" >&2
    exit 1
fi

source "${SNAKEMAKE_ACTIVATE}"

unset PYTHONPATH
unset PYTHONHOME

export TMPDIR="${TASK_TMP}"
export SINGULARITY_TMPDIR="${TASK_TMP}"
export SINGULARITY_CACHEDIR="${CACHE_DIR}"

# ============================================================
# Utility functions
# ============================================================

check_file() {
    local FILE="$1"
    local LABEL="$2"

    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing or empty ${LABEL}:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
}

record_status() {
    local STATUS="$1"
    local COMPLETION_TIME="$2"

    {
        printf "job_id\tspecies_code\tspecies_name\tevidence_mode\tstart_time\tcompletion_time\tstatus\n"

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "${SLURM_JOB_ID:-NA}" \
            "${SP_CODE}" \
            "${SP_NAME}" \
            "${EVIDENCE_MODE}" \
            "${START_TIME}" \
            "${COMPLETION_TIME}" \
            "${STATUS}"
    } > "${STATUS_FILE}"
}

# ============================================================
# Software checks
# ============================================================

for TOOL in python snakemake singularity samtools gzip awk grep
do
    if ! command -v "${TOOL}" >/dev/null 2>&1; then
        echo "ERROR: Required executable is unavailable: ${TOOL}" >&2
        exit 1
    fi
done

# ============================================================
# Input checks
# ============================================================

check_file "${SNAKEFILE}" "BRAKER4 Snakefile"
check_file "${GENOME}" "${SP_CODE} nuclear genome"
check_file "${PROTEINS}" "Viridiplantae protein database"

if ! grep -q '^>' "${GENOME}"; then
    echo "ERROR: No FASTA headers found in ${GENOME}" >&2
    exit 1
fi

if ! grep -q '^>' "${PROTEINS}"; then
    echo "ERROR: No FASTA headers found in ${PROTEINS}" >&2
    exit 1
fi

if grep '^>' "${GENOME}" | grep -q '[[:space:]]'; then
    echo "ERROR: Genome FASTA headers contain whitespace:" >&2
    echo "${GENOME}" >&2
    exit 1
fi


# ============================================================
# Record provenance
# ============================================================

{
    echo "date=$(date --iso-8601=seconds)"
    echo "hostname=$(hostname)"
    echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
    echo "slurm_array_job_id=${SLURM_ARRAY_JOB_ID:-NA}"
    echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
    echo "species_code=${SP_CODE}"
    echo "species_name=${SP_NAME}"
    echo "evidence_mode=${EVIDENCE_MODE}"
    echo "project_directory=${PROJECT_DIR}"
    echo "genome=${GENOME}"
    echo "proteins=${PROTEINS}"
    echo "run_directory=${RUN_DIR}"
    echo "threads=${THREADS}"
    echo
    echo "python=$(command -v python)"
    python --version
    echo
    echo "snakemake=$(command -v snakemake)"
    snakemake --version
    echo
    echo "singularity=$(command -v singularity)"
    singularity --version
    echo
    echo "samtools=$(command -v samtools)"
    samtools --version | head -n 2
} > "${VERSION_FILE}"

# ============================================================
# Create samples.csv
#
# Columns:
# sample_name
# genome
# genome_masked
# protein_fasta
# bam_files
# fastq_r1
# fastq_r2
# sra_ids
# varus_genus
# varus_species
# isoseq_bam
# isoseq_fastq
# busco_lineage
# reference_gtf
# ============================================================

python - "${SAMPLES_FILE}" "${SP_CODE}" "${GENOME}" "${PROTEINS}" <<'PY'
import csv
import sys

output_file, species_code, genome, proteins = sys.argv[1:]

header = [
    "sample_name",
    "genome",
    "genome_masked",
    "protein_fasta",
    "bam_files",
    "fastq_r1",
    "fastq_r2",
    "sra_ids",
    "varus_genus",
    "varus_species",
    "isoseq_bam",
    "isoseq_fastq",
    "busco_lineage",
    "reference_gtf",
]

row = [
    species_code,
    genome,
    "",
    proteins,
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "embryophyta_odb12",
    "",
]

with open(output_file, "w", newline="", encoding="utf-8") as handle:
    writer = csv.writer(handle)
    writer.writerow(header)
    writer.writerow(row)
PY

check_file "${SAMPLES_FILE}" "samples.csv"


COLUMN_COUNT=$(
    python - "${SAMPLES_FILE}" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="", encoding="utf-8") as handle:
    rows = list(csv.reader(handle))

if len(rows) != 2:
    raise SystemExit(
        f"ERROR: Expected two rows in samples.csv; observed {len(rows)}"
    )

if len(rows[0]) != 14:
    raise SystemExit(
        f"ERROR: Header contains {len(rows[0])} columns instead of 14"
    )

if len(rows[1]) != 14:
    raise SystemExit(
        f"ERROR: Sample row contains {len(rows[1])} columns instead of 14"
    )

print(len(rows[1]))
PY
)

if [[ "${COLUMN_COUNT}" -ne 14 ]]; then
    echo "ERROR: samples.csv did not pass column validation." >&2
    exit 1
fi

# ============================================================
# Create BRAKER4 config.ini
# ============================================================

cat > "${CONFIG_FILE}" <<CONFIG
[paths]
samples_file = ${SAMPLES_FILE}
augustus_config_path = ${AUGUSTUS_CONFIG}

[containers]
braker3_image = docker://teambraker/braker3:v3.0.10
isoseq_image = docker://teambraker/braker3:isoseq
minimap2_image = docker://katharinahoff/minimap-minisplice:v0.1
minisplice_image = docker://katharinahoff/minimap-minisplice:v0.1
red_image = docker://quay.io/biocontainers/red:2018.09.10--h9948957_3
gffcompare_image = docker://quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_1
agat_image = docker://quay.io/biocontainers/agat:1.4.1--pl5321hdfd78af_0
pybarrnap_image = docker://quay.io/biocontainers/pybarrnap:0.5.1--pyhdfd78af_0
busco_image = docker://ezlabgva/busco:v6.0.0_cv1
omark_image = docker://quay.io/biocontainers/omark:0.4.1--pyh7e72e81_0
tetools_image = docker://dfam/tetools:latest
varus_image = docker://katharinahoff/varus-notebook:v0.0.6

[PARAMS]
fungus = 0
min_contig = 10000
use_varus = 0
skip_optimize_augustus = 0
skip_single_exon_downsampling = 0
downsampling_lambda = 2
use_dev_shm = 0
use_compleasm_hints = 1
skip_busco = 0
run_omark = 0
translation_table = 1
gc_donor = 0.001
allow_hinted_splicesites = gcag,atac
augustus_chunksize = 3000000
augustus_overlap = 500000
run_ncrna = 0
run_best_by_compleasm = 1
masking_tool = repeatmasker
use_minisplice = 0
no_cleanup = 0

[fantasia]
enable = 0

[OMARK]

[SLURM_ARGS]
cpus_per_task = 16
mem_of_node = 120000
max_runtime = 4320
CONFIG

check_file "${CONFIG_FILE}" "BRAKER4 config.ini"


export BRAKER4_CONFIG="${CONFIG_FILE}"

# ============================================================
# Print resolved configuration
# ============================================================

echo
echo "============================================================"
echo "BRAKER4 annotation"
echo "============================================================"
echo "SLURM job:       ${SLURM_JOB_ID:-NA}"
echo "Array task:      ${TASK_ID}"
echo "Species code:    ${SP_CODE}"
echo "Species name:    ${SP_NAME}"
echo "Evidence mode:   ${EVIDENCE_MODE}"
echo "Genome:          ${GENOME}"
echo "Protein DB:      ${PROTEINS}"
echo "Run directory:   ${RUN_DIR}"
echo "Samples file:    ${SAMPLES_FILE}"
echo "Config file:     ${CONFIG_FILE}"
echo "Temporary dir:   ${TASK_TMP}"
echo "Threads:         ${THREADS}"
echo "Start time:      $(date)"
echo

echo "samples.csv:"
cat "${SAMPLES_FILE}"

echo
echo "config.ini:"
cat "${CONFIG_FILE}"

START_TIME=$(date --iso-8601=seconds)
record_status "RUNNING" "NA"

# ============================================================
# If complete results already exist, do not rerun the species
# ============================================================

EXISTING_GFF3=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.gff3" -o -name "braker.gff3.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

EXISTING_AA=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.aa" -o -name "braker.aa.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

EXISTING_LONGEST_AA=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.longest.aa" -o -name "braker.longest.aa.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

if [[ -n "${EXISTING_GFF3}" &&
      -n "${EXISTING_AA}" &&
      -n "${EXISTING_LONGEST_AA}" ]]
then
    COMPLETION_TIME=$(date --iso-8601=seconds)
    record_status "ALREADY_COMPLETE" "${COMPLETION_TIME}"

    echo
    echo "Annotation outputs already exist for ${SP_CODE}."
    echo "GFF3:            ${EXISTING_GFF3}"
    echo "Proteins:        ${EXISTING_AA}"
    echo "Longest proteins:${EXISTING_LONGEST_AA}"
    exit 0
fi

# ============================================================
# Resume or start Snakemake
# ============================================================

cd "${RUN_DIR}"

echo
echo "Removing a possible stale Snakemake lock..."

python -m snakemake \
    -s "${SNAKEFILE}" \
    --unlock \
    --cores 1 \
    2>/dev/null || true

echo
echo "Starting or resuming BRAKER4 for ${SP_CODE}..."

set +e

python -m snakemake \
    -s "${SNAKEFILE}" \
    --cores "${THREADS}" \
    --use-singularity \
    --singularity-prefix "${CACHE_DIR}" \
    --singularity-args "${SINGULARITY_ARGS}" \
    --latency-wait 120 \
    --restart-times 3 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds

SNAKEMAKE_EXIT=$?

set -e

# ============================================================
# Validate final compressed or uncompressed outputs
# ============================================================

FINAL_GFF3=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.gff3" -o -name "braker.gff3.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

FINAL_GTF=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.gtf" -o -name "braker.gtf.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

FINAL_AA=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.aa" -o -name "braker.aa.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

FINAL_CDS=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.codingseq" -o -name "braker.codingseq.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

FINAL_LONGEST_AA=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.longest.aa" -o -name "braker.longest.aa.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

FINAL_LONGEST_GTF=$(
    find "${RESULTS_DIR}" \
        -maxdepth 1 \
        -type f \
        \( -name "braker.longest.gtf" -o -name "braker.longest.gtf.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

COMPLETION_TIME=$(date --iso-8601=seconds)

if [[ -n "${FINAL_GFF3}" &&
      -n "${FINAL_GTF}" &&
      -n "${FINAL_AA}" &&
      -n "${FINAL_CDS}" &&
      -n "${FINAL_LONGEST_AA}" &&
      -n "${FINAL_LONGEST_GTF}" ]]
then
    # Verify compressed files where applicable.
    for FILE in \
        "${FINAL_GFF3}" \
        "${FINAL_GTF}" \
        "${FINAL_AA}" \
        "${FINAL_CDS}" \
        "${FINAL_LONGEST_AA}" \
        "${FINAL_LONGEST_GTF}"
    do
        if [[ "${FILE}" == *.gz ]]; then
            gzip -t "${FILE}"
        fi
    done

    record_status "SUCCESS" "${COMPLETION_TIME}"

    echo
    echo "============================================================"
    echo "${SP_CODE} annotation completed successfully"
    echo "============================================================"
    echo "GFF3:               ${FINAL_GFF3}"
    echo "GTF:                ${FINAL_GTF}"
    echo "All proteins:       ${FINAL_AA}"
    echo "Coding sequences:   ${FINAL_CDS}"
    echo "Longest proteins:   ${FINAL_LONGEST_AA}"
    echo "Longest GTF:        ${FINAL_LONGEST_GTF}"
    echo "Snakemake exit:     ${SNAKEMAKE_EXIT}"
    echo "Completion time:    ${COMPLETION_TIME}"
    echo

    exit 0
fi

record_status "INCOMPLETE_OR_FAILED" "${COMPLETION_TIME}"

echo
echo "============================================================"
echo "${SP_CODE} annotation is incomplete"
echo "============================================================"
echo "Snakemake exit code: ${SNAKEMAKE_EXIT}"
echo "GFF3:               ${FINAL_GFF3:-NOT_FOUND}"
echo "GTF:                ${FINAL_GTF:-NOT_FOUND}"
echo "All proteins:       ${FINAL_AA:-NOT_FOUND}"
echo "Coding sequences:   ${FINAL_CDS:-NOT_FOUND}"
echo "Longest proteins:   ${FINAL_LONGEST_AA:-NOT_FOUND}"
echo "Longest GTF:        ${FINAL_LONGEST_GTF:-NOT_FOUND}"
echo
echo "The run directory has been retained for safe resumption:"
echo "${RUN_DIR}"

exit "${SNAKEMAKE_EXIT}"
