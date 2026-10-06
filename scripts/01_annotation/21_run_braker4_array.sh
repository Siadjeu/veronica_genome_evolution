#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-0:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=braker4_array
#SBATCH --array=0-8%3
#SBATCH --output=07_annotation/logs/braker4_%A_%a.out
#SBATCH --error=07_annotation/logs/braker4_%A_%a.err

set -euo pipefail

# ============================================================
# Fixed canonical project path
#
# IMPORTANT:
# Use <CONTAINER_WORK_PREFIX> rather than the <HOST_WORK_PREFIX> alias because
# <CONTAINER_WORK_PREFIX> is the path mounted inside the BRAKER4 containers.
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SINGULARITY_ARGS="${SINGULARITY_ARGS:--B ${PROJECT_DIR}}" 
cd "${PROJECT_DIR}"

ARRAY_TABLE="${PROJECT_DIR}/00_metadata/braker4_annotation_array.tsv"

BRAKER4_DIR="${BRAKER4_DIR:?Set BRAKER4_DIR to the BRAKER4 installation directory}"
SNAKEFILE="${BRAKER4_DIR}/Snakefile"

SNAKEMAKE_ENV="${SNAKEMAKE_ENV:?Set SNAKEMAKE_ENV to the Snakemake environment activate script}"

PROTEINS="${PROTEINS:?Set PROTEINS to the Viridiplantae protein FASTA}"

GENOME_DIR="${PROJECT_DIR}/07_annotation/inputs/genomes"

VPAN_R1="${PROJECT_DIR}/06_rnaseq/clean/VPAN/SRR3491905_1.clean.fastq.gz"
VPAN_R2="${PROJECT_DIR}/06_rnaseq/clean/VPAN/SRR3491905_2.clean.fastq.gz"

RUN_ROOT="${PROJECT_DIR}/07_annotation/braker4_runs"
CACHE_DIR="${PROJECT_DIR}/07_annotation/.singularity_cache"
TMP_ROOT="${PROJECT_DIR}/07_annotation/tmp"
LOG_DIR="${PROJECT_DIR}/07_annotation/logs"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p \
    "${RUN_ROOT}" \
    "${CACHE_DIR}" \
    "${TMP_ROOT}" \
    "${LOG_DIR}"

# ============================================================
# Load the exact environment used successfully previously
# ============================================================

module purge
module load hpc-env/13.1
module load Python/3.11.3-GCCcore-13.1.0
module load SAMtools/1.18-GCC-13.1.0

source "${SNAKEMAKE_ENV}"

unset PYTHONPATH
unset PYTHONHOME

# ============================================================
# Validate executables
# ============================================================

for TOOL in python snakemake singularity samtools
do
    if ! command -v "${TOOL}" >/dev/null 2>&1
    then
        echo "ERROR: Required executable unavailable: ${TOOL}" >&2
        exit 1
    fi
done

# ============================================================
# Resolve current array sample
# ============================================================

if [[ ! -s "${ARRAY_TABLE}" ]]
then
    echo "ERROR: Missing array table: ${ARRAY_TABLE}" >&2
    exit 1
fi

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 2))
SAMPLE_LINE=$(sed -n "${LINE_NUMBER}p" "${ARRAY_TABLE}")

if [[ -z "${SAMPLE_LINE}" ]]
then
    echo "ERROR: No sample for array task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

IFS=$'\t' read -r ARRAY_ID CODE SPECIES MODE <<< "${SAMPLE_LINE}"

if [[ "${ARRAY_ID}" != "${SLURM_ARRAY_TASK_ID}" ]]
then
    echo "ERROR: Array index mismatch." >&2
    echo "SLURM task: ${SLURM_ARRAY_TASK_ID}" >&2
    echo "Table index: ${ARRAY_ID}" >&2
    exit 1
fi

# ============================================================
# Task-specific paths
# ============================================================

GENOME="${GENOME_DIR}/${CODE}.nuclear.fa"

RUN_DIR="${RUN_ROOT}/${CODE}"
SAMPLES_FILE="${RUN_DIR}/samples.csv"
CONFIG_FILE="${RUN_DIR}/config.ini"
AUGUSTUS_CONFIG="${RUN_DIR}/augustus_config"
TASK_TMP="${TMP_ROOT}/${CODE}"

STATUS_FILE="${RUN_DIR}/run_status.tsv"
VERSION_FILE="${RUN_DIR}/software_versions.txt"

mkdir -p \
    "${RUN_DIR}" \
    "${AUGUSTUS_CONFIG}" \
    "${TASK_TMP}"

export TMPDIR="${TASK_TMP}"
export SINGULARITY_TMPDIR="${TASK_TMP}"
export SINGULARITY_CACHEDIR="${CACHE_DIR}"

# ============================================================
# Input checks
# ============================================================

check_file() {
    local FILE="$1"
    local LABEL="$2"

    if [[ ! -s "${FILE}" ]]
    then
        echo "ERROR: Missing or empty ${LABEL}: ${FILE}" >&2
        exit 1
    fi
}

check_file "${SNAKEFILE}" "BRAKER4 Snakefile"
check_file "${SNAKEMAKE_ENV}" "Snakemake activation file"
check_file "${GENOME}" "${CODE} nuclear genome"
check_file "${PROTEINS}" "Viridiplantae protein database"

if [[ "${CODE}" == "VPAN" ]]
then
    check_file "${VPAN_R1}" "VPAN RNA-seq R1"
    check_file "${VPAN_R2}" "VPAN RNA-seq R2"

    gzip -t "${VPAN_R1}"
    gzip -t "${VPAN_R2}"
fi

if ! grep -q '^>' "${GENOME}"
then
    echo "ERROR: No FASTA headers in genome: ${GENOME}" >&2
    exit 1
fi

if ! grep -q '^>' "${PROTEINS}"
then
    echo "ERROR: No FASTA headers in protein database: ${PROTEINS}" >&2
    exit 1
fi

if grep '^>' "${GENOME}" | grep -q '[[:space:]]'
then
    echo "ERROR: Genome headers contain whitespace: ${GENOME}" >&2
    exit 1
fi

# ============================================================
# Confirm all important paths are visible inside Singularity
# ============================================================

echo
echo "Testing Singularity bind visibility..."

singularity exec \
    ${SINGULARITY_ARGS} \
    docker://teambraker/braker3:v3.0.10 \
    bash -c "
        set -e
        test -s '${GENOME}'
        test -s '${PROTEINS}'
    "

if [[ "${CODE}" == "VPAN" ]]
then
    singularity exec \
        ${SINGULARITY_ARGS} \
        docker://teambraker/braker3:v3.0.10 \
        bash -c "
            set -e
            test -s '${VPAN_R1}'
            test -s '${VPAN_R2}'
            gzip -t '${VPAN_R1}'
            gzip -t '${VPAN_R2}'
        "
fi

echo "PASS: Inputs are visible inside the container."

# ============================================================
# Record versions
# ============================================================

{
    echo "date=$(date --iso-8601=seconds)"
    echo "hostname=$(hostname)"
    echo "species_code=${CODE}"
    echo "species_name=${SPECIES}"
    echo "evidence_mode=${MODE}"
    echo "project_dir=${PROJECT_DIR}"
    echo
    echo "python=$(which python)"
    python --version
    echo
    echo "snakemake=$(which snakemake)"
    snakemake --version
    echo
    echo "singularity=$(which singularity)"
    singularity --version
    echo
    echo "samtools=$(which samtools)"
    samtools --version | head -n 2
} > "${VERSION_FILE}"

# ============================================================
# Recreate samples.csv using only <CONTAINER_WORK_PREFIX> paths
# ============================================================

printf '%s\n' \
"sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf" \
> "${SAMPLES_FILE}"

if [[ "${CODE}" == "VPAN" ]]
then
    printf '%s,%s,,%s,,%s,%s,,,,,,embryophyta_odb12,\n' \
        "${CODE}" \
        "${GENOME}" \
        "${PROTEINS}" \
        "${VPAN_R1}" \
        "${VPAN_R2}" \
        >> "${SAMPLES_FILE}"
else
    printf '%s,%s,,%s,,,,,,,,,embryophyta_odb12,\n' \
        "${CODE}" \
        "${GENOME}" \
        "${PROTEINS}" \
        >> "${SAMPLES_FILE}"
fi

SAMPLE_ROWS=$(
    awk 'NR>1 && NF>0 {n++} END {print n+0}' \
        "${SAMPLES_FILE}"
)

if [[ "${SAMPLE_ROWS}" -ne 1 ]]
then
    echo "ERROR: Expected one sample row in ${SAMPLES_FILE}" >&2
    exit 1
fi


# ============================================================
# Recreate task-specific config.ini
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
cpus_per_task = 48
mem_of_node = 120000
max_runtime = 4320
CONFIG

export BRAKER4_CONFIG="${CONFIG_FILE}"

# Ensure config also contains no inaccessible alias.

# ============================================================
# Display resolved task
# ============================================================

echo
echo "============================================================"
echo "BRAKER4 array task"
echo "============================================================"
echo "Array job ID:      ${SLURM_ARRAY_JOB_ID}"
echo "Array task ID:     ${SLURM_ARRAY_TASK_ID}"
echo "Species code:      ${CODE}"
echo "Species name:      ${SPECIES}"
echo "Evidence mode:     ${MODE}"
echo "Genome:            ${GENOME}"
echo "Proteins:          ${PROTEINS}"
echo "Samples file:      ${SAMPLES_FILE}"
echo "Configuration:     ${CONFIG_FILE}"
echo "Working directory: ${RUN_DIR}"
echo "Threads:           ${THREADS}"
echo "Started:           $(date)"
echo

echo "samples.csv:"
cat "${SAMPLES_FILE}"

# ============================================================
# Record status
# ============================================================

START_TIME=$(date --iso-8601=seconds)

printf "array_job_id\tarray_task_id\tspecies_code\tspecies_name\tevidence_mode\tstart_time\tcompletion_time\tstatus\n" \
    > "${STATUS_FILE}"

printf "%s\t%s\t%s\t%s\t%s\t%s\tNA\tRUNNING\n" \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${CODE}" \
    "${SPECIES}" \
    "${MODE}" \
    "${START_TIME}" \
    >> "${STATUS_FILE}"

# ============================================================
# Remove only incomplete HISAT2 outputs from failed VPAN run
#
# Repeat masking and other completed work remain untouched.
# ============================================================

if [[ "${CODE}" == "VPAN" ]]
then
    rm -f \
        "${RUN_DIR}/output/VPAN/hisat2_aligned/SRR3491905.sorted.bam" \
        "${RUN_DIR}/output/VPAN/hisat2_aligned/SRR3491905.sorted.bam.bai"

    rm -f \
        "${RUN_DIR}/output/VPAN/hisat2_aligned/SRR3491905.sorted.bam.tmp"* \
        2>/dev/null || true
fi

# ============================================================
# Run or resume BRAKER4
# ============================================================

cd "${RUN_DIR}"

python -m snakemake \
    -s "${SNAKEFILE}" \
    --unlock \
    --cores 1 \
    2>/dev/null || true

python -m snakemake \
    -s "${SNAKEFILE}" \
    --cores "${THREADS}" \
    --use-singularity \
    --singularity-prefix "${CACHE_DIR}" \
    --singularity-args "${SINGULARITY_ARGS}" \
    --latency-wait 120 \
    --restart-times 3 \
    --rerun-incomplete \
    --printshellcmds

# ============================================================
# Validate final outputs
# ============================================================

FINAL_GFF3=$(
    find "${RUN_DIR}" \
        -type f \
        -name "braker.gff3" \
        -size +0c \
        -print -quit 2>/dev/null || true
)

FINAL_PROTEINS=$(
    find "${RUN_DIR}" \
        -type f \
        \( -name "braker.aa" -o -name "braker.aa.gz" \) \
        -size +0c \
        -print -quit 2>/dev/null || true
)

if [[ -z "${FINAL_GFF3}" ]]
then
    echo "ERROR: No final braker.gff3 found for ${CODE}" >&2
    exit 1
fi

if [[ -z "${FINAL_PROTEINS}" ]]
then
    echo "ERROR: No final braker.aa file found for ${CODE}" >&2
    exit 1
fi

COMPLETION_TIME=$(date --iso-8601=seconds)

printf "array_job_id\tarray_task_id\tspecies_code\tspecies_name\tevidence_mode\tstart_time\tcompletion_time\tstatus\n" \
    > "${STATUS_FILE}"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS\n" \
    "${SLURM_ARRAY_JOB_ID}" \
    "${SLURM_ARRAY_TASK_ID}" \
    "${CODE}" \
    "${SPECIES}" \
    "${MODE}" \
    "${START_TIME}" \
    "${COMPLETION_TIME}" \
    >> "${STATUS_FILE}"

echo
echo "============================================================"
echo "BRAKER4 annotation completed successfully"
echo "============================================================"
echo "Species:        ${CODE}"
echo "Final GFF3:     ${FINAL_GFF3}"
echo "Final proteins: ${FINAL_PROTEINS}"
echo "Completed:      $(date)"
