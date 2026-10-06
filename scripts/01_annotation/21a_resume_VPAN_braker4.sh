#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4-0:00:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=braker4_VPAN_resume
#SBATCH --array=0
#SBATCH --output=07_annotation/logs/braker4_VPAN_resume_%A_%a.out
#SBATCH --error=07_annotation/logs/braker4_VPAN_resume_%A_%a.err

set -euo pipefail

# ============================================================
# Fixed project paths
#
# Use <CONTAINER_WORK_PREFIX> rather than <HOST_WORK_PREFIX> because <CONTAINER_WORK_PREFIX>
# is explicitly mounted inside the Singularity containers.
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SINGULARITY_ARGS="${SINGULARITY_ARGS:--B ${PROJECT_DIR}}" 
cd "${PROJECT_DIR}"

CODE="VPAN"
SPECIES="Veronica_panormitana"
MODE="ETP"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

BRAKER4_DIR="${BRAKER4_DIR:?Set BRAKER4_DIR to the BRAKER4 installation directory}"
SNAKEFILE="${BRAKER4_DIR}/Snakefile"

SNAKEMAKE_ENV="${SNAKEMAKE_ENV:?Set SNAKEMAKE_ENV to the Snakemake environment activate script}"

GENOME="${PROJECT_DIR}/07_annotation/inputs/genomes/VPAN.nuclear.fa"
PROTEINS="${PROTEINS:?Set PROTEINS to the Viridiplantae protein FASTA}"

VPAN_R1="${PROJECT_DIR}/06_rnaseq/clean/VPAN/SRR3491905_1.clean.fastq.gz"
VPAN_R2="${PROJECT_DIR}/06_rnaseq/clean/VPAN/SRR3491905_2.clean.fastq.gz"

RUN_DIR="${PROJECT_DIR}/07_annotation/braker4_runs/VPAN"
SAMPLES_FILE="${RUN_DIR}/samples.csv"
CONFIG_FILE="${RUN_DIR}/config.ini"
AUGUSTUS_CONFIG="${RUN_DIR}/augustus_config"

CACHE_DIR="${PROJECT_DIR}/07_annotation/.singularity_cache"
TASK_TMP="${PROJECT_DIR}/07_annotation/tmp/VPAN"
LOG_DIR="${PROJECT_DIR}/07_annotation/logs"

STATUS_FILE="${RUN_DIR}/run_status.tsv"
VERSION_FILE="${RUN_DIR}/software_versions.txt"

mkdir -p \
    "${RUN_DIR}" \
    "${AUGUSTUS_CONFIG}" \
    "${CACHE_DIR}" \
    "${TASK_TMP}" \
    "${LOG_DIR}"

# ============================================================
# Load the same software environment used successfully before
# ============================================================

module purge

module load hpc-env/13.1
module load Python/3.11.3-GCCcore-13.1.0
module load SAMtools/1.18-GCC-13.1.0

if [[ ! -s "${SNAKEMAKE_ENV}" ]]; then
    echo "ERROR: Snakemake activation file is missing:" >&2
    echo "${SNAKEMAKE_ENV}" >&2
    exit 1
fi

source "${SNAKEMAKE_ENV}"

# Prevent host Python settings from leaking into containers.
unset PYTHONPATH
unset PYTHONHOME

export TMPDIR="${TASK_TMP}"
export SINGULARITY_TMPDIR="${TASK_TMP}"
export SINGULARITY_CACHEDIR="${CACHE_DIR}"

# ============================================================
# Check required programs
# ============================================================

for TOOL in python snakemake singularity samtools gzip
do
    if ! command -v "${TOOL}" >/dev/null 2>&1; then
        echo "ERROR: Required executable is unavailable: ${TOOL}" >&2
        exit 1
    fi

    echo "${TOOL}: $(command -v "${TOOL}")"
done

# ============================================================
# Check required input files
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

check_file "${SNAKEFILE}" "BRAKER4 Snakefile"
check_file "${GENOME}" "VPAN nuclear genome"
check_file "${PROTEINS}" "Viridiplantae protein database"
check_file "${VPAN_R1}" "VPAN RNA-seq R1"
check_file "${VPAN_R2}" "VPAN RNA-seq R2"

gzip -t "${VPAN_R1}"
gzip -t "${VPAN_R2}"

if ! grep -q '^>' "${GENOME}"; then
    echo "ERROR: No FASTA headers found in ${GENOME}" >&2
    exit 1
fi

if ! grep -q '^>' "${PROTEINS}"; then
    echo "ERROR: No FASTA headers found in ${PROTEINS}" >&2
    exit 1
fi

if grep '^>' "${GENOME}" | grep -q '[[:space:]]'; then
    echo "ERROR: Genome FASTA headers contain whitespace." >&2
    exit 1
fi


# ============================================================
# Check that inputs are visible inside the BRAKER container
#
# Reuse the Singularity cache populated during the first run.
# ============================================================

BRAKER_IMAGE="${CACHE_DIR}/dbdad69110452498a7d3d6879e51d31a.simg"

if [[ -s "${BRAKER_IMAGE}" ]]; then
    echo "Testing input visibility with cached BRAKER container..."

    singularity exec \
        ${SINGULARITY_ARGS} \
        "${BRAKER_IMAGE}" \
        bash -c "
            set -euo pipefail

            test -s '${GENOME}'
            test -s '${PROTEINS}'
            test -s '${VPAN_R1}'
            test -s '${VPAN_R2}'

            gzip -t '${VPAN_R1}'
            gzip -t '${VPAN_R2}'
        "

    echo "PASS: All inputs are visible inside the container."
else
    echo "WARNING: Cached BRAKER container was not found at:"
    echo "${BRAKER_IMAGE}"
    echo "Snakemake will resolve the configured container image."
fi

# ============================================================
# Record software versions and provenance
# ============================================================

{
    echo "date=$(date --iso-8601=seconds)"
    echo "hostname=$(hostname)"
    echo "slurm_job_id=${SLURM_JOB_ID:-NA}"
    echo "slurm_array_job_id=${SLURM_ARRAY_JOB_ID:-NA}"
    echo "slurm_array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
    echo "species_code=${CODE}"
    echo "species_name=${SPECIES}"
    echo "evidence_mode=${MODE}"
    echo "project_directory=${PROJECT_DIR}"
    echo "genome=${GENOME}"
    echo "proteins=${PROTEINS}"
    echo "rna_r1=${VPAN_R1}"
    echo "rna_r2=${VPAN_R2}"
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
# Recreate samples.csv using canonical <CONTAINER_WORK_PREFIX> paths
# ============================================================

printf '%s\n' \
"sample_name,genome,genome_masked,protein_fasta,bam_files,fastq_r1,fastq_r2,sra_ids,varus_genus,varus_species,isoseq_bam,isoseq_fastq,busco_lineage,reference_gtf" \
> "${SAMPLES_FILE}"

printf '%s,%s,,%s,,%s,%s,,,,,,embryophyta_odb12,\n' \
    "${CODE}" \
    "${GENOME}" \
    "${PROTEINS}" \
    "${VPAN_R1}" \
    "${VPAN_R2}" \
    >> "${SAMPLES_FILE}"


SAMPLE_ROWS=$(
    awk '
        NR > 1 && NF > 0 {
            n++
        }

        END {
            print n+0
        }
    ' "${SAMPLES_FILE}"
)

if [[ "${SAMPLE_ROWS}" -ne 1 ]]; then
    echo "ERROR: Expected one sample row, found ${SAMPLE_ROWS}." >&2
    exit 1
fi

# ============================================================
# Recreate the proven BRAKER4 configuration
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
max_runtime = 120
CONFIG

export BRAKER4_CONFIG="${CONFIG_FILE}"


# ============================================================
# Print the resolved job configuration
# ============================================================

echo
echo "============================================================"
echo "BRAKER4 VPAN resume job"
echo "============================================================"
echo "SLURM job ID:     ${SLURM_JOB_ID:-NA}"
echo "Species:          ${CODE} — ${SPECIES}"
echo "Evidence mode:    ${MODE}"
echo "Genome:           ${GENOME}"
echo "Proteins:         ${PROTEINS}"
echo "RNA-seq R1:       ${VPAN_R1}"
echo "RNA-seq R2:       ${VPAN_R2}"
echo "Samples file:     ${SAMPLES_FILE}"
echo "Configuration:    ${CONFIG_FILE}"
echo "Working folder:   ${RUN_DIR}"
echo "Temporary folder: ${TASK_TMP}"
echo "Threads:          ${THREADS}"
echo "Started:          $(date)"
echo

echo "samples.csv:"
cat "${SAMPLES_FILE}"

echo
echo "config.ini:"
cat "${CONFIG_FILE}"

# ============================================================
# Record starting status
# ============================================================

START_TIME=$(date --iso-8601=seconds)

printf "job_id\tspecies_code\tspecies_name\tevidence_mode\tstart_time\tcompletion_time\tstatus\n" \
    > "${STATUS_FILE}"

printf "%s\t%s\t%s\t%s\t%s\tNA\tRUNNING\n" \
    "${SLURM_JOB_ID:-NA}" \
    "${CODE}" \
    "${SPECIES}" \
    "${MODE}" \
    "${START_TIME}" \
    >> "${STATUS_FILE}"

# ============================================================
# Remove only outputs from the failed HISAT2 alignment
#
# Do not remove:
#   - repeat-model library
#   - repeat-masked genome
#   - HISAT2 index
#   - AUGUSTUS configuration
#   - completed Snakemake rule outputs
# ============================================================

rm -f \
    "${RUN_DIR}/output/VPAN/hisat2_aligned/SRR3491905.sorted.bam" \
    "${RUN_DIR}/output/VPAN/hisat2_aligned/SRR3491905.sorted.bam.bai" \
    "${RUN_DIR}/output/VPAN/hisat2_aligned/SRR3491905.sorted.bam.tmp"*

# ============================================================
# Resume BRAKER4
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
echo "Resuming BRAKER4 workflow..."

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

# ============================================================
# Determine whether the full annotation has finished
# ============================================================

FINAL_GFF3=$(
    find "${RUN_DIR}" \
        -type f \
        \( -name "braker.gff3" -o -name "braker.gff3.gz" \) \
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

COMPLETION_TIME=$(date --iso-8601=seconds)

if [[ -n "${FINAL_GFF3}" && -n "${FINAL_PROTEINS}" ]]; then

    printf "job_id\tspecies_code\tspecies_name\tevidence_mode\tstart_time\tcompletion_time\tstatus\n" \
        > "${STATUS_FILE}"

    printf "%s\t%s\t%s\t%s\t%s\t%s\tSUCCESS\n" \
        "${SLURM_JOB_ID:-NA}" \
        "${CODE}" \
        "${SPECIES}" \
        "${MODE}" \
        "${START_TIME}" \
        "${COMPLETION_TIME}" \
        >> "${STATUS_FILE}"

    echo
    echo "============================================================"
    echo "VPAN BRAKER4 annotation completed successfully"
    echo "============================================================"
    echo "Final GFF3:     ${FINAL_GFF3}"
    echo "Final proteins: ${FINAL_PROTEINS}"
    echo "Completed:      $(date)"

else

    printf "job_id\tspecies_code\tspecies_name\tevidence_mode\tstart_time\tcompletion_time\tstatus\n" \
        > "${STATUS_FILE}"

    printf "%s\t%s\t%s\t%s\t%s\t%s\tINCOMPLETE_RESUBMIT\n" \
        "${SLURM_JOB_ID:-NA}" \
        "${CODE}" \
        "${SPECIES}" \
        "${MODE}" \
        "${START_TIME}" \
        "${COMPLETION_TIME}" \
        >> "${STATUS_FILE}"

    echo
    echo "============================================================"
    echo "BRAKER4 did not yet produce all final files"
    echo "============================================================"
    echo "The existing workflow directory has been retained."
    echo "Submit this same script again to continue."
    echo
    echo "Current candidate GFF3:     ${FINAL_GFF3:-NOT_FOUND}"
    echo "Current candidate proteins: ${FINAL_PROTEINS:-NOT_FOUND}"
fi
