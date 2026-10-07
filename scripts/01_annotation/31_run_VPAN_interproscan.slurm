#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --job-name=VPAN_IPR
#SBATCH --output=logs/VPAN_IPR_%j.out
#SBATCH --error=logs/VPAN_IPR_%j.err

set -euo pipefail


###############################################################################
# PROJECT
###############################################################################

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
cd "${PROJECT_ROOT}"

mkdir -p \
    logs \
    results/annotation/functional/interproscan/input \
    results/annotation/functional/interproscan/tmp


###############################################################################
# INPUTS
###############################################################################

ORIGINAL_QUERY="/dss/work/mixu4440/Veronica_genome_evolution/08_comparative_inputs/proteomes/VPAN.longest.aa.fa"

INTERPROSCAN="/dss/work/mixu4440/my_interproscan/interproscan-5.78-109.0/interproscan.sh"

CLEAN_QUERY="results/annotation/functional/interproscan/input/VPAN.longest.interproscan.fa"

OUT="results/annotation/functional/interproscan/VPAN.interproscan.selected_apps.tsv"

APPLICATIONS="Gene3D,PANTHER,Pfam,SMART,SUPERFAMILY"

EXPECTED_PROTEINS=28313


###############################################################################
# LOAD JAVA 11
###############################################################################

module purge
module load hpc-env/12.2
module load Java/11.0.20.1


###############################################################################
# JOB-SPECIFIC TEMP DIRECTORY
###############################################################################

TMPDIR_IPR="$(
    readlink -m \
        results/annotation/functional/interproscan/tmp/${SLURM_JOB_ID}
)"

mkdir -p "${TMPDIR_IPR}"


###############################################################################
# VALIDATE JAVA
###############################################################################

echo "============================================================"
echo "Java environment"
echo "============================================================"

which java
java -version

JAVA_VERSION="$(
    java -version 2>&1 \
    | awk -F '"' '/version/ {print $2}'
)"

JAVA_MAJOR="$(
    echo "${JAVA_VERSION}" \
    | cut -d. -f1
)"

if [[ "${JAVA_MAJOR}" != "11" ]]
then
    echo "ERROR: Java 11 is required."
    echo "Detected Java version: ${JAVA_VERSION}"
    exit 1
fi


###############################################################################
# VALIDATE INPUTS
###############################################################################

if [[ ! -s "${ORIGINAL_QUERY}" ]]
then
    echo "ERROR: original VPAN protein FASTA missing:"
    echo "${ORIGINAL_QUERY}"
    exit 1
fi


if [[ ! -x "${INTERPROSCAN}" ]]
then
    echo "ERROR: InterProScan executable missing:"
    echo "${INTERPROSCAN}"
    exit 1
fi


###############################################################################
# CLEAN FASTA FOR INTERPROSCAN
#
# InterProScan rejects "*" characters.
# Original BRAKER protein FASTA remains untouched.
###############################################################################

echo
echo "============================================================"
echo "Preparing InterProScan-compatible protein FASTA"
echo "============================================================"

awk '
    /^>/ {
        print
        next
    }
    {
        gsub(/\*/, "")
        print
    }
' "${ORIGINAL_QUERY}" > "${CLEAN_QUERY}"


###############################################################################
# FASTA QC
###############################################################################

ORIGINAL_COUNT="$(
    grep -c '^>' "${ORIGINAL_QUERY}"
)"

CLEAN_COUNT="$(
    grep -c '^>' "${CLEAN_QUERY}"
)"

ASTERISKS="$(
    grep -v '^>' "${CLEAN_QUERY}" \
    | tr -cd '*' \
    | wc -c
)"

DUPLICATES="$(
    grep '^>' "${CLEAN_QUERY}" \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | sort \
    | uniq -d \
    | wc -l
)"


echo "Original proteins:   ${ORIGINAL_COUNT}"
echo "Cleaned proteins:    ${CLEAN_COUNT}"
echo "Remaining *:         ${ASTERISKS}"
echo "Duplicate IDs:       ${DUPLICATES}"


if [[ "${ORIGINAL_COUNT}" -ne "${EXPECTED_PROTEINS}" ]]
then
    echo "ERROR: unexpected original protein count."
    echo "Expected: ${EXPECTED_PROTEINS}"
    echo "Observed: ${ORIGINAL_COUNT}"
    exit 1
fi


if [[ "${CLEAN_COUNT}" -ne "${EXPECTED_PROTEINS}" ]]
then
    echo "ERROR: cleaned FASTA changed the protein count."
    exit 1
fi


if [[ "${ASTERISKS}" -ne 0 ]]
then
    echo "ERROR: cleaned FASTA still contains * characters."
    exit 1
fi


if [[ "${DUPLICATES}" -ne 0 ]]
then
    echo "ERROR: duplicate protein IDs found."
    exit 1
fi


###############################################################################
# CHECK FOR EMPTY PROTEINS
###############################################################################

EMPTY_COUNT="$(
    awk '
        /^>/ {
            if (header != "" && length(seq) == 0) {
                empty++
            }

            header = $0
            seq = ""
            next
        }

        {
            seq = seq $0
        }

        END {
            if (header != "" && length(seq) == 0) {
                empty++
            }

            print empty + 0
        }
    ' "${CLEAN_QUERY}"
)"

echo "Empty proteins:      ${EMPTY_COUNT}"

if [[ "${EMPTY_COUNT}" -ne 0 ]]
then
    echo "ERROR: empty protein sequences detected."
    exit 1
fi


###############################################################################
# REMOVE STALE OUTPUT
###############################################################################

rm -f "${OUT}"


###############################################################################
# PROVENANCE
###############################################################################

echo
echo "============================================================"
echo "VPAN InterProScan functional annotation"
echo "============================================================"

echo "InterProScan executable:"
echo "${INTERPROSCAN}"

echo
echo "InterProScan version:"
"${INTERPROSCAN}" --version

echo
echo "Input:"
echo "${CLEAN_QUERY}"

echo
echo "Representative proteins:"
echo "${CLEAN_COUNT}"

echo
echo "Applications:"
echo "${APPLICATIONS}"

echo
echo "CPUs:"
echo "${SLURM_CPUS_PER_TASK}"

echo
echo "Temporary directory:"
echo "${TMPDIR_IPR}"

echo
echo "Output:"
echo "${OUT}"

echo
echo "Started:"
date --iso-8601=seconds


###############################################################################
# RUN INTERPROSCAN
###############################################################################

"${INTERPROSCAN}" \
    -i "${CLEAN_QUERY}" \
    -f TSV \
    -appl "${APPLICATIONS}" \
    -goterms \
    -iprlookup \
    -cpu "${SLURM_CPUS_PER_TASK}" \
    -T "${TMPDIR_IPR}" \
    -o "${OUT}"


###############################################################################
# VALIDATE OUTPUT
###############################################################################

if [[ ! -s "${OUT}" ]]
then
    echo "ERROR: InterProScan output missing or empty."
    exit 1
fi


###############################################################################
# RESULT SUMMARY
###############################################################################

ROWS="$(
    wc -l < "${OUT}"
)"

ANNOTATED="$(
    cut -f1 "${OUT}" \
    | sort -u \
    | wc -l
)"

PERCENT="$(
    awk \
        -v n="${ANNOTATED}" \
        -v total="${EXPECTED_PROTEINS}" \
        'BEGIN {
            printf "%.2f", 100*n/total
        }'
)"


echo
echo "============================================================"
echo "InterProScan results"
echo "============================================================"

echo "Output rows:"
echo "${ROWS}"

echo
echo "Proteins with >=1 InterProScan result:"
echo "${ANNOTATED}"

echo
echo "InterProScan-supported proteins:"
echo "${PERCENT}%"

echo
echo "Results by application:"

cut -f4 "${OUT}" \
    | sort \
    | uniq -c \
    | sort -nr


echo
echo "Output file:"
ls -lh "${OUT}"

echo
echo "Finished:"
date --iso-8601=seconds
