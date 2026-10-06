#!/bin/bash

set -euo pipefail

MANIFEST="${1:-config/genome_accessions.tsv}"
OUTDIR="${2:-data/genomes}"

mkdir -p "${OUTDIR}"

if ! command -v datasets >/dev/null 2>&1
then
    echo "ERROR: NCBI datasets CLI not found."
    exit 1
fi

tail -n +2 "${MANIFEST}" |
while IFS=$'\t' read -r code species accession ploidy chromosomes source
do
    if [[ "${accession}" == "LOCAL_ASSEMBLY" ]]
    then
        echo "Skipping ${code}: local assembly"
        continue
    fi

    echo "============================================================"
    echo "${code}  ${species}  ${accession}"
    echo "============================================================"

    ZIP="${OUTDIR}/${code}.${accession}.zip"
    TMP="${OUTDIR}/${code}_tmp"

    rm -rf "${TMP}"
    mkdir -p "${TMP}"

    datasets download genome accession "${accession}" \
        --include genome \
        --filename "${ZIP}"

    unzip -q \
        "${ZIP}" \
        -d "${TMP}"

    FASTA="$(
        find "${TMP}" \
            -type f \
            \( -name '*.fna' -o -name '*.fa' -o -name '*.fasta' \) \
            | head -1
    )"

    if [[ -z "${FASTA}" ]]
    then
        echo "ERROR: genome FASTA not found for ${code}"
        exit 1
    fi

    cp \
        "${FASTA}" \
        "${OUTDIR}/${code}.genome.fa"

    echo "Saved:"
    echo "${OUTDIR}/${code}.genome.fa"

    rm -rf "${TMP}" "${ZIP}"
done
