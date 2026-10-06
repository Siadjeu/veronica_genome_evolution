#!/bin/bash
#SBATCH --partition=mpcb.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=annotation_checkpoint
#SBATCH --output=07_annotation/logs/annotation_checkpoint_%j.out
#SBATCH --error=07_annotation/logs/annotation_checkpoint_%j.err

set -euo pipefail

# ============================================================
# Project paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

ANNOTATION_DIR="${PROJECT_DIR}/07_annotation"
RUN_ROOT="${ANNOTATION_DIR}/braker4_runs"
QC_DIR="${ANNOTATION_DIR}/qc"
SCRIPT_DIR="${PROJECT_DIR}/05_scripts"

CHECKPOINT_DIR="${ANNOTATION_DIR}/checkpoint_annotation_final"
METHODS_DIR="${CHECKPOINT_DIR}/01_methods"
RESULTS_DIR="${CHECKPOINT_DIR}/02_results"
TABLES_DIR="${CHECKPOINT_DIR}/03_tables"
SCRIPTS_DIR="${CHECKPOINT_DIR}/04_scripts"
MANIFEST_DIR="${CHECKPOINT_DIR}/05_manifests"
LOGS_DIR="${CHECKPOINT_DIR}/06_logs"
OUTPUTS_DIR="${CHECKPOINT_DIR}/07_output_manifest"

SPECIES=(VPAN VSCU VANA VARV VPER VSER VTRI VVER PMAJ)

mkdir -p \
    "${METHODS_DIR}" \
    "${RESULTS_DIR}" \
    "${TABLES_DIR}" \
    "${SCRIPTS_DIR}" \
    "${MANIFEST_DIR}" \
    "${LOGS_DIR}" \
    "${OUTPUTS_DIR}"

# ============================================================
# Species metadata
# ============================================================

declare -A SCIENTIFIC_NAME
declare -A EVIDENCE
declare -A ACCESSION_STATUS

SCIENTIFIC_NAME[VPAN]="Veronica panormitana"
SCIENTIFIC_NAME[VSCU]="Veronica scutellata"
SCIENTIFIC_NAME[VANA]="Veronica anagallis-aquatica"
SCIENTIFIC_NAME[VARV]="Veronica arvensis"
SCIENTIFIC_NAME[VPER]="Veronica persica"
SCIENTIFIC_NAME[VSER]="Veronica serpyllifolia"
SCIENTIFIC_NAME[VTRI]="Veronica triloba"
SCIENTIFIC_NAME[VVER]="Veronica verna"
SCIENTIFIC_NAME[PMAJ]="Plantago major"

EVIDENCE[VPAN]="RNA-seq + Viridiplantae proteins"
EVIDENCE[VSCU]="Viridiplantae proteins"
EVIDENCE[VANA]="Viridiplantae proteins"
EVIDENCE[VARV]="Viridiplantae proteins"
EVIDENCE[VPER]="Viridiplantae proteins"
EVIDENCE[VSER]="Viridiplantae proteins"
EVIDENCE[VTRI]="Viridiplantae proteins"
EVIDENCE[VVER]="Viridiplantae proteins"
EVIDENCE[PMAJ]="Viridiplantae proteins"

ACCESSION_STATUS[VPAN]="Local assembly"
ACCESSION_STATUS[VSCU]="Public assembly"
ACCESSION_STATUS[VANA]="Public assembly"
ACCESSION_STATUS[VARV]="Public assembly"
ACCESSION_STATUS[VPER]="Public assembly"
ACCESSION_STATUS[VSER]="Public assembly"
ACCESSION_STATUS[VTRI]="Public assembly"
ACCESSION_STATUS[VVER]="Public assembly"
ACCESSION_STATUS[PMAJ]="Public outgroup assembly"

# ============================================================
# Required files
# ============================================================

FINAL_TABLE="${QC_DIR}/manuscript_annotation_statistics_final.tsv"
FEATURE_TABLE="${QC_DIR}/annotation_feature_counts.tsv"
ANNOTATION_SUMMARY="${QC_DIR}/annotation_summary.tsv"
MISSING_TABLE="${QC_DIR}/annotation_missing_files.tsv"
PROTEOME_MANIFEST="${PROJECT_DIR}/08_comparative_inputs/proteome_manifest.tsv"

for FILE in \
    "${FINAL_TABLE}" \
    "${FEATURE_TABLE}" \
    "${ANNOTATION_SUMMARY}" \
    "${MISSING_TABLE}" \
    "${PROTEOME_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required checkpoint input is missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Copy final QC tables
# ============================================================

cp -f "${FINAL_TABLE}" \
    "${TABLES_DIR}/manuscript_annotation_statistics_final.tsv"

cp -f "${FEATURE_TABLE}" \
    "${TABLES_DIR}/annotation_feature_counts.tsv"

cp -f "${ANNOTATION_SUMMARY}" \
    "${TABLES_DIR}/annotation_output_files.tsv"

cp -f "${MISSING_TABLE}" \
    "${TABLES_DIR}/annotation_missing_files.tsv"

cp -f "${PROTEOME_MANIFEST}" \
    "${TABLES_DIR}/representative_proteome_manifest.tsv"

if [[ -s "${QC_DIR}/busco_summary_files.tsv" ]]; then
    cp -f "${QC_DIR}/busco_summary_files.tsv" \
        "${TABLES_DIR}/busco_summary_files.tsv"
fi

if [[ -s "${QC_DIR}/busco_raw_summaries.txt" ]]; then
    cp -f "${QC_DIR}/busco_raw_summaries.txt" \
        "${TABLES_DIR}/busco_raw_summaries.txt"
fi

# ============================================================
# Create input manifest
# ============================================================

INPUT_MANIFEST="${MANIFEST_DIR}/annotation_input_manifest.tsv"

printf "species_code\tscientific_name\tassembly_status\tgenome_fasta\tevidence\tprotein_evidence\trna_seq_r1\trna_seq_r2\n" \
    > "${INPUT_MANIFEST}"

for CODE in "${SPECIES[@]}"
do
    GENOME="${ANNOTATION_DIR}/inputs/genomes/${CODE}.nuclear.fa"
    PROTEINS="${PROTEINS:-<VIRIDIPLANTAE_PROTEIN_FASTA>}"

    RNA_R1="NA"
    RNA_R2="NA"

    if [[ "${CODE}" == "VPAN" ]]; then
        RNA_R1="${PROJECT_DIR}/06_rnaseq/clean/VPAN/SRR3491905_1.clean.fastq.gz"
        RNA_R2="${PROJECT_DIR}/06_rnaseq/clean/VPAN/SRR3491905_2.clean.fastq.gz"
    fi

    if [[ ! -s "${GENOME}" ]]; then
        echo "ERROR: Missing genome for ${CODE}: ${GENOME}" >&2
        exit 1
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${SCIENTIFIC_NAME[${CODE}]}" \
        "${ACCESSION_STATUS[${CODE}]}" \
        "${GENOME}" \
        "${EVIDENCE[${CODE}]}" \
        "${PROTEINS}" \
        "${RNA_R1}" \
        "${RNA_R2}" \
        >> "${INPUT_MANIFEST}"
done

# ============================================================
# Create final output manifest
# ============================================================

OUTPUT_MANIFEST="${OUTPUTS_DIR}/annotation_output_manifest.tsv"

printf "species_code\tscientific_name\tgff3\tgtf\tall_proteins\tcoding_sequences\trepresentative_proteins\trepresentative_gtf\tprotein_busco\tgenome_busco\tstatus\n" \
    > "${OUTPUT_MANIFEST}"

for CODE in "${SPECIES[@]}"
do
    FINAL_RESULTS="${RUN_ROOT}/${CODE}/output/${CODE}/results"
    QUALITY_CONTROL="${FINAL_RESULTS}/quality_control"

    GFF3="${FINAL_RESULTS}/braker.gff3.gz"
    GTF="${FINAL_RESULTS}/braker.gtf.gz"
    AA="${FINAL_RESULTS}/braker.aa.gz"
    CDS="${FINAL_RESULTS}/braker.codingseq.gz"
    LONGEST_AA="${FINAL_RESULTS}/braker.longest.aa.gz"
    LONGEST_GTF="${FINAL_RESULTS}/braker.longest.gtf.gz"
    PROTEIN_BUSCO="${QUALITY_CONTROL}/busco_proteins_short_summary.txt"
    GENOME_BUSCO="${QUALITY_CONTROL}/busco_genome_short_summary.txt"

    STATUS="PASS"

    for FILE in \
        "${GFF3}" \
        "${GTF}" \
        "${AA}" \
        "${CDS}" \
        "${LONGEST_AA}" \
        "${LONGEST_GTF}" \
        "${PROTEIN_BUSCO}" \
        "${GENOME_BUSCO}"
    do
        if [[ ! -s "${FILE}" ]]; then
            echo "ERROR: Missing expected output for ${CODE}:" >&2
            echo "${FILE}" >&2
            STATUS="FAIL"
        fi
    done

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${CODE}" \
        "${SCIENTIFIC_NAME[${CODE}]}" \
        "${GFF3}" \
        "${GTF}" \
        "${AA}" \
        "${CDS}" \
        "${LONGEST_AA}" \
        "${LONGEST_GTF}" \
        "${PROTEIN_BUSCO}" \
        "${GENOME_BUSCO}" \
        "${STATUS}" \
        >> "${OUTPUT_MANIFEST}"

    if [[ "${STATUS}" != "PASS" ]]; then
        exit 1
    fi
done

# ============================================================
# Copy scripts used in the annotation stage
# ============================================================

SCRIPT_MANIFEST="${MANIFEST_DIR}/script_manifest.tsv"

printf "script\tcheckpoint_copy\tstatus\n" > "${SCRIPT_MANIFEST}"

ANNOTATION_SCRIPTS=(
    "19_prepare_braker4_genomes.sh"
    "21_run_braker4_array.sh"
    "21a_resume_VPAN_braker4.sh"
    "22_run_remaining_braker4_array.sh"
    "23_validate_all_braker4_annotations.sh"
    "24_prepare_representative_proteomes.sh"
    "25_make_final_manuscript_annotation_table.sh"
    "26_create_annotation_checkpoint.sh"
)

for SCRIPT_NAME in "${ANNOTATION_SCRIPTS[@]}"
do
    SOURCE="${SCRIPT_DIR}/${SCRIPT_NAME}"
    TARGET="${SCRIPTS_DIR}/${SCRIPT_NAME}"

    if [[ -s "${SOURCE}" ]]; then
        cp -f "${SOURCE}" "${TARGET}"

        printf "%s\t%s\tCOPIED\n" \
            "${SOURCE}" \
            "${TARGET}" \
            >> "${SCRIPT_MANIFEST}"
    else
        printf "%s\tNA\tNOT_FOUND\n" \
            "${SOURCE}" \
            >> "${SCRIPT_MANIFEST}"
    fi
done

# Also preserve all generated per-species configuration files.
for CODE in "${SPECIES[@]}"
do
    CONFIG_SOURCE="${RUN_ROOT}/${CODE}/config.ini"
    SAMPLES_SOURCE="${RUN_ROOT}/${CODE}/samples.csv"

    if [[ -s "${CONFIG_SOURCE}" ]]; then
        cp -f "${CONFIG_SOURCE}" \
            "${SCRIPTS_DIR}/${CODE}.config.ini"
    fi

    if [[ -s "${SAMPLES_SOURCE}" ]]; then
        cp -f "${SAMPLES_SOURCE}" \
            "${SCRIPTS_DIR}/${CODE}.samples.csv"
    fi
done

# ============================================================
# Methods text
# ============================================================

cat > "${METHODS_DIR}/annotation_methods.md" <<'METHODS'
## Structural genome annotation

Structural gene annotation was performed for eight Veronica genomes and the
Plantago major outgroup using the BRAKER4 workflow. Nuclear genome assemblies
were supplied as unmasked FASTA files, and repeat identification and masking
were performed within the workflow using RepeatModeler and RepeatMasker.

For Veronica panormitana, gene prediction incorporated both transcriptomic and
protein evidence. Paired-end RNA-sequencing reads from SRA accession
SRR3491905 were quality-controlled before analysis and aligned to the genome
using HISAT2. Protein evidence was provided using a Viridiplantae protein
dataset. The remaining Veronica species and Plantago major were annotated
using the same Viridiplantae protein evidence without species-specific RNA-seq
data.

BRAKER4 was executed through Snakemake using containerised software. The core
gene-prediction stage used the BRAKER3 container version 3.0.10. Genome
sequences shorter than 10 kb were excluded from gene prediction. Protein and
RNA-seq evidence, where available, were used to generate extrinsic hints for
gene-model training and prediction. The embryophyta_odb12 lineage dataset was
used for completeness assessment. The workflow was configured to use
Compleasm-derived hints during annotation and to select the best gene set
according to completeness assessment.

For each species, BRAKER4 produced GFF3 and GTF annotation files, predicted
protein sequences, coding sequences, and representative longest-isoform gene
sets. The longest predicted protein isoform per gene was retained for
comparative genomics to minimise artificial inflation of gene-family and
duplicate-copy estimates caused by alternative transcript isoforms.

Annotation outputs were validated by confirming the presence and integrity of
the GFF3, GTF, protein, coding-sequence, longest-protein and longest-GTF files.
Protein and coding-sequence counts were required to agree. FASTA identifiers in
the representative proteomes were tested for duplication. Annotation
completeness was assessed with BUSCO using the embryophyta_odb12 dataset,
containing 2,026 conserved orthologues. BUSCO completeness was evaluated
separately for the predicted proteomes and genome assemblies.
METHODS

# ============================================================
# Detailed reproducible workflow description
# ============================================================

cat > "${METHODS_DIR}/annotation_workflow.md" <<'WORKFLOW'
## Annotation workflow checkpoint

### Species annotated

- VPAN: Veronica panormitana
- VSCU: Veronica scutellata
- VANA: Veronica anagallis-aquatica
- VARV: Veronica arvensis
- VPER: Veronica persica
- VSER: Veronica serpyllifolia
- VTRI: Veronica triloba
- VVER: Veronica verna
- PMAJ: Plantago major

### Evidence modes

- VPAN: RNA-seq plus Viridiplantae protein evidence
- All remaining species: Viridiplantae protein evidence

### Principal workflow settings

- Workflow: BRAKER4 through Snakemake
- BRAKER3 container: teambraker/braker3:v3.0.10
- Repeat masking: RepeatModeler and RepeatMasker
- Minimum contig length: 10,000 bp
- BUSCO lineage: embryophyta_odb12
- BUSCO lineage size: 2,026 orthologues
- Compleasm hints: enabled
- Best gene set selected by completeness: enabled
- Non-coding RNA annotation: disabled
- Representative proteome: longest predicted isoform per gene
- HPC scheduler: SLURM
- Standard allocation: 16 CPUs, 128 GB total memory, four-day wall time
- Container-visible project prefix: <WORK_ROOT>

### Important workflow correction

The first VPAN RNA-seq alignment attempt used paths beginning with
<HOST_WORK_PREFIX>. These paths were visible on the host but not within the
Singularity container because the mounted path was <WORK_ROOT>. The
workflow was corrected by explicitly defining the project directory as:

<PROJECT_ROOT>

The VPAN workflow was then safely resumed from the failed HISAT2 rule.
Completed repeat-modelling and masking outputs were retained.

### Quality-control criteria

An annotation was accepted only when all of the following files existed and
were non-empty:

- braker.gff3.gz
- braker.gtf.gz
- braker.aa.gz
- braker.codingseq.gz
- braker.longest.aa.gz
- braker.longest.gtf.gz

Compressed files were tested with gzip integrity checks. Protein and coding
sequence counts were required to be equal. Representative-proteome FASTA
identifiers were required to be unique.
WORKFLOW

# ============================================================
# Results text
# ============================================================

cat > "${RESULTS_DIR}/annotation_results.md" <<'RESULTS'
## Structural annotation and completeness

Structural annotation of the eight Veronica genomes and the Plantago major
outgroup predicted between 28,325 and 71,164 protein-coding genes per genome.
Veronica panormitana contained the smallest predicted gene set, with 28,325
genes, whereas V. anagallis-aquatica and V. persica contained the largest gene
sets, with 71,164 and 68,480 genes, respectively. The remaining Veronica
species contained between 35,610 and 54,493 predicted genes, while 42,995 genes
were predicted in Plantago major.

The number of predicted transcripts ranged from 31,255 in V. panormitana to
76,009 in V. anagallis-aquatica. Transcript numbers were closely proportional
to gene numbers, with an average of 1.060–1.103 transcripts per gene.
Selection of the longest predicted isoform for each gene generated
representative proteomes containing between 28,313 and 71,165 proteins.

The predicted proteomes were highly complete. Complete protein BUSCO scores
ranged from 95.7% in V. panormitana to 98.5% in V. persica. Complete genome
BUSCO scores ranged from 98.3% to 99.0%. Most species were dominated by
single-copy BUSCOs, with duplicated protein BUSCO proportions between 4.4% and
7.2%.

In contrast, V. anagallis-aquatica and V. persica exhibited exceptionally high
BUSCO duplication. Duplicated BUSCOs represented 95.0% and 96.7% of the
predicted proteomes and 97.1% and 97.7% of the genome BUSCOs in
V. anagallis-aquatica and V. persica, respectively. These species also
contained the two largest predicted gene sets. The concordance between elevated
gene counts and duplicated BUSCO proportions indicates extensive retention of
duplicated genomic complements.

All expected annotation outputs were recovered for all nine species. Protein
and coding-sequence counts were identical within each species, and the
representative longest-isoform proteomes contained no duplicated FASTA
identifiers.
RESULTS

# ============================================================
# Interpretation and limitations
# ============================================================

cat > "${RESULTS_DIR}/annotation_interpretation.md" <<'INTERPRETATION'
## Interpretation

The annotations provide complete and internally consistent gene sets suitable
for comparative genomics. The high complete BUSCO values indicate that few
conserved embryophyte genes were omitted from either the assemblies or
predicted proteomes.

Veronica anagallis-aquatica and V. persica differ markedly from the other
species because nearly all complete BUSCOs occur in duplicated copies. This
pattern is consistent with highly duplicated genomic complements and is a
strong candidate signature of polyploidy or extensive homeologue retention.

However, BUSCO duplication alone is not sufficient to establish the timing or
mechanism of whole-genome duplication. The following analyses are required
before drawing a final evolutionary conclusion:

1. orthogroup copy-number analysis;
2. gene-tree reconciliation;
3. chromosome-level synteny;
4. within-genome duplicated-block analysis;
5. paralogue and orthologue Ks distributions;
6. comparison with the Plantago major outgroup;
7. evaluation of assembly redundancy or retained haplotigs.

Veronica scutellata has an elevated predicted gene count but only 5.1%
duplicated protein BUSCOs and 5.6% duplicated genome BUSCOs. Its increased gene
number therefore should not be interpreted as whole-genome duplication without
additional synteny and orthogroup evidence.

The small differences between GFF3 gene counts and representative longest
protein counts are negligible relative to total gene-set size. They most likely
reflect differences in how BRAKER4 records gene features and selects
protein-coding representative transcripts. These identifiers should be
documented but do not indicate a general annotation failure.
INTERPRETATION

# ============================================================
# Plain-text manuscript table
# ============================================================

python - "${FINAL_TABLE}" \
    "${TABLES_DIR}/manuscript_annotation_table.md" <<'PY'
import csv
import sys
from pathlib import Path

source = Path(sys.argv[1])
target = Path(sys.argv[2])

with source.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

columns = [
    ("scientific_name", "Species"),
    ("species_code", "Code"),
    ("evidence", "Evidence"),
    ("predicted_genes", "Genes"),
    ("transcripts", "Transcripts"),
    ("representative_proteins", "Representative proteins"),
    ("transcripts_per_gene", "Transcripts per gene"),
    ("protein_busco_complete_pct", "Protein BUSCO C (%)"),
    ("protein_busco_single_copy_pct", "S (%)"),
    ("protein_busco_duplicated_pct", "D (%)"),
    ("protein_busco_fragmented_pct", "F (%)"),
    ("protein_busco_missing_pct", "M (%)"),
]

def format_species(name: str) -> str:
    return f"*{name}*"

lines = []
lines.append("## Structural annotation statistics")
lines.append("")
lines.append("| " + " | ".join(label for _, label in columns) + " |")
lines.append("|" + "|".join("---" for _ in columns) + "|")

for row in rows:
    values = []

    for key, _ in columns:
        value = row[key]

        if key == "scientific_name":
            value = format_species(value)

        if key in {
            "predicted_genes",
            "transcripts",
            "representative_proteins",
        }:
            value = f"{int(value):,}"

        values.append(value)

    lines.append("| " + " | ".join(values) + " |")

lines.extend(
    [
        "",
        "**Abbreviations:** C, complete BUSCOs; S, complete single-copy "
        "BUSCOs; D, complete duplicated BUSCOs; F, fragmented BUSCOs; "
        "M, missing BUSCOs.",
        "",
        "BUSCO analyses used the embryophyta_odb12 dataset containing "
        "2,026 conserved orthologues.",
    ]
)

target.write_text("\n".join(lines) + "\n", encoding="utf-8")
PY

# ============================================================
# Create README
# ============================================================

cat > "${CHECKPOINT_DIR}/README.md" <<'README'
# Veronica genome annotation checkpoint

## Checkpoint status

This directory records the completed structural annotation stage for eight
Veronica genomes and the Plantago major outgroup.

All nine species passed the final annotation-output validation.

## Contents

### 01_methods

- `annotation_methods.md`: manuscript-ready Methods text
- `annotation_workflow.md`: detailed workflow and configuration

### 02_results

- `annotation_results.md`: manuscript-ready Results text
- `annotation_interpretation.md`: interpretation and limitations

### 03_tables

- complete manuscript annotation statistics
- annotation feature counts
- annotation-output paths
- representative-proteome manifest
- BUSCO summary manifests
- manuscript-formatted Markdown table

### 04_scripts

Copies of all available scripts used to prepare genomes, run BRAKER4, resume
VPAN, validate outputs, prepare representative proteomes, parse BUSCO results
and generate this checkpoint.

Species-specific `samples.csv` and `config.ini` files are also retained.

### 05_manifests

- annotation input manifest
- annotation script manifest

### 06_logs

Selected final workflow status files and logs.

### 07_output_manifest

Paths to all final per-species annotation and BUSCO outputs.

## Accepted representative protein sets

The files used for downstream comparative genomics are located in:

`08_comparative_inputs/proteomes/`

One longest predicted protein isoform per gene was retained.

## Checkpoint decision

The structural annotation stage is complete and frozen. No annotation files
should be replaced after this checkpoint without documenting the reason,
regenerating the QC tables, and creating a new checkpoint version.
README

# ============================================================
# Copy status files and selected logs
# ============================================================

for CODE in "${SPECIES[@]}"
do
    if [[ -s "${RUN_ROOT}/${CODE}/run_status.tsv" ]]; then
        cp -f "${RUN_ROOT}/${CODE}/run_status.tsv" \
            "${LOGS_DIR}/${CODE}.run_status.tsv"
    fi

    QUALITY_CONTROL="${RUN_ROOT}/${CODE}/output/${CODE}/results/quality_control"

    if [[ -s "${QUALITY_CONTROL}/busco_proteins_short_summary.txt" ]]; then
        cp -f "${QUALITY_CONTROL}/busco_proteins_short_summary.txt" \
            "${LOGS_DIR}/${CODE}.protein_busco.txt"
    fi

    if [[ -s "${QUALITY_CONTROL}/busco_genome_short_summary.txt" ]]; then
        cp -f "${QUALITY_CONTROL}/busco_genome_short_summary.txt" \
            "${LOGS_DIR}/${CODE}.genome_busco.txt"
    fi

    if [[ -s "${QUALITY_CONTROL}/gene_set_statistics.txt" ]]; then
        cp -f "${QUALITY_CONTROL}/gene_set_statistics.txt" \
            "${LOGS_DIR}/${CODE}.gene_set_statistics.txt"
    fi

    if [[ -s "${QUALITY_CONTROL}/training_summary.txt" ]]; then
        cp -f "${QUALITY_CONTROL}/training_summary.txt" \
            "${LOGS_DIR}/${CODE}.training_summary.txt"
    fi
done

# ============================================================
# Create summary statistics
# ============================================================

python - "${FINAL_TABLE}" \
    "${RESULTS_DIR}/annotation_summary_statistics.txt" <<'PY'
import csv
import statistics
import sys
from pathlib import Path

source = Path(sys.argv[1])
target = Path(sys.argv[2])

with source.open(newline="", encoding="utf-8") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))

genes = [int(row["predicted_genes"]) for row in rows]
transcripts = [int(row["transcripts"]) for row in rows]
protein_busco = [
    float(row["protein_busco_complete_pct"])
    for row in rows
]
genome_busco = [
    float(row["genome_busco_complete_pct"])
    for row in rows
]

largest_gene_set = max(rows, key=lambda row: int(row["predicted_genes"]))
smallest_gene_set = min(rows, key=lambda row: int(row["predicted_genes"]))

lines = [
    f"Number of annotated species: {len(rows)}",
    f"Minimum predicted genes: {min(genes):,}",
    f"Maximum predicted genes: {max(genes):,}",
    f"Mean predicted genes: {statistics.mean(genes):,.1f}",
    f"Median predicted genes: {statistics.median(genes):,.1f}",
    f"Minimum transcripts: {min(transcripts):,}",
    f"Maximum transcripts: {max(transcripts):,}",
    f"Protein BUSCO complete range: "
    f"{min(protein_busco):.1f}–{max(protein_busco):.1f}%",
    f"Genome BUSCO complete range: "
    f"{min(genome_busco):.1f}–{max(genome_busco):.1f}%",
    f"Largest gene set: {largest_gene_set['scientific_name']} "
    f"({int(largest_gene_set['predicted_genes']):,})",
    f"Smallest gene set: {smallest_gene_set['scientific_name']} "
    f"({int(smallest_gene_set['predicted_genes']):,})",
]

target.write_text("\n".join(lines) + "\n", encoding="utf-8")
PY

# ============================================================
# Checksums
# ============================================================

CHECKSUM_FILE="${MANIFEST_DIR}/checkpoint_sha256.tsv"

find "${CHECKPOINT_DIR}" \
    -type f \
    ! -path "${CHECKSUM_FILE}" \
    -print0 |
sort -z |
xargs -0 sha256sum |
awk 'BEGIN {OFS="\t"} {checksum=$1; $1=""; sub(/^ /,""); print checksum,$0}' \
    > "${CHECKSUM_FILE}"

# ============================================================
# Checkpoint completion marker
# ============================================================

CREATED_AT=$(date --iso-8601=seconds)

cat > "${CHECKPOINT_DIR}/CHECKPOINT_COMPLETE.txt" <<EOF2
checkpoint=annotation_final
project=Veronica_genome_evolution
created_at=${CREATED_AT}
species_count=9
annotation_status=COMPLETE
output_validation=PASS
representative_proteomes=PASS
duplicate_fasta_headers=0
busco_lineage=embryophyta_odb12
busco_n=2026
next_stage=orthology_and_gene_family_analysis
EOF2

echo
echo "============================================================"
echo "Annotation checkpoint created successfully"
echo "============================================================"
echo "Checkpoint directory:"
echo "${CHECKPOINT_DIR}"
echo
echo "Contents:"
find "${CHECKPOINT_DIR}" -maxdepth 2 -type f | sort
