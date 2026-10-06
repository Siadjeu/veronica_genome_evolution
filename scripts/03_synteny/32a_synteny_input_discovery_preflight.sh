#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=synteny_preflight
#SBATCH --output=10_synteny/logs/synteny_preflight_%j.out
#SBATCH --error=10_synteny/logs/synteny_preflight_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

SYNTENY_DIR="${PROJECT_DIR}/10_synteny"
DISCOVERY_DIR="${SYNTENY_DIR}/input_discovery"
TABLE_DIR="${DISCOVERY_DIR}/tables"
MANIFEST_DIR="${SYNTENY_DIR}/manifests"
QC_DIR="${SYNTENY_DIR}/qc"
LOG_DIR="${SYNTENY_DIR}/logs"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_input_discovery"

PROTEOME_DIR="${PROJECT_DIR}/08_comparative_inputs/proteomes"

mkdir -p \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" \
    "${QC_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

cd "${PROJECT_DIR}"

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
    echo "ERROR: Conda was not found." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate genespace_env

PREVIOUS_CHECKPOINT="${PROJECT_DIR}/09_orthology/checkpoint_gene_tree_topology_recovered/GENE_TREE_RECOVERY_COMPLETE.txt"

if [[ ! -s "${PREVIOUS_CHECKPOINT}" ]] ||
   ! grep -q '^status=PASS$' "${PREVIOUS_CHECKPOINT}"; then
    echo "ERROR: Step 31B checkpoint is missing or invalid." >&2
    exit 1
fi

rm -f "${TABLE_DIR}"/*.tsv
rm -f "${MANIFEST_DIR}"/synteny_input_manifest.proposed.tsv
rm -f "${QC_DIR}"/synteny_input_discovery_summary.tsv
rm -f "${CHECKPOINT_DIR}"/*

python - \
    "${PROJECT_DIR}" \
    "${PROTEOME_DIR}" \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" \
    "${QC_DIR}" <<'PY'
from __future__ import annotations

import csv
import gzip
import os
import re
import sys
from collections import Counter
from pathlib import Path

(
    project_dir_name,
    proteome_dir_name,
    table_dir_name,
    manifest_dir_name,
    qc_dir_name,
) = sys.argv[1:]

project_dir = Path(project_dir_name).resolve()
proteome_dir = Path(proteome_dir_name).resolve()
table_dir = Path(table_dir_name)
manifest_dir = Path(manifest_dir_name)
qc_dir = Path(qc_dir_name)

species = {
    "VPAN": ("Veronica panormitana", "2x", "18", "diploid_Veronica"),
    "VSCU": ("Veronica scutellata", "2x", "18", "diploid_Veronica"),
    "VANA": ("Veronica anagallis-aquatica", "4x", "36", "tetraploid_Veronica"),
    "VARV": ("Veronica arvensis", "2x", "16", "diploid_Veronica"),
    "VPER": ("Veronica persica", "4x", "28", "tetraploid_Veronica"),
    "VSER": ("Veronica serpyllifolia", "2x", "14", "diploid_Veronica"),
    "VTRI": ("Veronica triloba", "2x", "18", "diploid_Veronica"),
    "VVER": ("Veronica verna", "2x", "16", "diploid_Veronica"),
    "PMAJ": ("Plantago major", "2x", "12", "diploid_outgroup"),
}

excluded_directory_names = {
    ".git",
    "09_orthology",
    "10_synteny",
    "08_comparative_inputs",
    "proteomes",
    "proteins",
    "orthofinder",
    "busco_downloads",
}

annotation_suffixes = (
    ".gff",
    ".gff3",
    ".gtf",
    ".gff.gz",
    ".gff3.gz",
    ".gtf.gz",
)

fasta_suffixes = (
    ".fa",
    ".fasta",
    ".fna",
    ".fas",
    ".fa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fas.gz",
)

protein_name_tokens = {
    "protein",
    "proteins",
    "proteome",
    "pep",
    "peptide",
    "amino",
    "longest.aa",
    ".aa.",
    "cds",
    "mrna",
    "transcript",
    "rna",
}

genome_positive_tokens = {
    "genomic",
    "genome",
    "assembly",
    "chromosome",
    "chromosomes",
    "primary",
    "final",
    "softmasked",
    "masked",
}

def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace",
        )

    return path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    )

def contains_species_code(path: Path, code: str) -> bool:
    text = str(path).upper()

    patterns = [
        f"/{code}/",
        f"/{code}_",
        f"_{code}_",
        f"/{code}.",
        f"{code}.",
        f"{code}_",
    ]

    return any(pattern in text for pattern in patterns)

def fasta_sequence_type(path: Path, maximum_letters: int = 100000) -> tuple[str, float]:
    nucleotide = set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv")
    amino_acid_only = set("EFILPQZefilpqz")

    nucleotide_count = 0
    amino_indicator_count = 0
    total = 0

    try:
        with open_text(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    continue

                sequence = re.sub(r"[^A-Za-z]", "", line)

                for character in sequence:
                    total += 1

                    if character in nucleotide:
                        nucleotide_count += 1

                    if character in amino_acid_only:
                        amino_indicator_count += 1

                    if total >= maximum_letters:
                        break

                if total >= maximum_letters:
                    break
    except Exception:
        return "UNREADABLE", 0.0

    if total == 0:
        return "EMPTY", 0.0

    nucleotide_fraction = nucleotide_count / total
    amino_indicator_fraction = amino_indicator_count / total

    if nucleotide_fraction >= 0.95 and amino_indicator_fraction <= 0.01:
        return "NUCLEOTIDE", nucleotide_fraction

    return "PROTEIN_OR_NONNUCLEOTIDE", nucleotide_fraction

def fasta_stats(path: Path) -> dict:
    count = 0
    lengths = []
    current_length = 0
    first_headers = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if count:
                    lengths.append(current_length)

                count += 1
                current_length = 0

                if len(first_headers) < 10:
                    first_headers.append(
                        line[1:].strip().split()[0]
                    )
            else:
                current_length += len(
                    re.sub(r"\s+", "", line)
                )

    if count:
        lengths.append(current_length)

    total = sum(lengths)
    ordered = sorted(lengths, reverse=True)

    cumulative = 0
    n50 = 0

    for length in ordered:
        cumulative += length

        if total and cumulative >= total / 2:
            n50 = length
            break

    return {
        "sequence_count": count,
        "total_bp": total,
        "n50_bp": n50,
        "max_bp": max(lengths, default=0),
        "first_headers": ",".join(first_headers),
    }

def annotation_stats(path: Path) -> dict:
    feature_counts = Counter()
    seqids = set()
    malformed = 0

    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                malformed += 1
                continue

            seqids.add(fields[0])
            feature_counts[fields[2].lower()] += 1

    return {
        "seqid_count": len(seqids),
        "gene_count": feature_counts["gene"],
        "mrna_count": (
            feature_counts["mrna"]
            + feature_counts["transcript"]
        ),
        "cds_count": feature_counts["cds"],
        "malformed_lines": malformed,
        "first_seqids": ",".join(sorted(seqids)[:10]),
    }

def genome_score(path: Path, code: str, stats: dict) -> int:
    lower = str(path).lower()
    name = path.name.lower()

    score = 0

    if contains_species_code(path, code):
        score += 50

    if code.lower() in name:
        score += 20

    if any(token in lower for token in genome_positive_tokens):
        score += 15

    if "genomic.fna" in name:
        score += 30

    if "chromosome" in lower:
        score += 20

    if "braker" in lower:
        score -= 20

    if stats["total_bp"] >= 100_000_000:
        score += 30
    elif stats["total_bp"] >= 10_000_000:
        score += 10
    else:
        score -= 30

    if stats["max_bp"] >= 1_000_000:
        score += 30

    if stats["n50_bp"] >= 1_000_000:
        score += 20

    return score

def annotation_score(path: Path, code: str) -> int:
    text = str(path).lower()
    name = path.name.lower()

    score = 0

    if contains_species_code(path, code):
        score += 50

    if code.lower() in name:
        score += 20

    if "braker.longest.gtf" in name:
        score += 40
    elif "braker.gtf" in name:
        score += 30
    elif "braker.gff3" in name:
        score += 25

    if "old" in text or "backup" in text:
        score -= 20

    return score

# ============================================================
# Discover nucleotide genome FASTAs
# ============================================================

genome_candidates = []
annotation_candidates = []

for root, directory_names, file_names in os.walk(project_dir):
    directory_names[:] = [
        directory
        for directory in directory_names
        if directory.lower() not in excluded_directory_names
        and not directory.startswith(".")
    ]

    root_path = Path(root)

    for file_name in file_names:
        path = (root_path / file_name).resolve()
        lower_name = path.name.lower()

        if lower_name.endswith(annotation_suffixes):
            for code in species:
                if contains_species_code(path, code):
                    annotation_candidates.append(
                        {
                            "species_code": code,
                            "score": annotation_score(path, code),
                            "size_bytes": path.stat().st_size,
                            "path": str(path),
                        }
                    )
            continue

        if not lower_name.endswith(fasta_suffixes):
            continue

        if proteome_dir in path.parents:
            continue

        if any(token in lower_name for token in protein_name_tokens):
            continue

        sequence_type, nucleotide_fraction = fasta_sequence_type(path)

        if sequence_type != "NUCLEOTIDE":
            continue

        try:
            stats = fasta_stats(path)
        except Exception:
            continue

        # Genome assemblies should normally be substantially larger than
        # transcript/CDS datasets and include at least one long sequence.
        if stats["total_bp"] < 50_000_000:
            continue

        if stats["max_bp"] < 100_000:
            continue

        for code in species:
            if not contains_species_code(path, code):
                continue

            genome_candidates.append(
                {
                    "species_code": code,
                    "score": genome_score(path, code, stats),
                    "size_bytes": path.stat().st_size,
                    "nucleotide_fraction": f"{nucleotide_fraction:.5f}",
                    "sequence_count": stats["sequence_count"],
                    "total_bp": stats["total_bp"],
                    "n50_bp": stats["n50_bp"],
                    "max_bp": stats["max_bp"],
                    "path": str(path),
                }
            )

genome_candidates.sort(
    key=lambda row: (
        row["species_code"],
        -row["score"],
        -row["n50_bp"],
        -row["total_bp"],
        row["path"],
    )
)

annotation_candidates.sort(
    key=lambda row: (
        row["species_code"],
        -row["score"],
        -row["size_bytes"],
        row["path"],
    )
)

with (
    table_dir / "discovered_genome_candidates.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    fields = [
        "species_code",
        "score",
        "size_bytes",
        "nucleotide_fraction",
        "sequence_count",
        "total_bp",
        "n50_bp",
        "max_bp",
        "path",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(genome_candidates)

with (
    table_dir / "discovered_annotation_candidates.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    fields = [
        "species_code",
        "score",
        "size_bytes",
        "path",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(annotation_candidates)

# ============================================================
# Select highest-ranked inputs
# ============================================================

selected_genomes = {}
selected_annotations = {}

for code in species:
    genome_rows = [
        row
        for row in genome_candidates
        if row["species_code"] == code
    ]

    annotation_rows = [
        row
        for row in annotation_candidates
        if row["species_code"] == code
    ]

    selected_genomes[code] = (
        Path(genome_rows[0]["path"])
        if genome_rows
        else None
    )

    selected_annotations[code] = (
        Path(annotation_rows[0]["path"])
        if annotation_rows
        else None
    )

preflight_rows = []

for code, metadata in species.items():
    scientific_name, ploidy, expected_2n, group = metadata

    genome = selected_genomes[code]
    annotation = selected_annotations[code]
    protein = proteome_dir / f"{code}.longest.aa.fa"

    warnings = []

    genome_values = {
        "sequence_count": 0,
        "total_bp": 0,
        "n50_bp": 0,
        "max_bp": 0,
        "first_headers": "",
    }

    annotation_values = {
        "seqid_count": 0,
        "gene_count": 0,
        "mrna_count": 0,
        "cds_count": 0,
        "malformed_lines": 0,
        "first_seqids": "",
    }

    if genome and genome.is_file():
        genome_values = fasta_stats(genome)
    else:
        warnings.append("MISSING_GENOME")

    if annotation and annotation.is_file():
        annotation_values = annotation_stats(annotation)
    else:
        warnings.append("MISSING_ANNOTATION")

    protein_count = 0

    if protein.is_file():
        protein_count = fasta_stats(protein)["sequence_count"]
    else:
        warnings.append("MISSING_PROTEIN")

    if genome:
        sequence_type, _ = fasta_sequence_type(genome)

        if sequence_type != "NUCLEOTIDE":
            warnings.append("GENOME_NOT_NUCLEOTIDE")

        if genome_values["total_bp"] < 50_000_000:
            warnings.append("GENOME_TOO_SMALL")

        if genome_values["max_bp"] < 100_000:
            warnings.append("NO_LONG_GENOMIC_SEQUENCE")

        if proteome_dir in genome.resolve().parents:
            warnings.append("GENOME_IS_PROTEOME_FILE")

    if annotation and annotation_values["cds_count"] == 0:
        warnings.append("NO_CDS_FEATURES")

    if annotation_values["malformed_lines"] > 0:
        warnings.append("MALFORMED_ANNOTATION_LINES")

    status = "PASS" if not warnings else "REVIEW"

    preflight_rows.append(
        {
            "species_code": code,
            "scientific_name": scientific_name,
            "ploidy": ploidy,
            "expected_2n": expected_2n,
            "analysis_group": group,
            "genome_path": str(genome) if genome else "",
            "annotation_path": str(annotation) if annotation else "",
            "protein_path": str(protein.resolve()) if protein.exists() else "",
            "genome_sequence_count": genome_values["sequence_count"],
            "genome_size_bp": genome_values["total_bp"],
            "genome_N50_bp": genome_values["n50_bp"],
            "genome_max_sequence_bp": genome_values["max_bp"],
            "annotation_seqid_count": annotation_values["seqid_count"],
            "gene_feature_count": annotation_values["gene_count"],
            "mrna_feature_count": annotation_values["mrna_count"],
            "cds_feature_count": annotation_values["cds_count"],
            "protein_sequence_count": protein_count,
            "status": status,
            "warnings": ";".join(warnings),
        }
    )

fields = [
    "species_code",
    "scientific_name",
    "ploidy",
    "expected_2n",
    "analysis_group",
    "genome_path",
    "annotation_path",
    "protein_path",
    "genome_sequence_count",
    "genome_size_bp",
    "genome_N50_bp",
    "genome_max_sequence_bp",
    "annotation_seqid_count",
    "gene_feature_count",
    "mrna_feature_count",
    "cds_feature_count",
    "protein_sequence_count",
    "status",
    "warnings",
]

with (
    table_dir / "synteny_input_preflight.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(preflight_rows)

with (
    manifest_dir / "synteny_input_manifest.proposed.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow(
        [
            "species_code",
            "scientific_name",
            "ploidy",
            "analysis_group",
            "genome_fasta",
            "annotation_gff",
            "protein_fasta",
            "approved",
            "notes",
        ]
    )

    for row in preflight_rows:
        writer.writerow(
            [
                row["species_code"],
                row["scientific_name"],
                row["ploidy"],
                row["analysis_group"],
                row["genome_path"],
                row["annotation_path"],
                row["protein_path"],
                "NO",
                row["warnings"],
            ]
        )

missing_genomes = sum(
    not row["genome_path"]
    for row in preflight_rows
)

missing_annotations = sum(
    not row["annotation_path"]
    for row in preflight_rows
)

missing_proteins = sum(
    not row["protein_path"]
    for row in preflight_rows
)

review_species = sum(
    row["status"] != "PASS"
    for row in preflight_rows
)

with (
    qc_dir / "synteny_input_discovery_summary.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow(["metric", "value"])
    writer.writerow(["species_expected", len(species)])
    writer.writerow(["genome_candidates_found", len(genome_candidates)])
    writer.writerow(["annotation_candidates_found", len(annotation_candidates)])
    writer.writerow(["missing_genomes", missing_genomes])
    writer.writerow(["missing_annotations", missing_annotations])
    writer.writerow(["missing_proteins", missing_proteins])
    writer.writerow(["species_requiring_review", review_species])

print("Corrected synteny input discovery completed.")
print("Genome candidates:", len(genome_candidates))
print("Annotation candidates:", len(annotation_candidates))
print("Missing genomes:", missing_genomes)
print("Missing annotations:", missing_annotations)
print("Species requiring review:", review_species)
PY

echo
echo "============================================================"
echo "Corrected discovery summary"
echo "============================================================"

column -t -s $'\t' \
    "${QC_DIR}/synteny_input_discovery_summary.tsv"

echo
echo "============================================================"
echo "Corrected proposed inputs"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/synteny_input_preflight.tsv"

cat > "${CHECKPOINT_DIR}/SYNTENY_INPUT_DISCOVERY_COMPLETE.txt" <<EOF2
checkpoint=synteny_input_discovery
date=$(date --iso-8601=seconds)
status=PASS
next_step=review_corrected_genome_and_annotation_paths
EOF2

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "Corrected Step 32A completed."
