#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=synteny_inputs
#SBATCH --output=10_synteny/logs/synteny_inputs_%j.out
#SBATCH --error=10_synteny/logs/synteny_inputs_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

PREFLIGHT_TABLE="${SYNTENY_DIR}/input_discovery/tables/synteny_input_preflight.tsv"

OUTPUT_DIR="${SYNTENY_DIR}/prepared_inputs"
GENOME_DIR="${OUTPUT_DIR}/genomes"
GTF_DIR="${OUTPUT_DIR}/annotations"
PROTEIN_DIR="${OUTPUT_DIR}/proteomes"
TABLE_DIR="${OUTPUT_DIR}/tables"
LOG_DIR="${SYNTENY_DIR}/logs"

MANIFEST_DIR="${SYNTENY_DIR}/manifests"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_prepared_inputs"

mkdir -p \
    "${GENOME_DIR}" \
    "${GTF_DIR}" \
    "${PROTEIN_DIR}" \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

cd "${PROJECT_DIR}"

if [[ ! -s "${PREFLIGHT_TABLE}" ]]; then
    echo "ERROR: Step 32A preflight table is missing:" >&2
    echo "${PREFLIGHT_TABLE}" >&2
    exit 1
fi

STEP32A_CHECKPOINT="${SYNTENY_DIR}/checkpoint_input_discovery/SYNTENY_INPUT_DISCOVERY_COMPLETE.txt"

if [[ ! -s "${STEP32A_CHECKPOINT}" ]] ||
   ! grep -q '^status=PASS$' "${STEP32A_CHECKPOINT}"; then
    echo "ERROR: Step 32A checkpoint is missing or invalid." >&2
    exit 1
fi

rm -f "${GENOME_DIR}"/*
rm -f "${GTF_DIR}"/*
rm -f "${PROTEIN_DIR}"/*
rm -f "${TABLE_DIR}"/*
rm -f "${CHECKPOINT_DIR}"/*

python - \
    "${PREFLIGHT_TABLE}" \
    "${GENOME_DIR}" \
    "${GTF_DIR}" \
    "${PROTEIN_DIR}" \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" <<'PY'
from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

(
    preflight_name,
    genome_dir_name,
    gtf_dir_name,
    protein_dir_name,
    table_dir_name,
    manifest_dir_name,
) = sys.argv[1:]

preflight_file = Path(preflight_name)
genome_dir = Path(genome_dir_name)
gtf_dir = Path(gtf_dir_name)
protein_dir = Path(protein_dir_name)
table_dir = Path(table_dir_name)
manifest_dir = Path(manifest_dir_name)

for directory in [
    genome_dir,
    gtf_dir,
    protein_dir,
    table_dir,
    manifest_dir,
]:
    directory.mkdir(parents=True, exist_ok=True)

expected_haploid_chromosomes = {
    "VPAN": 9,
    "VSCU": 9,
    "VANA": 18,
    "VARV": 8,
    "VPER": 14,
    "VSER": 7,
    "VTRI": 9,
    "VVER": 8,
    "PMAJ": 6,
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

def write_text(path: Path, compress: bool = False):
    if compress:
        return gzip.open(
            path,
            "wt",
            encoding="utf-8",
        )

    return path.open(
        "w",
        encoding="utf-8",
    )

def fasta_headers(path: Path) -> list[str]:
    headers = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                headers.append(
                    line[1:].strip().split()[0]
                )

    return headers

def parse_gtf_attributes(text: str) -> dict[str, str]:
    attributes = {}

    for match in re.finditer(
        r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"',
        text,
    ):
        attributes[match.group(1)] = match.group(2)

    for part in text.strip().strip(";").split(";"):
        part = part.strip()

        if not part or " " in part:
            continue

        if "=" in part:
            key, value = part.split("=", 1)
            attributes.setdefault(key, value)

    return attributes

def normalize_identifier(identifier: str) -> str:
    value = identifier.strip()

    value = re.sub(
        r"^(transcript:|protein:|gene:)",
        "",
        value,
    )

    return value

def identifier_variants(identifier: str) -> set[str]:
    value = normalize_identifier(identifier)

    variants = {
        value,
        value.split()[0],
    }

    if "|" in value:
        parts = value.split("|")
        variants.update(parts)
        variants.add(parts[-1])

    variants.add(
        re.sub(r"\.p\d+$", "", value)
    )

    return {
        item
        for item in variants
        if item
    }

def read_proteins(path: Path):
    records = []
    current_header = None
    current_sequence = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_header is not None:
                    records.append(
                        (
                            current_header,
                            "".join(current_sequence),
                        )
                    )

                current_header = line[1:].rstrip("\n")
                current_sequence = []
            else:
                current_sequence.append(line.strip())

    if current_header is not None:
        records.append(
            (
                current_header,
                "".join(current_sequence),
            )
        )

    return records

with preflight_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    preflight_rows = list(
        csv.DictReader(handle, delimiter="\t")
    )

summary_rows = []
seqid_rows = []
manifest_rows = []

for row in preflight_rows:
    code = row["species_code"]

    genome_source = Path(row["genome_path"])
    annotation_source = Path(row["annotation_path"])
    protein_source = Path(row["protein_path"])

    if not genome_source.is_file():
        raise SystemExit(
            f"ERROR: Missing genome for {code}: "
            f"{genome_source}"
        )

    if not annotation_source.is_file():
        raise SystemExit(
            f"ERROR: Missing annotation for {code}: "
            f"{annotation_source}"
        )

    if not protein_source.is_file():
        raise SystemExit(
            f"ERROR: Missing proteome for {code}: "
            f"{protein_source}"
        )

    chromosome_headers = fasta_headers(
        genome_source
    )
    chromosome_set = set(chromosome_headers)

    genome_output = (
        genome_dir / f"{code}.chromosomes.fa.gz"
    )

    if genome_source.resolve() != genome_output.resolve():
        with open_text(genome_source) as source, gzip.open(
            genome_output,
            "wt",
            encoding="utf-8",
        ) as destination:
            shutil.copyfileobj(source, destination)

    all_annotation_seqids = set()
    retained_annotation_seqids = set()

    all_gene_ids = set()
    retained_gene_ids = set()

    all_transcript_ids = set()
    retained_transcript_ids = set()

    total_feature_count = 0
    retained_feature_count = 0
    total_cds_count = 0
    retained_cds_count = 0

    annotation_output = (
        gtf_dir / f"{code}.chromosomes.longest.gtf.gz"
    )

    with open_text(annotation_source) as source, gzip.open(
        annotation_output,
        "wt",
        encoding="utf-8",
    ) as destination:
        for line in source:
            if line.startswith("#") or not line.strip():
                destination.write(line)
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            seqid = fields[0]
            feature = fields[2].lower()
            attributes = parse_gtf_attributes(
                fields[8]
            )

            all_annotation_seqids.add(seqid)
            total_feature_count += 1

            gene_id = (
                attributes.get("gene_id")
                or attributes.get("gene")
                or attributes.get("Parent")
                or ""
            )

            transcript_id = (
                attributes.get("transcript_id")
                or attributes.get("transcript")
                or attributes.get("ID")
                or attributes.get("Parent")
                or ""
            )

            if gene_id:
                all_gene_ids.add(
                    normalize_identifier(gene_id)
                )

            if feature in {
                "transcript",
                "mrna",
                "cds",
                "exon",
            } and transcript_id:
                retained_candidate = normalize_identifier(
                    transcript_id.split(",")[0]
                )
                all_transcript_ids.add(
                    retained_candidate
                )

            if feature == "cds":
                total_cds_count += 1

            if seqid not in chromosome_set:
                continue

            destination.write(line)
            retained_annotation_seqids.add(seqid)
            retained_feature_count += 1

            if gene_id:
                retained_gene_ids.add(
                    normalize_identifier(gene_id)
                )

            if feature in {
                "transcript",
                "mrna",
                "cds",
                "exon",
            } and transcript_id:
                retained_transcript_ids.add(
                    normalize_identifier(
                        transcript_id.split(",")[0]
                    )
                )

            if feature == "cds":
                retained_cds_count += 1

    protein_records = read_proteins(
        protein_source
    )

    protein_lookup = {}

    for header, sequence in protein_records:
        primary_id = header.split()[0]

        for variant in identifier_variants(
            primary_id
        ):
            protein_lookup.setdefault(
                variant,
                (header, sequence),
            )

    selected_proteins = {}
    unmatched_transcripts = []

    for transcript_id in sorted(
        retained_transcript_ids
    ):
        match = None

        for variant in identifier_variants(
            transcript_id
        ):
            if variant in protein_lookup:
                match = protein_lookup[variant]
                break

        if match is None:
            unmatched_transcripts.append(
                transcript_id
            )
            continue

        header, sequence = match
        selected_proteins[header.split()[0]] = (
            header,
            sequence,
        )

    protein_output = (
        protein_dir / f"{code}.chromosomes.longest.aa.fa"
    )

    with protein_output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for primary_id in sorted(
            selected_proteins
        ):
            header, sequence = selected_proteins[
                primary_id
            ]

            handle.write(f">{header}\n")

            for start in range(
                0,
                len(sequence),
                60,
            ):
                handle.write(
                    sequence[start:start + 60]
                    + "\n"
                )

    missing_chromosome_annotations = (
        chromosome_set
        - retained_annotation_seqids
    )

    annotation_only_seqids = (
        all_annotation_seqids
        - chromosome_set
    )

    transcript_match_pct = (
        100.0
        * len(selected_proteins)
        / len(retained_transcript_ids)
        if retained_transcript_ids
        else 0.0
    )

    feature_retention_pct = (
        100.0
        * retained_feature_count
        / total_feature_count
        if total_feature_count
        else 0.0
    )

    warnings = []

    expected_chromosomes = (
        expected_haploid_chromosomes[code]
    )

    if len(chromosome_headers) != expected_chromosomes:
        warnings.append(
            "CHROMOSOME_COUNT_DIFFERS_FROM_EXPECTED"
        )

    if missing_chromosome_annotations:
        warnings.append(
            "CHROMOSOME_WITHOUT_ANNOTATION"
        )

    if transcript_match_pct < 95.0:
        warnings.append(
            "LOW_TRANSCRIPT_PROTEIN_ID_MATCH"
        )

    if len(selected_proteins) == 0:
        warnings.append(
            "NO_FILTERED_PROTEINS"
        )

    status = (
        "PASS"
        if not warnings
        else "REVIEW"
    )

    for seqid in chromosome_headers:
        seqid_rows.append(
            {
                "species_code": code,
                "genome_seqid": seqid,
                "present_in_annotation": (
                    seqid in retained_annotation_seqids
                ),
            }
        )

    summary_rows.append(
        {
            "species_code": code,
            "expected_haploid_chromosomes": expected_chromosomes,
            "observed_genome_sequences": len(chromosome_headers),
            "annotation_seqids_total": len(all_annotation_seqids),
            "annotation_seqids_on_chromosomes": len(retained_annotation_seqids),
            "annotation_only_seqids_removed": len(annotation_only_seqids),
            "genome_seqids_without_annotation": len(missing_chromosome_annotations),
            "features_total": total_feature_count,
            "features_retained": retained_feature_count,
            "feature_retention_pct": f"{feature_retention_pct:.2f}",
            "cds_total": total_cds_count,
            "cds_retained": retained_cds_count,
            "transcripts_retained": len(retained_transcript_ids),
            "proteins_retained": len(selected_proteins),
            "unmatched_retained_transcripts": len(unmatched_transcripts),
            "transcript_protein_match_pct": f"{transcript_match_pct:.2f}",
            "status": status,
            "warnings": ";".join(warnings),
        }
    )

    unmatched_output = (
        table_dir
        / f"{code}.unmatched_transcript_ids.txt"
    )

    unmatched_output.write_text(
        "\n".join(unmatched_transcripts)
        + ("\n" if unmatched_transcripts else ""),
        encoding="utf-8",
    )

    manifest_rows.append(
        {
            "species_code": code,
            "scientific_name": row["scientific_name"],
            "ploidy": row["ploidy"],
            "analysis_group": row["analysis_group"],
            "genome_fasta": str(genome_output),
            "annotation_gtf": str(annotation_output),
            "protein_fasta": str(protein_output),
            "chromosome_count": len(chromosome_headers),
            "gene_or_transcript_count": len(selected_proteins),
            "status": status,
            "approved": "NO",
            "notes": ";".join(warnings),
        }
    )

summary_fields = [
    "species_code",
    "expected_haploid_chromosomes",
    "observed_genome_sequences",
    "annotation_seqids_total",
    "annotation_seqids_on_chromosomes",
    "annotation_only_seqids_removed",
    "genome_seqids_without_annotation",
    "features_total",
    "features_retained",
    "feature_retention_pct",
    "cds_total",
    "cds_retained",
    "transcripts_retained",
    "proteins_retained",
    "unmatched_retained_transcripts",
    "transcript_protein_match_pct",
    "status",
    "warnings",
]

with (
    table_dir / "chromosome_input_compatibility.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(summary_rows)

with (
    table_dir / "chromosome_seqid_overlap.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "genome_seqid",
            "present_in_annotation",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(seqid_rows)

manifest_fields = [
    "species_code",
    "scientific_name",
    "ploidy",
    "analysis_group",
    "genome_fasta",
    "annotation_gtf",
    "protein_fasta",
    "chromosome_count",
    "gene_or_transcript_count",
    "status",
    "approved",
    "notes",
]

with (
    manifest_dir / "synteny_input_manifest.chromosome_filtered.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=manifest_fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)

pass_count = sum(
    row["status"] == "PASS"
    for row in summary_rows
)

review_count = len(summary_rows) - pass_count

with (
    table_dir / "chromosome_input_summary.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["metric", "value"])
    writer.writerow(["species_total", len(summary_rows)])
    writer.writerow(["species_pass", pass_count])
    writer.writerow(["species_review", review_count])

print("Chromosome-compatible inputs prepared.")
print("Species total:", len(summary_rows))
print("PASS:", pass_count)
print("REVIEW:", review_count)
PY

echo
echo "============================================================"
echo "Chromosome-input summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/chromosome_input_summary.tsv"

echo
echo "============================================================"
echo "Coordinate and protein compatibility"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/chromosome_input_compatibility.tsv"

echo
echo "============================================================"
echo "Prepared manifest"
echo "============================================================"

column -t -s $'\t' \
    "${MANIFEST_DIR}/synteny_input_manifest.chromosome_filtered.tsv"

cat > "${CHECKPOINT_DIR}/CHROMOSOME_SYNTENY_INPUTS_PREPARED.txt" <<EOF2
checkpoint=chromosome_synteny_inputs_prepared
date=$(date --iso-8601=seconds)
status=PASS
next_step=review_coordinate_and_identifier_compatibility
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
echo "Step 32B completed successfully."
