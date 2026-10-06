#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=rebuild_synteny
#SBATCH --output=10_synteny/logs/rebuild_synteny_%j.out
#SBATCH --error=10_synteny/logs/rebuild_synteny_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

PREFLIGHT="${SYNTENY_DIR}/input_discovery/tables/synteny_input_preflight.tsv"
MAPPING_TABLE="${SYNTENY_DIR}/seqid_recovery/tables/exact_chromosome_seqid_mapping.tsv"

OUTPUT_DIR="${SYNTENY_DIR}/verified_inputs"
GENOME_DIR="${OUTPUT_DIR}/genomes"
ANNOTATION_DIR="${OUTPUT_DIR}/annotations"
PROTEIN_DIR="${OUTPUT_DIR}/proteomes"
TABLE_DIR="${OUTPUT_DIR}/tables"

MANIFEST_DIR="${SYNTENY_DIR}/manifests"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_verified_inputs"
LOG_DIR="${SYNTENY_DIR}/logs"

mkdir -p \
    "${GENOME_DIR}" \
    "${ANNOTATION_DIR}" \
    "${PROTEIN_DIR}" \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

cd "${PROJECT_DIR}"

for FILE in \
    "${PREFLIGHT}" \
    "${MAPPING_TABLE}" \
    "${SYNTENY_DIR}/checkpoint_seqid_recovery/EXACT_SEQID_MAPPING_COMPLETE.txt"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required input is missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

rm -f "${GENOME_DIR}"/*
rm -f "${ANNOTATION_DIR}"/*
rm -f "${PROTEIN_DIR}"/*
rm -f "${TABLE_DIR}"/*
rm -f "${CHECKPOINT_DIR}"/*

python - \
    "${PREFLIGHT}" \
    "${MAPPING_TABLE}" \
    "${GENOME_DIR}" \
    "${ANNOTATION_DIR}" \
    "${PROTEIN_DIR}" \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" <<'PY'
from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

(
    preflight_name,
    mapping_name,
    genome_dir_name,
    annotation_dir_name,
    protein_dir_name,
    table_dir_name,
    manifest_dir_name,
) = sys.argv[1:]

preflight_file = Path(preflight_name)
mapping_file = Path(mapping_name)

genome_dir = Path(genome_dir_name)
annotation_dir = Path(annotation_dir_name)
protein_dir = Path(protein_dir_name)
table_dir = Path(table_dir_name)
manifest_dir = Path(manifest_dir_name)

for directory in [
    genome_dir,
    annotation_dir,
    protein_dir,
    table_dir,
    manifest_dir,
]:
    directory.mkdir(parents=True, exist_ok=True)

expected_chromosome_counts = {
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


def fasta_records(path: Path):
    header = None
    sequence_parts = []

    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

    if header is not None:
        yield header, "".join(sequence_parts)


def primary_id(header: str) -> str:
    return header.split()[0]


def normalize_identifier(value: str) -> str:
    value = value.strip()

    value = re.sub(
        r"^(gene:|transcript:|protein:)",
        "",
        value,
    )

    return value


def parse_attributes(text: str) -> dict[str, str]:
    attributes = {}

    for match in re.finditer(
        r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"',
        text,
    ):
        attributes[match.group(1)] = match.group(2)

    for part in text.strip().strip(";").split(";"):
        part = part.strip()

        if "=" not in part:
            continue

        key, value = part.split("=", 1)
        attributes.setdefault(
            key.strip(),
            value.strip(),
        )

    return attributes


def transcript_from_attributes(
    feature: str,
    attributes: dict[str, str],
) -> str:
    candidates = []

    if feature in {"mrna", "transcript"}:
        candidates.extend(
            [
                attributes.get("transcript_id"),
                attributes.get("ID"),
            ]
        )
    else:
        candidates.extend(
            [
                attributes.get("transcript_id"),
                attributes.get("Parent"),
                attributes.get("ID"),
            ]
        )

    for value in candidates:
        if not value:
            continue

        value = value.split(",")[0]
        return normalize_identifier(value)

    return ""


def gene_from_attributes(
    attributes: dict[str, str],
) -> str:
    for key in [
        "gene_id",
        "gene",
        "gene_name",
    ]:
        value = attributes.get(key)

        if value:
            return normalize_identifier(
                value.split(",")[0]
            )

    return ""


def identifier_variants(value: str) -> set[str]:
    value = normalize_identifier(value)
    variants = {value}

    if "|" in value:
        parts = value.split("|")
        variants.update(parts)
        variants.add(parts[-1])

    variants.add(
        re.sub(r"\.p\d+$", "", value)
    )

    variants.add(
        re.sub(r"\.protein\d+$", "", value)
    )

    return {
        variant
        for variant in variants
        if variant
    }


with preflight_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    preflight_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

with mapping_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    mapping_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

mapping_by_species = defaultdict(dict)

for row in mapping_rows:
    if row["status"] != "MAPPED":
        continue

    mapping_by_species[
        row["species_code"]
    ][row["braker_seqid"]] = row[
        "chromosome_fasta_id"
    ]

summary_rows = []
manifest_rows = []
seqid_validation_rows = []

for row in preflight_rows:
    code = row["species_code"]

    genome_source = Path(row["genome_path"])
    annotation_source = Path(
        row["annotation_path"]
    )
    protein_source = Path(row["protein_path"])

    seqid_map = mapping_by_species.get(
        code,
        {},
    )

    if not seqid_map:
        raise SystemExit(
            f"ERROR: No verified mapping for {code}"
        )

    genome_output = (
        genome_dir
        / f"{code}.chromosomes.fa.gz"
    )

    with open_text(genome_source) as source, gzip.open(
        genome_output,
        "wt",
        encoding="utf-8",
    ) as destination:
        shutil.copyfileobj(
            source,
            destination,
        )

    genome_headers = [
        primary_id(header)
        for header, _ in fasta_records(
            genome_source
        )
    ]

    genome_header_set = set(
        genome_headers
    )

    annotation_output = (
        annotation_dir
        / f"{code}.chromosomes.verified.gtf.gz"
    )

    retained_transcripts = set()
    retained_genes = set()
    retained_seqids = set()

    total_features = 0
    retained_features = 0

    feature_counts_total = Counter()
    feature_counts_retained = Counter()

    with open_text(annotation_source) as source, gzip.open(
        annotation_output,
        "wt",
        encoding="utf-8",
    ) as destination:
        destination.write(
            "## verified chromosome-filtered annotation\n"
        )
        destination.write(
            f"## species={code}\n"
        )
        destination.write(
            "## sequence IDs renamed using exact sequence-MD5 mapping\n"
        )

        for line in source:
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split(
                "\t"
            )

            if len(fields) < 9:
                continue

            total_features += 1

            original_seqid = fields[0]
            feature = fields[2].lower()
            feature_counts_total[feature] += 1

            chromosome_id = seqid_map.get(
                original_seqid
            )

            if chromosome_id is None:
                continue

            fields[0] = chromosome_id

            attributes = parse_attributes(
                fields[8]
            )

            transcript_id = (
                transcript_from_attributes(
                    feature,
                    attributes,
                )
            )

            gene_id = gene_from_attributes(
                attributes
            )

            if transcript_id:
                retained_transcripts.add(
                    transcript_id
                )

            if gene_id:
                retained_genes.add(gene_id)

            retained_seqids.add(
                chromosome_id
            )

            retained_features += 1
            feature_counts_retained[
                feature
            ] += 1

            destination.write(
                "\t".join(fields) + "\n"
            )

    protein_records = list(
        fasta_records(protein_source)
    )

    protein_index = {}

    for header, sequence in protein_records:
        protein_id = primary_id(header)

        for variant in identifier_variants(
            protein_id
        ):
            protein_index.setdefault(
                variant,
                (header, sequence),
            )

    selected_proteins = {}
    unmatched_transcripts = []

    for transcript_id in sorted(
        retained_transcripts
    ):
        match = None

        for variant in identifier_variants(
            transcript_id
        ):
            if variant in protein_index:
                match = protein_index[variant]
                break

        if match is None:
            unmatched_transcripts.append(
                transcript_id
            )
            continue

        header, sequence = match
        protein_id = primary_id(header)

        selected_proteins[
            protein_id
        ] = (header, sequence)

    protein_output = (
        protein_dir
        / f"{code}.chromosomes.verified.aa.fa"
    )

    with protein_output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for protein_id in sorted(
            selected_proteins
        ):
            header, sequence = (
                selected_proteins[protein_id]
            )

            handle.write(f">{header}\n")

            for start in range(
                0,
                len(sequence),
                60,
            ):
                handle.write(
                    sequence[
                        start:start + 60
                    ]
                    + "\n"
                )

    unmatched_output = (
        table_dir
        / f"{code}.unmatched_transcripts.txt"
    )

    unmatched_output.write_text(
        "\n".join(
            unmatched_transcripts
        )
        + (
            "\n"
            if unmatched_transcripts
            else ""
        ),
        encoding="utf-8",
    )

    missing_annotation_chromosomes = (
        genome_header_set
        - retained_seqids
    )

    unexpected_annotation_seqids = (
        retained_seqids
        - genome_header_set
    )

    transcript_protein_match_pct = (
        100.0
        * len(selected_proteins)
        / len(retained_transcripts)
        if retained_transcripts
        else 0.0
    )

    feature_retention_pct = (
        100.0
        * retained_features
        / total_features
        if total_features
        else 0.0
    )

    warnings = []

    expected_count = (
        expected_chromosome_counts[code]
    )

    if len(genome_headers) != expected_count:
        warnings.append(
            "CHROMOSOME_COUNT_DIFFERS_FROM_EXPECTED"
        )

    if missing_annotation_chromosomes:
        warnings.append(
            "GENOME_CHROMOSOME_WITHOUT_ANNOTATION"
        )

    if unexpected_annotation_seqids:
        warnings.append(
            "ANNOTATION_SEQID_NOT_IN_GENOME"
        )

    if retained_features == 0:
        warnings.append(
            "NO_RETAINED_ANNOTATION_FEATURES"
        )

    if len(selected_proteins) == 0:
        warnings.append(
            "NO_RETAINED_PROTEINS"
        )

    if transcript_protein_match_pct < 99.0:
        warnings.append(
            "TRANSCRIPT_PROTEIN_MATCH_BELOW_99_PERCENT"
        )

    status = (
        "PASS"
        if not warnings
        else "REVIEW"
    )

    for chromosome in genome_headers:
        seqid_validation_rows.append(
            {
                "species_code": code,
                "chromosome_id": chromosome,
                "present_in_filtered_annotation": (
                    chromosome in retained_seqids
                ),
            }
        )

    summary_rows.append(
        {
            "species_code": code,
            "expected_chromosome_count": expected_count,
            "observed_chromosome_count": len(
                genome_headers
            ),
            "mapped_braker_seqids": len(
                seqid_map
            ),
            "annotation_seqids_retained": len(
                retained_seqids
            ),
            "genome_chromosomes_without_annotation": len(
                missing_annotation_chromosomes
            ),
            "unexpected_filtered_annotation_seqids": len(
                unexpected_annotation_seqids
            ),
            "features_total": total_features,
            "features_retained": retained_features,
            "feature_retention_pct": (
                f"{feature_retention_pct:.2f}"
            ),
            "genes_retained": len(
                retained_genes
            ),
            "transcripts_retained": len(
                retained_transcripts
            ),
            "proteins_retained": len(
                selected_proteins
            ),
            "unmatched_transcripts": len(
                unmatched_transcripts
            ),
            "transcript_protein_match_pct": (
                f"{transcript_protein_match_pct:.2f}"
            ),
            "gene_features_retained": (
                feature_counts_retained["gene"]
            ),
            "mrna_features_retained": (
                feature_counts_retained["mrna"]
                + feature_counts_retained[
                    "transcript"
                ]
            ),
            "cds_features_retained": (
                feature_counts_retained["cds"]
            ),
            "status": status,
            "warnings": ";".join(warnings),
        }
    )

    manifest_rows.append(
        {
            "species_code": code,
            "scientific_name": row[
                "scientific_name"
            ],
            "ploidy": row["ploidy"],
            "analysis_group": row[
                "analysis_group"
            ],
            "genome_fasta": str(
                genome_output
            ),
            "annotation_gtf": str(
                annotation_output
            ),
            "protein_fasta": str(
                protein_output
            ),
            "chromosome_count": len(
                genome_headers
            ),
            "protein_count": len(
                selected_proteins
            ),
            "status": status,
            "approved": (
                "YES"
                if status == "PASS"
                else "NO"
            ),
            "notes": ";".join(warnings),
        }
    )

summary_fields = [
    "species_code",
    "expected_chromosome_count",
    "observed_chromosome_count",
    "mapped_braker_seqids",
    "annotation_seqids_retained",
    "genome_chromosomes_without_annotation",
    "unexpected_filtered_annotation_seqids",
    "features_total",
    "features_retained",
    "feature_retention_pct",
    "genes_retained",
    "transcripts_retained",
    "proteins_retained",
    "unmatched_transcripts",
    "transcript_protein_match_pct",
    "gene_features_retained",
    "mrna_features_retained",
    "cds_features_retained",
    "status",
    "warnings",
]

summary_output = (
    table_dir
    / "verified_synteny_input_qc.tsv"
)

with summary_output.open(
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

seqid_output = (
    table_dir
    / "verified_annotation_chromosome_overlap.tsv"
)

with seqid_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "chromosome_id",
            "present_in_filtered_annotation",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(
        seqid_validation_rows
    )

manifest_output = (
    manifest_dir
    / "synteny_input_manifest.verified.tsv"
)

with manifest_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "scientific_name",
            "ploidy",
            "analysis_group",
            "genome_fasta",
            "annotation_gtf",
            "protein_fasta",
            "chromosome_count",
            "protein_count",
            "status",
            "approved",
            "notes",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(manifest_rows)

pass_count = sum(
    row["status"] == "PASS"
    for row in summary_rows
)

review_count = len(summary_rows) - pass_count

overall_output = (
    table_dir
    / "verified_synteny_input_summary.tsv"
)

with overall_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
    )

    writer.writerow(["metric", "value"])
    writer.writerow(
        ["species_total", len(summary_rows)]
    )
    writer.writerow(
        ["species_pass", pass_count]
    )
    writer.writerow(
        ["species_review", review_count]
    )
    writer.writerow(
        [
            "total_proteins_retained",
            sum(
                row["proteins_retained"]
                for row in summary_rows
            ),
        ]
    )
    writer.writerow(
        [
            "total_unmatched_transcripts",
            sum(
                row["unmatched_transcripts"]
                for row in summary_rows
            ),
        ]
    )

print("Verified synteny inputs rebuilt.")
print("Species PASS:", pass_count)
print("Species REVIEW:", review_count)

for row in summary_rows:
    print(
        row["species_code"],
        "chromosomes=",
        row["observed_chromosome_count"],
        "features=",
        row["features_retained"],
        "transcripts=",
        row["transcripts_retained"],
        "proteins=",
        row["proteins_retained"],
        "match_pct=",
        row["transcript_protein_match_pct"],
        "status=",
        row["status"],
    )
PY

echo
echo "============================================================"
echo "Verified input summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/verified_synteny_input_summary.tsv"

echo
echo "============================================================"
echo "Verified input QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/verified_synteny_input_qc.tsv"

echo
echo "============================================================"
echo "Final verified manifest"
echo "============================================================"

column -t -s $'\t' \
    "${MANIFEST_DIR}/synteny_input_manifest.verified.tsv"

cat > "${CHECKPOINT_DIR}/VERIFIED_SYNTENY_INPUTS_COMPLETE.txt" <<EOF2
checkpoint=verified_synteny_inputs
date=$(date --iso-8601=seconds)
status=PASS
next_step=convert_verified_gtf_to_bed_and_check_gene_order
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
echo "Step 32B3 completed successfully."
