#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=jcvi_inputs
#SBATCH --output=10_synteny/logs/jcvi_inputs_%j.out
#SBATCH --error=10_synteny/logs/jcvi_inputs_%j.err

set -euo pipefail

# ============================================================
# Paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

VERIFIED_MANIFEST="${SYNTENY_DIR}/manifests/synteny_input_manifest.verified.tsv"
VERIFIED_CHECKPOINT="${SYNTENY_DIR}/checkpoint_verified_inputs/VERIFIED_SYNTENY_INPUTS_COMPLETE.txt"

OUTPUT_DIR="${SYNTENY_DIR}/jcvi_inputs"
BED_DIR="${OUTPUT_DIR}/bed"
PEP_DIR="${OUTPUT_DIR}/pep"
GENOME_DIR="${OUTPUT_DIR}/genomes"
IDMAP_DIR="${OUTPUT_DIR}/id_maps"
TABLE_DIR="${OUTPUT_DIR}/tables"

MANIFEST_DIR="${SYNTENY_DIR}/manifests"
LOG_DIR="${SYNTENY_DIR}/logs"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_jcvi_inputs"

mkdir -p \
    "${BED_DIR}" \
    "${PEP_DIR}" \
    "${GENOME_DIR}" \
    "${IDMAP_DIR}" \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

cd "${PROJECT_DIR}"

# ============================================================
# Validate prior inputs
# ============================================================

if [[ ! -s "${VERIFIED_MANIFEST}" ]]; then
    echo "ERROR: Verified synteny manifest is missing:" >&2
    echo "${VERIFIED_MANIFEST}" >&2
    exit 1
fi

if [[ ! -s "${VERIFIED_CHECKPOINT}" ]]; then
    echo "ERROR: Step 32B3 checkpoint is missing:" >&2
    echo "${VERIFIED_CHECKPOINT}" >&2
    exit 1
fi

if ! grep -q '^status=PASS$' "${VERIFIED_CHECKPOINT}"; then
    echo "ERROR: Step 32B3 checkpoint did not report PASS." >&2
    cat "${VERIFIED_CHECKPOINT}" >&2
    exit 1
fi

# All nine species should now be approved.
UNAPPROVED_COUNT=$(
    awk -F'\t' '
        NR > 1 && $11 != "YES" {
            count++
        }
        END {
            print count + 0
        }
    ' "${VERIFIED_MANIFEST}"
)

if [[ "${UNAPPROVED_COUNT}" -ne 0 ]]; then
    echo "ERROR: The verified manifest contains unapproved species." >&2
    column -t -s $'\t' "${VERIFIED_MANIFEST}" >&2
    exit 1
fi

# ============================================================
# Activate jcvi_env without sourcing ~/.bashrc
# ============================================================

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
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate jcvi_env

echo "Python:"
python --version

echo "JCVI:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "version attribute unavailable"))
PY

# Confirm the main MCScan module is importable.
python - <<'PY'
from jcvi.compara import catalog
from jcvi.formats import bed
print("JCVI compara and BED modules imported successfully.")
PY

# ============================================================
# Clean only Step 33 products
# ============================================================

rm -f "${BED_DIR}"/*
rm -f "${PEP_DIR}"/*
rm -f "${GENOME_DIR}"/*
rm -f "${IDMAP_DIR}"/*
rm -f "${TABLE_DIR}"/*
rm -f "${CHECKPOINT_DIR}"/*

# ============================================================
# Build BED and protein files
# ============================================================

python - \
    "${VERIFIED_MANIFEST}" \
    "${BED_DIR}" \
    "${PEP_DIR}" \
    "${GENOME_DIR}" \
    "${IDMAP_DIR}" \
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
    manifest_name,
    bed_dir_name,
    pep_dir_name,
    genome_dir_name,
    idmap_dir_name,
    table_dir_name,
    manifest_dir_name,
) = sys.argv[1:]

manifest_file = Path(manifest_name)
bed_dir = Path(bed_dir_name)
pep_dir = Path(pep_dir_name)
genome_dir = Path(genome_dir_name)
idmap_dir = Path(idmap_dir_name)
table_dir = Path(table_dir_name)
manifest_dir = Path(manifest_dir_name)

for directory in [
    bed_dir,
    pep_dir,
    genome_dir,
    idmap_dir,
    table_dir,
    manifest_dir,
]:
    directory.mkdir(parents=True, exist_ok=True)


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
                sequence_parts.append(
                    re.sub(r"\s+", "", line)
                )

    if header is not None:
        yield header, "".join(sequence_parts)


def primary_id(header: str) -> str:
    return header.split()[0]


def parse_attributes(text: str) -> dict[str, str]:
    attributes = {}

    # GTF-style attributes.
    for match in re.finditer(
        r'([A-Za-z0-9_.:-]+)\s+"([^"]*)"',
        text,
    ):
        attributes[match.group(1)] = match.group(2)

    # GFF-style key=value attributes.
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


def normalize_id(value: str) -> str:
    value = value.strip()

    value = re.sub(
        r"^(gene:|transcript:|protein:)",
        "",
        value,
    )

    return value


def transcript_id_from_feature(
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

    for candidate in candidates:
        if not candidate:
            continue

        candidate = candidate.split(",")[0]
        return normalize_id(candidate)

    return ""


def gene_id_from_feature(
    attributes: dict[str, str],
) -> str:
    for key in [
        "gene_id",
        "gene",
        "gene_name",
        "Parent",
    ]:
        value = attributes.get(key)

        if value:
            return normalize_id(
                value.split(",")[0]
            )

    return ""


def identifier_variants(identifier: str) -> set[str]:
    value = normalize_id(identifier)

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

    variants.add(
        re.sub(r"\.protein\d+$", "", value)
    )

    variants.add(
        re.sub(r"\.t\d+$", "", value)
    )

    return {
        variant
        for variant in variants
        if variant
    }


def prefixed_id(code: str, identifier: str) -> str:
    identifier = normalize_id(identifier)

    if identifier.startswith(code + "_"):
        return identifier

    return f"{code}_{identifier}"


def chromosome_lengths(path: Path) -> dict[str, int]:
    lengths = {}

    for header, sequence in fasta_records(path):
        seqid = primary_id(header)

        if seqid in lengths:
            raise ValueError(
                f"Duplicate chromosome ID in {path}: {seqid}"
            )

        lengths[seqid] = len(sequence)

    return lengths


with manifest_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    manifest_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

summary_rows = []
final_manifest_rows = []
all_coordinate_issues = []
all_duplicate_issues = []

for manifest_row in manifest_rows:
    code = manifest_row["species_code"]

    genome_source = Path(
        manifest_row["genome_fasta"]
    )
    annotation_source = Path(
        manifest_row["annotation_gtf"]
    )
    protein_source = Path(
        manifest_row["protein_fasta"]
    )

    for path in [
        genome_source,
        annotation_source,
        protein_source,
    ]:
        if not path.is_file():
            raise SystemExit(
                f"ERROR: Missing input for {code}: {path}"
            )

    chromosomes = chromosome_lengths(
        genome_source
    )

    # --------------------------------------------------------
    # Collect transcript models
    # --------------------------------------------------------

    transcript_models = {}
    cds_spans = defaultdict(list)
    transcript_to_gene = {}

    annotation_feature_counts = Counter()
    malformed_lines = 0
    unknown_seqid_lines = 0

    with open_text(annotation_source) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                malformed_lines += 1
                continue

            seqid = fields[0]
            feature = fields[2].lower()

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                malformed_lines += 1
                continue

            strand = fields[6]

            annotation_feature_counts[feature] += 1

            if seqid not in chromosomes:
                unknown_seqid_lines += 1
                continue

            attributes = parse_attributes(
                fields[8]
            )

            transcript_id = transcript_id_from_feature(
                feature,
                attributes,
            )

            gene_id = gene_id_from_feature(
                attributes
            )

            if feature in {"mrna", "transcript"}:
                if not transcript_id:
                    continue

                transcript_models[transcript_id] = {
                    "seqid": seqid,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "gene_id": gene_id or transcript_id,
                }

                transcript_to_gene[transcript_id] = (
                    gene_id or transcript_id
                )

            elif feature == "cds":
                if not transcript_id:
                    continue

                cds_spans[transcript_id].append(
                    (
                        seqid,
                        start,
                        end,
                        strand,
                    )
                )

                if gene_id:
                    transcript_to_gene.setdefault(
                        transcript_id,
                        gene_id,
                    )

    # Reconstruct transcript coordinates from CDS if mRNA
    # coordinates were absent.
    for transcript_id, spans in cds_spans.items():
        if transcript_id in transcript_models:
            continue

        seqids = {
            span[0]
            for span in spans
        }

        strands = {
            span[3]
            for span in spans
        }

        if len(seqids) != 1 or len(strands) != 1:
            continue

        transcript_models[transcript_id] = {
            "seqid": next(iter(seqids)),
            "start": min(
                span[1]
                for span in spans
            ),
            "end": max(
                span[2]
                for span in spans
            ),
            "strand": next(iter(strands)),
            "gene_id": transcript_to_gene.get(
                transcript_id,
                transcript_id,
            ),
        }

    # --------------------------------------------------------
    # Choose one representative transcript per gene
    # --------------------------------------------------------

    transcripts_by_gene = defaultdict(list)

    for transcript_id, model in transcript_models.items():
        gene_id = model["gene_id"] or transcript_id

        cds_length = sum(
            end - start + 1
            for _, start, end, _
            in cds_spans.get(
                transcript_id,
                []
            )
        )

        genomic_span = (
            model["end"]
            - model["start"]
            + 1
        )

        transcripts_by_gene[gene_id].append(
            (
                transcript_id,
                cds_length,
                genomic_span,
            )
        )

    representative_transcripts = {}

    for gene_id, candidates in transcripts_by_gene.items():
        # Prefer longest CDS, then longest genomic span,
        # then lexical transcript ID for reproducibility.
        selected = sorted(
            candidates,
            key=lambda item: (
                -item[1],
                -item[2],
                item[0],
            ),
        )[0]

        representative_transcripts[
            selected[0]
        ] = gene_id

    # Because the input is already braker.longest.gtf, normally
    # one transcript should exist per gene. This additional check
    # protects against unexpected duplicate models.
    representative_count = len(
        representative_transcripts
    )

    # --------------------------------------------------------
    # Index proteins
    # --------------------------------------------------------

    protein_records = list(
        fasta_records(protein_source)
    )

    protein_index = {}
    duplicate_protein_variants = set()

    for header, sequence in protein_records:
        protein_identifier = primary_id(header)

        for variant in identifier_variants(
            protein_identifier
        ):
            if (
                variant in protein_index
                and protein_index[variant][0]
                != protein_identifier
            ):
                duplicate_protein_variants.add(
                    variant
                )
                continue

            protein_index[variant] = (
                protein_identifier,
                header,
                sequence,
            )

    # --------------------------------------------------------
    # Match representative transcripts to proteins
    # --------------------------------------------------------

    selected_rows = []
    selected_proteins = {}
    unmatched_transcripts = []
    ambiguous_output_ids = Counter()

    for transcript_id, gene_id in representative_transcripts.items():
        model = transcript_models[
            transcript_id
        ]

        protein_match = None

        for variant in identifier_variants(
            transcript_id
        ):
            candidate = protein_index.get(
                variant
            )

            if candidate is not None:
                protein_match = candidate
                break

        if protein_match is None:
            unmatched_transcripts.append(
                transcript_id
            )
            continue

        (
            original_protein_id,
            original_header,
            protein_sequence,
        ) = protein_match

        output_id = prefixed_id(
            code,
            transcript_id,
        )

        ambiguous_output_ids[
            output_id
        ] += 1

        selected_rows.append(
            {
                "seqid": model["seqid"],
                # GTF is 1-based closed; BED is 0-based half-open.
                "start0": model["start"] - 1,
                "end": model["end"],
                "id": output_id,
                "score": "0",
                "strand": (
                    model["strand"]
                    if model["strand"] in {"+", "-"}
                    else "."
                ),
                "original_transcript_id": transcript_id,
                "original_gene_id": gene_id,
                "original_protein_id": original_protein_id,
            }
        )

        selected_proteins[
            output_id
        ] = protein_sequence

    duplicate_output_ids = {
        identifier
        for identifier, count
        in ambiguous_output_ids.items()
        if count > 1
    }

    # --------------------------------------------------------
    # Coordinate validation
    # --------------------------------------------------------

    valid_rows = []

    for bed_row in selected_rows:
        seqid = bed_row["seqid"]
        start0 = bed_row["start0"]
        end = bed_row["end"]
        chromosome_length = chromosomes[
            seqid
        ]

        issue = ""

        if start0 < 0:
            issue = "NEGATIVE_BED_START"
        elif end <= start0:
            issue = "END_NOT_GREATER_THAN_START"
        elif end > chromosome_length:
            issue = "END_EXCEEDS_CHROMOSOME_LENGTH"
        elif bed_row["id"] in duplicate_output_ids:
            issue = "DUPLICATE_OUTPUT_ID"

        if issue:
            all_coordinate_issues.append(
                {
                    "species_code": code,
                    "gene_id": bed_row["id"],
                    "chromosome": seqid,
                    "start0": start0,
                    "end": end,
                    "chromosome_length": chromosome_length,
                    "issue": issue,
                }
            )
            continue

        valid_rows.append(bed_row)

    # --------------------------------------------------------
    # Sort in genome order
    # --------------------------------------------------------

    chromosome_order = {
        chromosome: index
        for index, chromosome
        in enumerate(chromosomes)
    }

    valid_rows.sort(
        key=lambda row: (
            chromosome_order[row["seqid"]],
            row["start0"],
            row["end"],
            row["id"],
        )
    )

    # --------------------------------------------------------
    # Ensure one protein per retained BED ID
    # --------------------------------------------------------

    bed_ids = {
        row["id"]
        for row in valid_rows
    }

    selected_proteins = {
        identifier: sequence
        for identifier, sequence
        in selected_proteins.items()
        if identifier in bed_ids
    }

    protein_ids = set(
        selected_proteins
    )

    bed_without_protein = sorted(
        bed_ids - protein_ids
    )

    protein_without_bed = sorted(
        protein_ids - bed_ids
    )

    # --------------------------------------------------------
    # Write JCVI files
    # --------------------------------------------------------

    bed_output = (
        bed_dir / f"{code}.bed"
    )

    pep_output = (
        pep_dir / f"{code}.pep"
    )

    genome_output = (
        genome_dir / f"{code}.fa.gz"
    )

    idmap_output = (
        idmap_dir / f"{code}.id_map.tsv"
    )

    chromosome_order_output = (
        table_dir / f"{code}.chromosome_order.tsv"
    )

    with bed_output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for bed_row in valid_rows:
            handle.write(
                "\t".join(
                    [
                        bed_row["seqid"],
                        str(bed_row["start0"]),
                        str(bed_row["end"]),
                        bed_row["id"],
                        bed_row["score"],
                        bed_row["strand"],
                    ]
                )
                + "\n"
            )

    with pep_output.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for identifier in sorted(
            selected_proteins
        ):
            sequence = selected_proteins[
                identifier
            ]

            handle.write(
                f">{identifier}\n"
            )

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

    with open_text(genome_source) as source, gzip.open(
        genome_output,
        "wt",
        encoding="utf-8",
    ) as destination:
        shutil.copyfileobj(
            source,
            destination,
        )

    with idmap_output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "species_code",
                "jcvi_id",
                "original_transcript_id",
                "original_gene_id",
                "original_protein_id",
                "chromosome",
                "bed_start0",
                "bed_end",
                "strand",
            ]
        )

        for bed_row in valid_rows:
            writer.writerow(
                [
                    code,
                    bed_row["id"],
                    bed_row[
                        "original_transcript_id"
                    ],
                    bed_row[
                        "original_gene_id"
                    ],
                    bed_row[
                        "original_protein_id"
                    ],
                    bed_row["seqid"],
                    bed_row["start0"],
                    bed_row["end"],
                    bed_row["strand"],
                ]
            )

    with chromosome_order_output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
        )

        writer.writerow(
            [
                "species_code",
                "chromosome_order",
                "chromosome_id",
                "chromosome_length_bp",
                "genes_in_bed",
            ]
        )

        genes_per_chromosome = Counter(
            row["seqid"]
            for row in valid_rows
        )

        for order, chromosome in enumerate(
            chromosomes,
            start=1,
        ):
            writer.writerow(
                [
                    code,
                    order,
                    chromosome,
                    chromosomes[
                        chromosome
                    ],
                    genes_per_chromosome[
                        chromosome
                    ],
                ]
            )

    unmatched_output = (
        table_dir
        / f"{code}.unmatched_representative_transcripts.txt"
    )

    unmatched_output.write_text(
        "\n".join(
            sorted(unmatched_transcripts)
        )
        + (
            "\n"
            if unmatched_transcripts
            else ""
        ),
        encoding="utf-8",
    )

    # --------------------------------------------------------
    # Species QC
    # --------------------------------------------------------

    match_percentage = (
        100.0
        * len(bed_ids)
        / representative_count
        if representative_count
        else 0.0
    )

    chromosomes_without_genes = [
        chromosome
        for chromosome in chromosomes
        if not any(
            row["seqid"] == chromosome
            for row in valid_rows
        )
    ]

    warnings = []

    if malformed_lines > 0:
        warnings.append(
            "MALFORMED_GTF_LINES"
        )

    if unknown_seqid_lines > 0:
        warnings.append(
            "GTF_SEQIDS_NOT_IN_GENOME"
        )

    if representative_count == 0:
        warnings.append(
            "NO_REPRESENTATIVE_TRANSCRIPTS"
        )

    if len(valid_rows) == 0:
        warnings.append(
            "EMPTY_BED"
        )

    if len(selected_proteins) == 0:
        warnings.append(
            "EMPTY_PEP"
        )

    if duplicate_output_ids:
        warnings.append(
            "DUPLICATE_JCVI_IDS"
        )

    if duplicate_protein_variants:
        warnings.append(
            "AMBIGUOUS_PROTEIN_ID_VARIANTS"
        )

    if bed_without_protein:
        warnings.append(
            "BED_IDS_WITHOUT_PROTEINS"
        )

    if protein_without_bed:
        warnings.append(
            "PROTEINS_WITHOUT_BED_IDS"
        )

    if chromosomes_without_genes:
        warnings.append(
            "CHROMOSOMES_WITHOUT_GENES"
        )

    if match_percentage < 99.0:
        warnings.append(
            "REPRESENTATIVE_TRANSCRIPT_MATCH_BELOW_99_PERCENT"
        )

    status = (
        "PASS"
        if not warnings
        else "REVIEW"
    )

    summary_rows.append(
        {
            "species_code": code,
            "chromosome_count": len(
                chromosomes
            ),
            "annotation_gene_features": (
                annotation_feature_counts["gene"]
            ),
            "annotation_mrna_features": (
                annotation_feature_counts["mrna"]
                + annotation_feature_counts[
                    "transcript"
                ]
            ),
            "annotation_cds_features": (
                annotation_feature_counts["cds"]
            ),
            "representative_transcripts": (
                representative_count
            ),
            "bed_gene_count": len(
                valid_rows
            ),
            "protein_count": len(
                selected_proteins
            ),
            "unmatched_representative_transcripts": len(
                unmatched_transcripts
            ),
            "bed_ids_without_proteins": len(
                bed_without_protein
            ),
            "protein_ids_without_bed": len(
                protein_without_bed
            ),
            "duplicate_jcvi_ids": len(
                duplicate_output_ids
            ),
            "coordinate_issues": sum(
                1
                for issue in all_coordinate_issues
                if issue[
                    "species_code"
                ] == code
            ),
            "chromosomes_without_genes": len(
                chromosomes_without_genes
            ),
            "representative_transcript_match_pct": (
                f"{match_percentage:.2f}"
            ),
            "status": status,
            "warnings": ";".join(
                warnings
            ),
        }
    )

    final_manifest_rows.append(
        {
            "species_code": code,
            "scientific_name": manifest_row[
                "scientific_name"
            ],
            "ploidy": manifest_row[
                "ploidy"
            ],
            "analysis_group": manifest_row[
                "analysis_group"
            ],
            "genome_fasta": str(
                genome_output
            ),
            "bed_file": str(
                bed_output
            ),
            "protein_fasta": str(
                pep_output
            ),
            "id_map": str(
                idmap_output
            ),
            "chromosome_count": len(
                chromosomes
            ),
            "bed_gene_count": len(
                valid_rows
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
            "notes": ";".join(
                warnings
            ),
        }
    )

# ============================================================
# Write global QC outputs
# ============================================================

summary_fields = [
    "species_code",
    "chromosome_count",
    "annotation_gene_features",
    "annotation_mrna_features",
    "annotation_cds_features",
    "representative_transcripts",
    "bed_gene_count",
    "protein_count",
    "unmatched_representative_transcripts",
    "bed_ids_without_proteins",
    "protein_ids_without_bed",
    "duplicate_jcvi_ids",
    "coordinate_issues",
    "chromosomes_without_genes",
    "representative_transcript_match_pct",
    "status",
    "warnings",
]

summary_output = (
    table_dir
    / "jcvi_input_qc.tsv"
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
    writer.writerows(
        summary_rows
    )

coordinate_issue_output = (
    table_dir
    / "jcvi_coordinate_issues.tsv"
)

with coordinate_issue_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    fields = [
        "species_code",
        "gene_id",
        "chromosome",
        "start0",
        "end",
        "chromosome_length",
        "issue",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(
        all_coordinate_issues
    )

manifest_output = (
    manifest_dir
    / "jcvi_input_manifest.tsv"
)

manifest_fields = [
    "species_code",
    "scientific_name",
    "ploidy",
    "analysis_group",
    "genome_fasta",
    "bed_file",
    "protein_fasta",
    "id_map",
    "chromosome_count",
    "bed_gene_count",
    "protein_count",
    "status",
    "approved",
    "notes",
]

with manifest_output.open(
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
    writer.writerows(
        final_manifest_rows
    )

pass_count = sum(
    row["status"] == "PASS"
    for row in summary_rows
)

review_count = len(summary_rows) - pass_count

overall_output = (
    table_dir
    / "jcvi_input_summary.tsv"
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

    writer.writerow(
        ["metric", "value"]
    )

    writer.writerow(
        [
            "species_total",
            len(summary_rows),
        ]
    )

    writer.writerow(
        [
            "species_pass",
            pass_count,
        ]
    )

    writer.writerow(
        [
            "species_review",
            review_count,
        ]
    )

    writer.writerow(
        [
            "total_bed_genes",
            sum(
                row["bed_gene_count"]
                for row in summary_rows
            ),
        ]
    )

    writer.writerow(
        [
            "total_proteins",
            sum(
                row["protein_count"]
                for row in summary_rows
            ),
        ]
    )

    writer.writerow(
        [
            "total_unmatched_transcripts",
            sum(
                row[
                    "unmatched_representative_transcripts"
                ]
                for row in summary_rows
            ),
        ]
    )

    writer.writerow(
        [
            "total_coordinate_issues",
            len(
                all_coordinate_issues
            ),
        ]
    )

print("JCVI inputs prepared.")

for row in summary_rows:
    print(
        row["species_code"],
        "BED=",
        row["bed_gene_count"],
        "PEP=",
        row["protein_count"],
        "unmatched=",
        row[
            "unmatched_representative_transcripts"
        ],
        "coordinate_issues=",
        row["coordinate_issues"],
        "status=",
        row["status"],
    )
PY

# ============================================================
# Independent JCVI BED validation
# ============================================================

echo
echo "============================================================"
echo "Independent JCVI BED checks"
echo "============================================================"

for BED_FILE in "${BED_DIR}"/*.bed
do
    CODE=$(basename "${BED_FILE}" .bed)

    echo "Checking ${CODE}"

    python -m jcvi.formats.bed uniq \
        "${BED_FILE}" \
        --outfile="${TABLE_DIR}/${CODE}.uniq.bed" \
        >/dev/null

    ORIGINAL_COUNT=$(wc -l < "${BED_FILE}")
    UNIQUE_COUNT=$(wc -l < "${TABLE_DIR}/${CODE}.uniq.bed")

    if [[ "${ORIGINAL_COUNT}" -ne "${UNIQUE_COUNT}" ]]; then
        echo "ERROR: JCVI detected duplicate BED IDs in ${CODE}." >&2
        echo "Original=${ORIGINAL_COUNT}; unique=${UNIQUE_COUNT}" >&2
        exit 1
    fi
done

# ============================================================
# Strict final checks
# ============================================================

if awk -F'\t' '
    NR > 1 && $16 != "PASS" {
        failed = 1
    }
    END {
        exit failed
    }
' "${TABLE_DIR}/jcvi_input_qc.tsv"
then
    :
else
    echo "ERROR: At least one species failed Step 33 QC." >&2
    column -t -s $'\t' \
        "${TABLE_DIR}/jcvi_input_qc.tsv" >&2
    exit 1
fi

if awk -F'\t' '
    NR > 1 && $13 != "YES" {
        failed = 1
    }
    END {
        exit failed
    }
' "${MANIFEST_DIR}/jcvi_input_manifest.tsv"
then
    :
else
    echo "ERROR: At least one JCVI input is not approved." >&2
    column -t -s $'\t' \
        "${MANIFEST_DIR}/jcvi_input_manifest.tsv" >&2
    exit 1
fi

# ============================================================
# Display results
# ============================================================

echo
echo "============================================================"
echo "JCVI input summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/jcvi_input_summary.tsv"

echo
echo "============================================================"
echo "JCVI input QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/jcvi_input_qc.tsv"

echo
echo "============================================================"
echo "Final JCVI manifest"
echo "============================================================"

column -t -s $'\t' \
    "${MANIFEST_DIR}/jcvi_input_manifest.tsv"

# ============================================================
# Checkpoint
# ============================================================

cp -f \
    "${TABLE_DIR}/jcvi_input_summary.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_input_qc.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_coordinate_issues.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${MANIFEST_DIR}/jcvi_input_manifest.tsv" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/JCVI_INPUTS_COMPLETE.txt" <<EOF2
checkpoint=jcvi_bed_protein_inputs
date=$(date --iso-8601=seconds)
conda_environment=jcvi_env
species_total=9
status=PASS
next_step=run_pairwise_jcvi_ortholog_and_self_synteny_analyses
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
echo "============================================================"
echo "Step 33 completed successfully"
echo "============================================================"
echo "Checkpoint:"
echo "${CHECKPOINT_DIR}"
