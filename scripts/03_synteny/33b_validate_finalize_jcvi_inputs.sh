#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=jcvi_validate
#SBATCH --output=10_synteny/logs/jcvi_validate_%j.out
#SBATCH --error=10_synteny/logs/jcvi_validate_%j.err

set -euo pipefail

# ============================================================
# Project paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

VERIFIED_MANIFEST="${SYNTENY_DIR}/manifests/synteny_input_manifest.verified.tsv"

INPUT_DIR="${SYNTENY_DIR}/jcvi_inputs"
BED_DIR="${INPUT_DIR}/bed"
PEP_DIR="${INPUT_DIR}/pep"
GENOME_DIR="${INPUT_DIR}/genomes"
IDMAP_DIR="${INPUT_DIR}/id_maps"
TABLE_DIR="${INPUT_DIR}/tables"

MANIFEST_DIR="${SYNTENY_DIR}/manifests"
OUTPUT_MANIFEST="${MANIFEST_DIR}/jcvi_input_manifest.tsv"

CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_jcvi_inputs"
LOG_DIR="${SYNTENY_DIR}/logs"

mkdir -p \
    "${TABLE_DIR}" \
    "${MANIFEST_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

cd "${PROJECT_DIR}"

# ============================================================
# Check required inputs
# ============================================================

if [[ ! -s "${VERIFIED_MANIFEST}" ]]; then
    echo "ERROR: Missing verified manifest:" >&2
    echo "${VERIFIED_MANIFEST}" >&2
    exit 1
fi

for CODE in VPAN VSCU VANA VARV VPER VSER VTRI VVER PMAJ
do
    for FILE in \
        "${BED_DIR}/${CODE}.bed" \
        "${PEP_DIR}/${CODE}.pep" \
        "${GENOME_DIR}/${CODE}.fa.gz" \
        "${IDMAP_DIR}/${CODE}.id_map.tsv"
    do
        if [[ ! -s "${FILE}" ]]; then
            echo "ERROR: Missing or empty Step 33 product:" >&2
            echo "${FILE}" >&2
            exit 1
        fi
    done
done

# ============================================================
# Activate jcvi_env
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

echo "Python version:"
python --version

echo "JCVI version:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY

python - <<'PY'
from jcvi.compara import catalog
from jcvi.formats import bed
print("JCVI modules imported successfully.")
PY

# ============================================================
# Remove only old validation products
# ============================================================

rm -f \
    "${TABLE_DIR}/jcvi_input_qc.tsv" \
    "${TABLE_DIR}/jcvi_input_summary.tsv" \
    "${TABLE_DIR}/jcvi_coordinate_issues.tsv" \
    "${TABLE_DIR}/jcvi_duplicate_ids.tsv" \
    "${TABLE_DIR}/jcvi_id_mismatches.tsv" \
    "${TABLE_DIR}/jcvi_gene_order_issues.tsv" \
    "${OUTPUT_MANIFEST}"

rm -f "${CHECKPOINT_DIR}"/*

# ============================================================
# Independent strict validation
# ============================================================

python - \
    "${VERIFIED_MANIFEST}" \
    "${BED_DIR}" \
    "${PEP_DIR}" \
    "${GENOME_DIR}" \
    "${IDMAP_DIR}" \
    "${TABLE_DIR}" \
    "${OUTPUT_MANIFEST}" <<'PY'
from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

(
    verified_manifest_name,
    bed_dir_name,
    pep_dir_name,
    genome_dir_name,
    idmap_dir_name,
    table_dir_name,
    output_manifest_name,
) = sys.argv[1:]

verified_manifest = Path(verified_manifest_name)
bed_dir = Path(bed_dir_name)
pep_dir = Path(pep_dir_name)
genome_dir = Path(genome_dir_name)
idmap_dir = Path(idmap_dir_name)
table_dir = Path(table_dir_name)
output_manifest = Path(output_manifest_name)

table_dir.mkdir(parents=True, exist_ok=True)
output_manifest.parent.mkdir(parents=True, exist_ok=True)


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


def fasta_ids(path: Path):
    return [
        header.split()[0]
        for header, _ in fasta_records(path)
    ]


def fasta_lengths(path: Path):
    lengths = {}

    for header, sequence in fasta_records(path):
        seqid = header.split()[0]

        if seqid in lengths:
            raise ValueError(
                f"Duplicate genome sequence ID in {path}: {seqid}"
            )

        lengths[seqid] = len(sequence)

    return lengths


with verified_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    verified_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

expected_species = {
    "VPAN",
    "VSCU",
    "VANA",
    "VARV",
    "VPER",
    "VSER",
    "VTRI",
    "VVER",
    "PMAJ",
}

observed_species = {
    row["species_code"]
    for row in verified_rows
}

if observed_species != expected_species:
    missing = sorted(
        expected_species - observed_species
    )
    unexpected = sorted(
        observed_species - expected_species
    )

    raise SystemExit(
        "ERROR: Species mismatch in verified manifest. "
        f"Missing={missing}; unexpected={unexpected}"
    )

qc_rows = []
manifest_rows = []

coordinate_issues = []
duplicate_issues = []
id_mismatch_rows = []
gene_order_issues = []

for metadata in verified_rows:
    code = metadata["species_code"]

    bed_file = bed_dir / f"{code}.bed"
    pep_file = pep_dir / f"{code}.pep"
    genome_file = genome_dir / f"{code}.fa.gz"
    idmap_file = idmap_dir / f"{code}.id_map.tsv"

    chromosome_lengths = fasta_lengths(
        genome_file
    )

    chromosome_order = {
        chromosome: index
        for index, chromosome in enumerate(
            chromosome_lengths
        )
    }

    # --------------------------------------------------------
    # Read BED
    # --------------------------------------------------------

    bed_rows = []
    malformed_bed_lines = 0

    with bed_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 6:
                malformed_bed_lines += 1
                coordinate_issues.append(
                    {
                        "species_code": code,
                        "line_number": line_number,
                        "gene_id": "",
                        "chromosome": (
                            fields[0]
                            if fields
                            else ""
                        ),
                        "start0": "",
                        "end": "",
                        "chromosome_length": "",
                        "issue": (
                            "BED_COLUMN_COUNT_NOT_SIX"
                        ),
                    }
                )
                continue

            chromosome, start_text, end_text, gene_id, score, strand = fields

            try:
                start0 = int(start_text)
                end = int(end_text)
            except ValueError:
                coordinate_issues.append(
                    {
                        "species_code": code,
                        "line_number": line_number,
                        "gene_id": gene_id,
                        "chromosome": chromosome,
                        "start0": start_text,
                        "end": end_text,
                        "chromosome_length": "",
                        "issue": (
                            "NON_INTEGER_BED_COORDINATE"
                        ),
                    }
                )
                continue

            bed_rows.append(
                {
                    "line_number": line_number,
                    "chromosome": chromosome,
                    "start0": start0,
                    "end": end,
                    "gene_id": gene_id,
                    "score": score,
                    "strand": strand,
                }
            )

    bed_ids = [
        row["gene_id"]
        for row in bed_rows
    ]

    bed_id_counts = Counter(
        bed_ids
    )

    duplicate_bed_ids = sorted(
        gene_id
        for gene_id, count in bed_id_counts.items()
        if count > 1
    )

    for gene_id in duplicate_bed_ids:
        duplicate_issues.append(
            {
                "species_code": code,
                "file_type": "BED",
                "identifier": gene_id,
                "occurrences": bed_id_counts[
                    gene_id
                ],
            }
        )

    # --------------------------------------------------------
    # Read protein FASTA
    # --------------------------------------------------------

    protein_ids = fasta_ids(
        pep_file
    )

    protein_id_counts = Counter(
        protein_ids
    )

    duplicate_protein_ids = sorted(
        protein_id
        for protein_id, count
        in protein_id_counts.items()
        if count > 1
    )

    for protein_id in duplicate_protein_ids:
        duplicate_issues.append(
            {
                "species_code": code,
                "file_type": "PEP",
                "identifier": protein_id,
                "occurrences": protein_id_counts[
                    protein_id
                ],
            }
        )

    bed_id_set = set(
        bed_ids
    )

    protein_id_set = set(
        protein_ids
    )

    bed_without_protein = sorted(
        bed_id_set - protein_id_set
    )

    protein_without_bed = sorted(
        protein_id_set - bed_id_set
    )

    for identifier in bed_without_protein:
        id_mismatch_rows.append(
            {
                "species_code": code,
                "identifier": identifier,
                "issue": "BED_ID_WITHOUT_PROTEIN",
            }
        )

    for identifier in protein_without_bed:
        id_mismatch_rows.append(
            {
                "species_code": code,
                "identifier": identifier,
                "issue": "PROTEIN_ID_WITHOUT_BED",
            }
        )

    # --------------------------------------------------------
    # Validate coordinates and species prefixes
    # --------------------------------------------------------

    species_coordinate_issue_count = 0
    invalid_prefix_count = 0

    for row in bed_rows:
        chromosome = row["chromosome"]
        start0 = row["start0"]
        end = row["end"]
        gene_id = row["gene_id"]

        issue = ""

        if chromosome not in chromosome_lengths:
            issue = "BED_CHROMOSOME_NOT_IN_GENOME"

        elif start0 < 0:
            issue = "NEGATIVE_BED_START"

        elif end <= start0:
            issue = "BED_END_NOT_GREATER_THAN_START"

        elif end > chromosome_lengths[
            chromosome
        ]:
            issue = "BED_END_EXCEEDS_CHROMOSOME"

        elif row["strand"] not in {
            "+",
            "-",
            ".",
        }:
            issue = "INVALID_BED_STRAND"

        if issue:
            species_coordinate_issue_count += 1

            coordinate_issues.append(
                {
                    "species_code": code,
                    "line_number": row[
                        "line_number"
                    ],
                    "gene_id": gene_id,
                    "chromosome": chromosome,
                    "start0": start0,
                    "end": end,
                    "chromosome_length": (
                        chromosome_lengths.get(
                            chromosome,
                            "",
                        )
                    ),
                    "issue": issue,
                }
            )

        if not gene_id.startswith(
            code + "_"
        ):
            invalid_prefix_count += 1

            id_mismatch_rows.append(
                {
                    "species_code": code,
                    "identifier": gene_id,
                    "issue": (
                        "BED_ID_MISSING_SPECIES_PREFIX"
                    ),
                }
            )

    invalid_protein_prefix_count = 0

    for protein_id in protein_ids:
        if not protein_id.startswith(
            code + "_"
        ):
            invalid_protein_prefix_count += 1

            id_mismatch_rows.append(
                {
                    "species_code": code,
                    "identifier": protein_id,
                    "issue": (
                        "PROTEIN_ID_MISSING_SPECIES_PREFIX"
                    ),
                }
            )

    # --------------------------------------------------------
    # Validate BED gene order
    # --------------------------------------------------------

    previous_key = None
    species_order_issue_count = 0

    for row in bed_rows:
        chromosome = row["chromosome"]

        if chromosome not in chromosome_order:
            continue

        current_key = (
            chromosome_order[
                chromosome
            ],
            row["start0"],
            row["end"],
            row["gene_id"],
        )

        if (
            previous_key is not None
            and current_key < previous_key
        ):
            species_order_issue_count += 1

            gene_order_issues.append(
                {
                    "species_code": code,
                    "line_number": row[
                        "line_number"
                    ],
                    "gene_id": row[
                        "gene_id"
                    ],
                    "chromosome": chromosome,
                    "start0": row[
                        "start0"
                    ],
                    "issue": (
                        "BED_NOT_SORTED_IN_GENOME_ORDER"
                    ),
                }
            )

        previous_key = current_key

    # --------------------------------------------------------
    # Chromosome representation
    # --------------------------------------------------------

    genes_per_chromosome = Counter(
        row["chromosome"]
        for row in bed_rows
        if row["chromosome"]
        in chromosome_lengths
    )

    chromosomes_without_genes = sorted(
        chromosome
        for chromosome in chromosome_lengths
        if genes_per_chromosome[
            chromosome
        ] == 0
    )

    # --------------------------------------------------------
    # Validate ID map
    # --------------------------------------------------------

    with idmap_file.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        idmap_rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    idmap_ids = {
        row["jcvi_id"]
        for row in idmap_rows
    }

    bed_without_idmap = sorted(
        bed_id_set - idmap_ids
    )

    idmap_without_bed = sorted(
        idmap_ids - bed_id_set
    )

    for identifier in bed_without_idmap:
        id_mismatch_rows.append(
            {
                "species_code": code,
                "identifier": identifier,
                "issue": "BED_ID_WITHOUT_ID_MAP",
            }
        )

    for identifier in idmap_without_bed:
        id_mismatch_rows.append(
            {
                "species_code": code,
                "identifier": identifier,
                "issue": "ID_MAP_ID_WITHOUT_BED",
            }
        )

    # --------------------------------------------------------
    # Species status
    # --------------------------------------------------------

    warnings = []

    if malformed_bed_lines > 0:
        warnings.append(
            "MALFORMED_BED_LINES"
        )

    if duplicate_bed_ids:
        warnings.append(
            "DUPLICATE_BED_IDS"
        )

    if duplicate_protein_ids:
        warnings.append(
            "DUPLICATE_PROTEIN_IDS"
        )

    if bed_without_protein:
        warnings.append(
            "BED_IDS_WITHOUT_PROTEINS"
        )

    if protein_without_bed:
        warnings.append(
            "PROTEIN_IDS_WITHOUT_BED"
        )

    if species_coordinate_issue_count > 0:
        warnings.append(
            "BED_COORDINATE_ISSUES"
        )

    if species_order_issue_count > 0:
        warnings.append(
            "BED_GENE_ORDER_ISSUES"
        )

    if invalid_prefix_count > 0:
        warnings.append(
            "BED_IDS_MISSING_SPECIES_PREFIX"
        )

    if invalid_protein_prefix_count > 0:
        warnings.append(
            "PROTEIN_IDS_MISSING_SPECIES_PREFIX"
        )

    if chromosomes_without_genes:
        warnings.append(
            "CHROMOSOMES_WITHOUT_GENES"
        )

    if bed_without_idmap:
        warnings.append(
            "BED_IDS_WITHOUT_ID_MAP"
        )

    if idmap_without_bed:
        warnings.append(
            "ID_MAP_IDS_WITHOUT_BED"
        )

    if len(bed_rows) == 0:
        warnings.append(
            "EMPTY_BED"
        )

    if len(protein_ids) == 0:
        warnings.append(
            "EMPTY_PROTEIN_FASTA"
        )

    match_percentage = (
        100.0
        * len(
            bed_id_set
            & protein_id_set
        )
        / len(bed_id_set)
        if bed_id_set
        else 0.0
    )

    status = (
        "PASS"
        if not warnings
        else "FAIL"
    )

    qc_rows.append(
        {
            "species_code": code,
            "chromosome_count": len(
                chromosome_lengths
            ),
            "bed_gene_count": len(
                bed_rows
            ),
            "unique_bed_ids": len(
                bed_id_set
            ),
            "protein_count": len(
                protein_ids
            ),
            "unique_protein_ids": len(
                protein_id_set
            ),
            "id_map_count": len(
                idmap_rows
            ),
            "duplicate_bed_ids": len(
                duplicate_bed_ids
            ),
            "duplicate_protein_ids": len(
                duplicate_protein_ids
            ),
            "bed_ids_without_proteins": len(
                bed_without_protein
            ),
            "protein_ids_without_bed": len(
                protein_without_bed
            ),
            "bed_ids_without_id_map": len(
                bed_without_idmap
            ),
            "id_map_ids_without_bed": len(
                idmap_without_bed
            ),
            "coordinate_issues": (
                species_coordinate_issue_count
            ),
            "gene_order_issues": (
                species_order_issue_count
            ),
            "chromosomes_without_genes": len(
                chromosomes_without_genes
            ),
            "bed_ids_missing_species_prefix": (
                invalid_prefix_count
            ),
            "protein_ids_missing_species_prefix": (
                invalid_protein_prefix_count
            ),
            "bed_protein_match_pct": (
                f"{match_percentage:.2f}"
            ),
            "status": status,
            "warnings": ";".join(
                warnings
            ),
        }
    )

    manifest_rows.append(
        {
            "species_code": code,
            "scientific_name": metadata[
                "scientific_name"
            ],
            "ploidy": metadata[
                "ploidy"
            ],
            "analysis_group": metadata[
                "analysis_group"
            ],
            "genome_fasta": str(
                genome_file
            ),
            "bed_file": str(
                bed_file
            ),
            "protein_fasta": str(
                pep_file
            ),
            "id_map": str(
                idmap_file
            ),
            "chromosome_count": len(
                chromosome_lengths
            ),
            "bed_gene_count": len(
                bed_rows
            ),
            "protein_count": len(
                protein_ids
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
# Write QC tables
# ============================================================

qc_fields = [
    "species_code",
    "chromosome_count",
    "bed_gene_count",
    "unique_bed_ids",
    "protein_count",
    "unique_protein_ids",
    "id_map_count",
    "duplicate_bed_ids",
    "duplicate_protein_ids",
    "bed_ids_without_proteins",
    "protein_ids_without_bed",
    "bed_ids_without_id_map",
    "id_map_ids_without_bed",
    "coordinate_issues",
    "gene_order_issues",
    "chromosomes_without_genes",
    "bed_ids_missing_species_prefix",
    "protein_ids_missing_species_prefix",
    "bed_protein_match_pct",
    "status",
    "warnings",
]

qc_output = (
    table_dir / "jcvi_input_qc.tsv"
)

with qc_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=qc_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(qc_rows)

coordinate_fields = [
    "species_code",
    "line_number",
    "gene_id",
    "chromosome",
    "start0",
    "end",
    "chromosome_length",
    "issue",
]

with (
    table_dir
    / "jcvi_coordinate_issues.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=coordinate_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(
        coordinate_issues
    )

with (
    table_dir
    / "jcvi_duplicate_ids.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "file_type",
            "identifier",
            "occurrences",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(
        duplicate_issues
    )

with (
    table_dir
    / "jcvi_id_mismatches.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "identifier",
            "issue",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(
        id_mismatch_rows
    )

with (
    table_dir
    / "jcvi_gene_order_issues.tsv"
).open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "line_number",
            "gene_id",
            "chromosome",
            "start0",
            "issue",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(
        gene_order_issues
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

with output_manifest.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=manifest_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(
        manifest_rows
    )

species_pass = sum(
    row["status"] == "PASS"
    for row in qc_rows
)

species_fail = len(
    qc_rows
) - species_pass

summary_output = (
    table_dir / "jcvi_input_summary.tsv"
)

with summary_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        ["metric", "value"]
    )
    writer.writerow(
        ["species_total", len(qc_rows)]
    )
    writer.writerow(
        ["species_pass", species_pass]
    )
    writer.writerow(
        ["species_fail", species_fail]
    )
    writer.writerow(
        [
            "total_bed_genes",
            sum(
                row["bed_gene_count"]
                for row in qc_rows
            ),
        ]
    )
    writer.writerow(
        [
            "total_proteins",
            sum(
                row["protein_count"]
                for row in qc_rows
            ),
        ]
    )
    writer.writerow(
        [
            "total_duplicate_bed_ids",
            sum(
                row["duplicate_bed_ids"]
                for row in qc_rows
            ),
        ]
    )
    writer.writerow(
        [
            "total_duplicate_protein_ids",
            sum(
                row["duplicate_protein_ids"]
                for row in qc_rows
            ),
        ]
    )
    writer.writerow(
        [
            "total_id_mismatches",
            len(id_mismatch_rows),
        ]
    )
    writer.writerow(
        [
            "total_coordinate_issues",
            len(coordinate_issues),
        ]
    )
    writer.writerow(
        [
            "total_gene_order_issues",
            len(gene_order_issues),
        ]
    )

print("Corrected Step 33 validation completed.")

for row in qc_rows:
    print(
        row["species_code"],
        "BED=",
        row["bed_gene_count"],
        "PEP=",
        row["protein_count"],
        "ID_match=",
        row["bed_protein_match_pct"],
        "duplicates=",
        row["duplicate_bed_ids"],
        "coordinates=",
        row["coordinate_issues"],
        "order=",
        row["gene_order_issues"],
        "status=",
        row["status"],
    )

if species_fail > 0:
    raise SystemExit(
        f"ERROR: {species_fail} species failed corrected Step 33 validation."
    )
PY

# ============================================================
# Display corrected results
# ============================================================

echo
echo "============================================================"
echo "Corrected JCVI input summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/jcvi_input_summary.tsv"

echo
echo "============================================================"
echo "Corrected JCVI input QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/jcvi_input_qc.tsv"

echo
echo "============================================================"
echo "Final JCVI input manifest"
echo "============================================================"

column -t -s $'\t' \
    "${OUTPUT_MANIFEST}"

# ============================================================
# Final strict shell checks
# ============================================================

SPECIES_PASS=$(
    awk -F'\t' '
        $1 == "species_pass" {
            print $2
        }
    ' "${TABLE_DIR}/jcvi_input_summary.tsv"
)

SPECIES_FAIL=$(
    awk -F'\t' '
        $1 == "species_fail" {
            print $2
        }
    ' "${TABLE_DIR}/jcvi_input_summary.tsv"
)

if [[ "${SPECIES_PASS}" -ne 9 ]]; then
    echo "ERROR: Expected 9 PASS species; observed ${SPECIES_PASS}." >&2
    exit 1
fi

if [[ "${SPECIES_FAIL}" -ne 0 ]]; then
    echo "ERROR: Expected 0 failed species; observed ${SPECIES_FAIL}." >&2
    exit 1
fi

UNAPPROVED=$(
    awk -F'\t' '
        NR > 1 && $13 != "YES" {
            count++
        }
        END {
            print count + 0
        }
    ' "${OUTPUT_MANIFEST}"
)

if [[ "${UNAPPROVED}" -ne 0 ]]; then
    echo "ERROR: Final JCVI manifest contains unapproved inputs." >&2
    exit 1
fi

# ============================================================
# Create final checkpoint
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
    "${TABLE_DIR}/jcvi_duplicate_ids.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_id_mismatches.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_gene_order_issues.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${OUTPUT_MANIFEST}" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/JCVI_INPUTS_COMPLETE.txt" <<EOF2
checkpoint=jcvi_bed_protein_inputs
date=$(date --iso-8601=seconds)
conda_environment=jcvi_env
jcvi_version=1.6.5
species_total=9
species_pass=9
species_fail=0
status=PASS
next_step=run_pairwise_and_self_jcvi_synteny
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
echo "Step 33B completed successfully"
echo "============================================================"
echo "Checkpoint:"
echo "${CHECKPOINT_DIR}/JCVI_INPUTS_COMPLETE.txt"
