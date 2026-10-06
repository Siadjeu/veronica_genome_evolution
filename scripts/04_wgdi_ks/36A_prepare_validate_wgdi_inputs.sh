#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36A
#SBATCH --output=11_wgdi/logs/step36A_%j.out
#SBATCH --error=11_wgdi/logs/step36A_%j.err

set -euo pipefail

###############################################################################
# USER CONFIGURATION
###############################################################################

# Keep the real private path here during analysis.
PROJECT_ROOT="$(pwd)"

# The script must be launched from the project root.
EXPECTED_MARKERS=(
    "05_scripts"
    "10_synteny"
)

SPECIES=(
    PMAJ
    VPAN
    VSCU
    VANA
    VARV
    VPER
    VSER
    VTRI
    VVER
)

###############################################################################
# PROJECT INPUTS
###############################################################################

SYNTENY_ROOT="${PROJECT_ROOT}/10_synteny"

BED_DIR="${SYNTENY_ROOT}/jcvi_inputs/bed"
PEP_DIR="${SYNTENY_ROOT}/jcvi_inputs/pep"
GENOME_DIR="${SYNTENY_ROOT}/jcvi_inputs/genomes"

###############################################################################
# WGDI OUTPUTS
###############################################################################

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

ADMIN_DIR="${WGDI_ROOT}/00_admin"
INPUT_DIR="${WGDI_ROOT}/01_inputs"
GFF_DIR="${INPUT_DIR}/gff"
LENS_DIR="${INPUT_DIR}/lens"
PEP_OUT_DIR="${INPUT_DIR}/pep"
CDS_OUT_DIR="${INPUT_DIR}/cds"
MAP_DIR="${INPUT_DIR}/maps"

QC_DIR="${WGDI_ROOT}/02_qc/input_preparation"
CONF_DIR="${WGDI_ROOT}/03_config_templates"
LOG_DIR="${WGDI_ROOT}/logs"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36A"

mkdir -p \
    "${ADMIN_DIR}" \
    "${GFF_DIR}" \
    "${LENS_DIR}" \
    "${PEP_OUT_DIR}" \
    "${CDS_OUT_DIR}" \
    "${MAP_DIR}" \
    "${QC_DIR}" \
    "${CONF_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

###############################################################################
# ENVIRONMENT
###############################################################################

for MARKER in "${EXPECTED_MARKERS[@]}"
do
    if [[ ! -e "${PROJECT_ROOT}/${MARKER}" ]]; then
        echo "ERROR: PROJECT_ROOT does not appear to be correct." >&2
        echo "Missing marker: ${PROJECT_ROOT}/${MARKER}" >&2
        echo "Launch this script from the project root." >&2
        exit 1
    fi
done

cd "${PROJECT_ROOT}"

module purge
module load hpc-env/13.1 2>/dev/null || true

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
else
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate wgdi_env

for PROGRAM in \
    python \
    wgdi
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: Required program not found in wgdi_env: ${PROGRAM}" >&2
        exit 1
    fi
done

echo "============================================================"
echo "Step 36A: prepare and validate WGDI inputs"
echo "============================================================"
echo "Start: $(date --iso-8601=seconds)"
echo "Host: $(hostname)"
echo "Python: $(python --version 2>&1)"
echo "WGDI executable: $(command -v wgdi)"

WGDI_VERSION="$(
    conda list 2>/dev/null |
    awk '$1 == "wgdi" {print $2; exit}'
)"

if [[ -z "${WGDI_VERSION}" ]]; then
    WGDI_VERSION="unknown"
fi

echo "WGDI package version: ${WGDI_VERSION}"

###############################################################################
# CLEAN ONLY STEP 36A OUTPUTS
###############################################################################

rm -f \
    "${ADMIN_DIR}/wgdi_input_manifest.tsv" \
    "${ADMIN_DIR}/wgdi_input_source_manifest.tsv" \
    "${ADMIN_DIR}/wgdi_environment_packages.tsv" \
    "${QC_DIR}/step36A_species_summary.tsv" \
    "${QC_DIR}/step36A_overall_summary.tsv" \
    "${QC_DIR}/step36A_missing_cds_candidates.tsv" \
    "${QC_DIR}/step36A_cds_candidate_inventory.tsv" \
    "${CHECKPOINT_DIR}/STEP36A_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/STEP36A_INCOMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

for CODE in "${SPECIES[@]}"
do
    rm -f \
        "${GFF_DIR}/${CODE}.wgdi.gff" \
        "${LENS_DIR}/${CODE}.wgdi.lens" \
        "${PEP_OUT_DIR}/${CODE}.wgdi.pep.fa" \
        "${CDS_OUT_DIR}/${CODE}.wgdi.cds.fa" \
        "${MAP_DIR}/${CODE}.gene_id_map.tsv" \
        "${QC_DIR}/${CODE}.wgdi_input_qc.tsv" \
        "${QC_DIR}/${CODE}.cds_candidate_scores.tsv"
done

###############################################################################
# RECORD ENVIRONMENT
###############################################################################

conda list |
awk 'BEGIN {
         OFS="\t";
         print "package","version","build","channel"
     }
     !/^#/ && NF >= 3 {
         channel=(NF >= 4 ? $4 : "");
         print $1,$2,$3,channel
     }' \
> "${ADMIN_DIR}/wgdi_environment_packages.tsv"

###############################################################################
# BUILD SOURCE MANIFEST
###############################################################################

SOURCE_MANIFEST="${ADMIN_DIR}/wgdi_input_source_manifest.tsv"

printf '%s\t%s\t%s\t%s\t%s\n' \
    "species_code" \
    "bed_source" \
    "protein_source" \
    "genome_source" \
    "cds_source" \
    > "${SOURCE_MANIFEST}"

resolve_one_file()
{
    local DESCRIPTION="$1"
    shift

    local MATCHES=()
    local CANDIDATE

    for CANDIDATE in "$@"
    do
        if [[ -s "${CANDIDATE}" ]]; then
            MATCHES+=("${CANDIDATE}")
        fi
    done

    if [[ "${#MATCHES[@]}" -ne 1 ]]; then
        echo "ERROR: Expected exactly one ${DESCRIPTION};" >&2
        echo "found ${#MATCHES[@]}." >&2

        if [[ "${#MATCHES[@]}" -gt 0 ]]; then
            printf '  %s\n' "${MATCHES[@]}" >&2
        fi

        return 1
    fi

    printf '%s\n' "${MATCHES[0]}"
}

for CODE in "${SPECIES[@]}"
do
    BED_SOURCE="$(
        resolve_one_file \
            "${CODE} BED file" \
            "${BED_DIR}/${CODE}.bed"
    )"

    PEP_SOURCE="$(
        resolve_one_file \
            "${CODE} protein file" \
            "${PEP_DIR}/${CODE}.pep" \
            "${PEP_DIR}/${CODE}.pep.fa" \
            "${PEP_DIR}/${CODE}.faa"
    )"

    GENOME_SOURCE="$(
        resolve_one_file \
            "${CODE} genome file" \
            "${GENOME_DIR}/${CODE}.fa" \
            "${GENOME_DIR}/${CODE}.fa.gz" \
            "${GENOME_DIR}/${CODE}.fasta" \
            "${GENOME_DIR}/${CODE}.fasta.gz" \
            "${GENOME_DIR}/${CODE}.fna" \
            "${GENOME_DIR}/${CODE}.fna.gz"
    )"

    # CDS is resolved later by exact sequence-ID matching.
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${CODE}" \
        "${BED_SOURCE}" \
        "${PEP_SOURCE}" \
        "${GENOME_SOURCE}" \
        "AUTO_DISCOVER" \
        >> "${SOURCE_MANIFEST}"
done

###############################################################################
# PREPARE WGDI GFF, LENS AND PROTEIN FILES
#
# WGDI GFF columns written here:
#   chromosome
#   gene_id
#   start
#   end
#   strand
#   order
#   original_gene_id
#
# WGDI lens columns:
#   chromosome
#   chromosome_length
#   number_of_genes
###############################################################################

python - \
    "${SOURCE_MANIFEST}" \
    "${GFF_DIR}" \
    "${LENS_DIR}" \
    "${PEP_OUT_DIR}" \
    "${MAP_DIR}" \
    "${QC_DIR}" <<'PY'
from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

source_manifest = Path(sys.argv[1])
gff_dir = Path(sys.argv[2])
lens_dir = Path(sys.argv[3])
pep_out_dir = Path(sys.argv[4])
map_dir = Path(sys.argv[5])
qc_dir = Path(sys.argv[6])


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path):
    header = None
    seq_parts = []

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)

                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise SystemExit(
                        f"ERROR: Sequence before FASTA header in {path}"
                    )

                seq_parts.append(
                    re.sub(r"\s+", "", line)
                )

    if header is not None:
        yield header, "".join(seq_parts)


def fasta_primary_id(header: str) -> str:
    return header.split()[0]


def natural_key(value: str):
    chunks = re.split(r"(\d+)", value)
    return [
        int(chunk) if chunk.isdigit() else chunk
        for chunk in chunks
    ]


with source_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    source_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

required_source_columns = {
    "species_code",
    "bed_source",
    "protein_source",
    "genome_source",
    "cds_source",
}

if not source_rows:
    raise SystemExit("ERROR: Source manifest is empty.")

missing_source_columns = required_source_columns.difference(
    source_rows[0]
)

if missing_source_columns:
    raise SystemExit(
        "ERROR: Source manifest lacks columns: "
        + ",".join(sorted(missing_source_columns))
    )

for source_row in source_rows:
    species = source_row["species_code"]
    bed_path = Path(source_row["bed_source"])
    pep_path = Path(source_row["protein_source"])
    genome_path = Path(source_row["genome_source"])

    for required_path in [
        bed_path,
        pep_path,
        genome_path,
    ]:
        if (
            not required_path.is_file()
            or required_path.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing input for {species}: {required_path}"
            )

    chromosome_lengths = {}

    for header, sequence in read_fasta(genome_path):
        sequence_id = fasta_primary_id(header)

        if sequence_id in chromosome_lengths:
            raise SystemExit(
                f"ERROR: Duplicate genome FASTA ID for "
                f"{species}: {sequence_id}"
            )

        chromosome_lengths[sequence_id] = len(sequence)

    if not chromosome_lengths:
        raise SystemExit(
            f"ERROR: No genome sequences read for {species}"
        )

    protein_records = {}

    for header, sequence in read_fasta(pep_path):
        protein_id = fasta_primary_id(header)

        if protein_id in protein_records:
            raise SystemExit(
                f"ERROR: Duplicate protein ID for "
                f"{species}: {protein_id}"
            )

        sequence = sequence.upper().rstrip("*")

        if not sequence:
            raise SystemExit(
                f"ERROR: Empty protein sequence for "
                f"{species}: {protein_id}"
            )

        invalid = set(sequence).difference(
            set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-.")
        )

        if invalid:
            raise SystemExit(
                f"ERROR: Invalid protein characters for "
                f"{species}:{protein_id}: {sorted(invalid)}"
            )

        protein_records[protein_id] = sequence

    bed_records = []
    seen_gene_ids = set()
    chromosome_gene_counts = Counter()

    with bed_path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line_number, raw_line in enumerate(
            handle,
            start=1,
        ):
            line = raw_line.rstrip("\n\r")

            if (
                not line
                or line.startswith("#")
            ):
                continue

            fields = line.split("\t")

            if len(fields) < 4:
                fields = line.split()

            if len(fields) < 4:
                raise SystemExit(
                    f"ERROR: BED line has fewer than four fields: "
                    f"{bed_path}:{line_number}"
                )

            chromosome = fields[0]
            start0 = int(fields[1])
            end = int(fields[2])
            gene_id = fields[3]

            strand = "+"

            if len(fields) >= 6 and fields[5] in {"+", "-"}:
                strand = fields[5]
            elif len(fields) >= 5 and fields[4] in {"+", "-"}:
                strand = fields[4]

            if chromosome not in chromosome_lengths:
                raise SystemExit(
                    f"ERROR: BED chromosome absent from genome for "
                    f"{species}: {chromosome}"
                )

            if gene_id in seen_gene_ids:
                raise SystemExit(
                    f"ERROR: Duplicate BED gene ID for "
                    f"{species}: {gene_id}"
                )

            if gene_id not in protein_records:
                raise SystemExit(
                    f"ERROR: BED gene lacks exact protein ID for "
                    f"{species}: {gene_id}"
                )

            if start0 < 0:
                raise SystemExit(
                    f"ERROR: Negative BED start for "
                    f"{species}:{gene_id}"
                )

            if end <= start0:
                raise SystemExit(
                    f"ERROR: Invalid BED interval for "
                    f"{species}:{gene_id}: {start0}-{end}"
                )

            if end > chromosome_lengths[chromosome]:
                raise SystemExit(
                    f"ERROR: BED interval exceeds chromosome length for "
                    f"{species}:{gene_id}"
                )

            # BED is 0-based, half-open.
            # WGDI coordinates are written as 1-based inclusive.
            start1 = start0 + 1

            bed_records.append(
                {
                    "chromosome": chromosome,
                    "start": start1,
                    "end": end,
                    "gene_id": gene_id,
                    "strand": strand,
                }
            )

            seen_gene_ids.add(gene_id)
            chromosome_gene_counts[chromosome] += 1

    if not bed_records:
        raise SystemExit(
            f"ERROR: No BED genes read for {species}"
        )

    extra_proteins = set(protein_records).difference(
        seen_gene_ids
    )

    if extra_proteins:
        raise SystemExit(
            f"ERROR: {species} has {len(extra_proteins)} proteins "
            "not represented in the validated BED. "
            f"Examples: {sorted(extra_proteins)[:10]}"
        )

    bed_records.sort(
        key=lambda row: (
            natural_key(row["chromosome"]),
            row["start"],
            row["end"],
            row["gene_id"],
        )
    )

    chromosome_order = defaultdict(int)

    gff_path = gff_dir / f"{species}.wgdi.gff"
    lens_path = lens_dir / f"{species}.wgdi.lens"
    pep_output = pep_out_dir / f"{species}.wgdi.pep.fa"
    map_output = map_dir / f"{species}.gene_id_map.tsv"
    qc_output = qc_dir / f"{species}.wgdi_input_qc.tsv"

    with gff_path.open(
        "w",
        encoding="utf-8",
    ) as gff_handle, map_output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as map_handle:
        map_fields = [
            "species_code",
            "chromosome",
            "wgdi_gene_id",
            "original_gene_id",
            "start",
            "end",
            "strand",
            "order",
        ]

        map_writer = csv.DictWriter(
            map_handle,
            fieldnames=map_fields,
            delimiter="\t",
            lineterminator="\n",
        )

        map_writer.writeheader()

        for row in bed_records:
            chromosome = row["chromosome"]
            chromosome_order[chromosome] += 1
            order = chromosome_order[chromosome]

            # We retain exact validated IDs rather than renaming them.
            wgdi_gene_id = row["gene_id"]

            gff_handle.write(
                "\t".join(
                    [
                        chromosome,
                        wgdi_gene_id,
                        str(row["start"]),
                        str(row["end"]),
                        row["strand"],
                        str(order),
                        row["gene_id"],
                    ]
                )
                + "\n"
            )

            map_writer.writerow(
                {
                    "species_code": species,
                    "chromosome": chromosome,
                    "wgdi_gene_id": wgdi_gene_id,
                    "original_gene_id": row["gene_id"],
                    "start": row["start"],
                    "end": row["end"],
                    "strand": row["strand"],
                    "order": order,
                }
            )

    used_chromosomes = sorted(
        chromosome_gene_counts,
        key=natural_key,
    )

    with lens_path.open(
        "w",
        encoding="utf-8",
    ) as lens_handle:
        for chromosome in used_chromosomes:
            lens_handle.write(
                "\t".join(
                    [
                        chromosome,
                        str(chromosome_lengths[chromosome]),
                        str(chromosome_gene_counts[chromosome]),
                    ]
                )
                + "\n"
            )

    with pep_output.open(
        "w",
        encoding="utf-8",
    ) as pep_handle:
        for row in bed_records:
            gene_id = row["gene_id"]
            sequence = protein_records[gene_id]

            pep_handle.write(f">{gene_id}\n")

            for start in range(0, len(sequence), 60):
                pep_handle.write(
                    sequence[start:start + 60] + "\n"
                )

    total_genome_bp = sum(
        chromosome_lengths[chromosome]
        for chromosome in used_chromosomes
    )

    total_gene_span = sum(
        row["end"] - row["start"] + 1
        for row in bed_records
    )

    qc_rows = [
        ("species_code", species),
        ("input_bed", str(bed_path)),
        ("input_protein", str(pep_path)),
        ("input_genome", str(genome_path)),
        ("wgdi_gff", str(gff_path)),
        ("wgdi_lens", str(lens_path)),
        ("wgdi_protein", str(pep_output)),
        ("genome_sequence_count", str(len(chromosome_lengths))),
        ("chromosomes_with_genes", str(len(used_chromosomes))),
        ("gene_count", str(len(bed_records))),
        ("protein_count", str(len(protein_records))),
        ("bed_protein_exact_match", "PASS"),
        ("duplicate_gene_ids", "0"),
        ("duplicate_protein_ids", "0"),
        ("coordinate_failures", "0"),
        ("total_used_chromosome_bp", str(total_genome_bp)),
        ("total_gene_span_bp", str(total_gene_span)),
        ("status_without_cds", "PASS"),
    ]

    with qc_output.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as qc_handle:
        writer = csv.writer(
            qc_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(["metric", "value"])
        writer.writerows(qc_rows)

    print(
        f"{species}: genes={len(bed_records):,}; "
        f"chromosomes={len(used_chromosomes):,}; "
        f"BED/PEP=PASS"
    )
PY

###############################################################################
# DISCOVER CDS CANDIDATES
#
# The script does not automatically trust filenames. It finds plausible CDS
# FASTAs and scores them by exact overlap with each validated WGDI protein-ID
# set. A CDS file is accepted only when it provides an exact one-to-one ID
# match for all genes of one species.
###############################################################################

CDS_INVENTORY="${QC_DIR}/step36A_cds_candidate_inventory.tsv"
MISSING_CDS="${QC_DIR}/step36A_missing_cds_candidates.tsv"

printf '%s\t%s\t%s\n' \
    "candidate_path" \
    "size_bytes" \
    "basename" \
    > "${CDS_INVENTORY}"

# Search likely annotation/input directories while excluding WGDI outputs.
# This searches filenames only and does not parse every file yet.
while IFS= read -r -d '' FILE
do
    SIZE="$(stat -c '%s' "${FILE}")"

    printf '%s\t%s\t%s\n' \
        "${FILE}" \
        "${SIZE}" \
        "$(basename "${FILE}")" \
        >> "${CDS_INVENTORY}"
done < <(
    find "${PROJECT_ROOT}" \
        -type f \
        \( \
            -iname '*cds*.fa' \
            -o -iname '*cds*.fa.gz' \
            -o -iname '*cds*.fasta' \
            -o -iname '*cds*.fasta.gz' \
            -o -iname '*cds*.fna' \
            -o -iname '*cds*.fna.gz' \
            -o -iname '*codingseq*' \
            -o -iname '*codingseq*.gz' \
        \) \
        ! -path "${WGDI_ROOT}/*" \
        -print0
)

python - \
    "${SOURCE_MANIFEST}" \
    "${CDS_INVENTORY}" \
    "${PEP_OUT_DIR}" \
    "${CDS_OUT_DIR}" \
    "${QC_DIR}" \
    "${MISSING_CDS}" <<'PY'
from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
from pathlib import Path

source_manifest = Path(sys.argv[1])
inventory_file = Path(sys.argv[2])
pep_dir = Path(sys.argv[3])
cds_out_dir = Path(sys.argv[4])
qc_dir = Path(sys.argv[5])
missing_file = Path(sys.argv[6])


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_fasta(path: Path):
    header = None
    seq_parts = []

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)

                header = line[1:].strip()
                seq_parts = []
            else:
                if header is None:
                    raise ValueError(
                        f"Sequence before FASTA header in {path}"
                    )

                seq_parts.append(
                    re.sub(r"\s+", "", line)
                )

    if header is not None:
        yield header, "".join(seq_parts)


def primary_id(header: str) -> str:
    return header.split()[0]


def load_ids(path: Path):
    ids = []
    seen = set()

    for header, _sequence in read_fasta(path):
        sequence_id = primary_id(header)

        if sequence_id in seen:
            raise ValueError(
                f"Duplicate FASTA ID in {path}: {sequence_id}"
            )

        seen.add(sequence_id)
        ids.append(sequence_id)

    return ids, seen


with source_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    source_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

with inventory_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    inventory_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

candidate_paths = [
    Path(row["candidate_path"])
    for row in inventory_rows
    if Path(row["candidate_path"]).is_file()
]

missing_rows = []
selected_by_species = {}

for source_row in source_rows:
    species = source_row["species_code"]
    pep_path = pep_dir / f"{species}.wgdi.pep.fa"

    _pep_order, pep_ids = load_ids(pep_path)

    score_rows = []

    for candidate in candidate_paths:
        try:
            _cds_order, cds_ids = load_ids(candidate)
        except Exception as exc:
            score_rows.append(
                {
                    "species_code": species,
                    "candidate_path": str(candidate),
                    "protein_ids": len(pep_ids),
                    "cds_ids": 0,
                    "intersection": 0,
                    "protein_missing_from_cds": len(pep_ids),
                    "extra_cds_ids": 0,
                    "protein_recall": "0.00000000",
                    "exact_id_set_match": "NO",
                    "parse_status": f"FAIL:{type(exc).__name__}",
                }
            )
            continue

        intersection = pep_ids.intersection(cds_ids)
        missing = pep_ids.difference(cds_ids)
        extra = cds_ids.difference(pep_ids)

        exact = (
            not missing
            and not extra
            and len(cds_ids) == len(pep_ids)
        )

        score_rows.append(
            {
                "species_code": species,
                "candidate_path": str(candidate),
                "protein_ids": len(pep_ids),
                "cds_ids": len(cds_ids),
                "intersection": len(intersection),
                "protein_missing_from_cds": len(missing),
                "extra_cds_ids": len(extra),
                "protein_recall": (
                    f"{len(intersection) / len(pep_ids):.8f}"
                ),
                "exact_id_set_match": "YES" if exact else "NO",
                "parse_status": "PASS",
            }
        )

    score_rows.sort(
        key=lambda row: (
            row["exact_id_set_match"] == "YES",
            float(row["protein_recall"]),
            row["intersection"],
            -row["extra_cds_ids"],
        ),
        reverse=True,
    )

    score_file = qc_dir / f"{species}.cds_candidate_scores.tsv"

    with score_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        fieldnames = [
            "species_code",
            "candidate_path",
            "protein_ids",
            "cds_ids",
            "intersection",
            "protein_missing_from_cds",
            "extra_cds_ids",
            "protein_recall",
            "exact_id_set_match",
            "parse_status",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(score_rows)

    exact_matches = [
        Path(row["candidate_path"])
        for row in score_rows
        if row["exact_id_set_match"] == "YES"
    ]

    if len(exact_matches) != 1:
        best = score_rows[0] if score_rows else None

        missing_rows.append(
            {
                "species_code": species,
                "exact_candidate_count": len(exact_matches),
                "best_candidate": (
                    best["candidate_path"]
                    if best
                    else "NONE"
                ),
                "best_protein_recall": (
                    best["protein_recall"]
                    if best
                    else "0.00000000"
                ),
                "best_missing_ids": (
                    best["protein_missing_from_cds"]
                    if best
                    else len(pep_ids)
                ),
                "best_extra_ids": (
                    best["extra_cds_ids"]
                    if best
                    else 0
                ),
                "status": "REVIEW_REQUIRED",
            }
        )
        continue

    selected_cds = exact_matches[0]
    selected_by_species[species] = selected_cds

    protein_order, _pep_set = load_ids(pep_path)

    cds_sequences = {}

    for header, sequence in read_fasta(selected_cds):
        sequence_id = primary_id(header)
        sequence = sequence.upper().replace("U", "T")

        invalid = set(sequence).difference(
            set("ACGTRYSWKMBDHVN")
        )

        if invalid:
            raise SystemExit(
                f"ERROR: Invalid CDS characters in "
                f"{species}:{sequence_id}: {sorted(invalid)}"
            )

        if not sequence:
            raise SystemExit(
                f"ERROR: Empty CDS for {species}:{sequence_id}"
            )

        cds_sequences[sequence_id] = sequence

    output_path = cds_out_dir / f"{species}.wgdi.cds.fa"

    with output_path.open(
        "w",
        encoding="utf-8",
    ) as handle:
        for gene_id in protein_order:
            sequence = cds_sequences[gene_id]

            handle.write(f">{gene_id}\n")

            for start in range(0, len(sequence), 60):
                handle.write(
                    sequence[start:start + 60] + "\n"
                )

if missing_rows:
    with missing_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        fieldnames = [
            "species_code",
            "exact_candidate_count",
            "best_candidate",
            "best_protein_recall",
            "best_missing_ids",
            "best_extra_ids",
            "status",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(missing_rows)

print(
    f"Exact CDS datasets resolved: "
    f"{len(selected_by_species)}/{len(source_rows)}"
)
PY

###############################################################################
# FINAL VALIDATION AND SUMMARY
###############################################################################

python - \
    "${SOURCE_MANIFEST}" \
    "${GFF_DIR}" \
    "${LENS_DIR}" \
    "${PEP_OUT_DIR}" \
    "${CDS_OUT_DIR}" \
    "${QC_DIR}" \
    "${ADMIN_DIR}/wgdi_input_manifest.tsv" \
    "${QC_DIR}/step36A_species_summary.tsv" <<'PY'
from __future__ import annotations

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

source_manifest = Path(sys.argv[1])
gff_dir = Path(sys.argv[2])
lens_dir = Path(sys.argv[3])
pep_dir = Path(sys.argv[4])
cds_dir = Path(sys.argv[5])
qc_dir = Path(sys.argv[6])
manifest_output = Path(sys.argv[7])
summary_output = Path(sys.argv[8])


def read_fasta_ids(path: Path):
    ids = []
    seen = set()

    with path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            sequence_id = line[1:].split()[0]

            if sequence_id in seen:
                raise SystemExit(
                    f"ERROR: Duplicate FASTA ID in {path}: "
                    f"{sequence_id}"
                )

            seen.add(sequence_id)
            ids.append(sequence_id)

    return ids, seen


with source_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    source_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

manifest_rows = []
summary_rows = []

for source_row in source_rows:
    species = source_row["species_code"]

    gff_path = gff_dir / f"{species}.wgdi.gff"
    lens_path = lens_dir / f"{species}.wgdi.lens"
    pep_path = pep_dir / f"{species}.wgdi.pep.fa"
    cds_path = cds_dir / f"{species}.wgdi.cds.fa"

    required_non_cds = [
        gff_path,
        lens_path,
        pep_path,
    ]

    for path in required_non_cds:
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing prepared file for {species}: {path}"
            )

    gff_ids = []
    gff_seen = set()
    chromosome_orders = defaultdict(list)
    chromosome_gene_counts = Counter()

    with gff_path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            fields = line.rstrip("\n").split("\t")

            if len(fields) != 7:
                raise SystemExit(
                    f"ERROR: WGDI GFF must have seven columns: "
                    f"{gff_path}:{line_number}"
                )

            chromosome, gene_id, start, end, strand, order, original = fields

            start = int(start)
            end = int(end)
            order = int(order)

            if gene_id in gff_seen:
                raise SystemExit(
                    f"ERROR: Duplicate WGDI GFF ID: "
                    f"{species}:{gene_id}"
                )

            if start < 1 or end < start:
                raise SystemExit(
                    f"ERROR: Invalid WGDI coordinates: "
                    f"{species}:{gene_id}"
                )

            if strand not in {"+", "-"}:
                raise SystemExit(
                    f"ERROR: Invalid strand: {species}:{gene_id}"
                )

            if original != gene_id:
                raise SystemExit(
                    f"ERROR: Original ID differs unexpectedly: "
                    f"{species}:{gene_id}:{original}"
                )

            gff_seen.add(gene_id)
            gff_ids.append(gene_id)
            chromosome_orders[chromosome].append(order)
            chromosome_gene_counts[chromosome] += 1

    for chromosome, orders in chromosome_orders.items():
        expected = list(
            range(1, len(orders) + 1)
        )

        if orders != expected:
            raise SystemExit(
                f"ERROR: Nonconsecutive gene order for "
                f"{species}:{chromosome}"
            )

    lens_counts = {}
    lens_lengths = {}

    with lens_path.open(
        "r",
        encoding="utf-8",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            fields = line.rstrip("\n").split("\t")

            if len(fields) != 3:
                raise SystemExit(
                    f"ERROR: WGDI lens must have three columns: "
                    f"{lens_path}:{line_number}"
                )

            chromosome = fields[0]
            length = int(fields[1])
            gene_count = int(fields[2])

            if chromosome in lens_counts:
                raise SystemExit(
                    f"ERROR: Duplicate lens chromosome: "
                    f"{species}:{chromosome}"
                )

            if length <= 0 or gene_count <= 0:
                raise SystemExit(
                    f"ERROR: Invalid lens values: "
                    f"{species}:{chromosome}"
                )

            lens_lengths[chromosome] = length
            lens_counts[chromosome] = gene_count

    if set(lens_counts) != set(chromosome_gene_counts):
        raise SystemExit(
            f"ERROR: GFF/lens chromosome mismatch for {species}"
        )

    for chromosome in lens_counts:
        if lens_counts[chromosome] != chromosome_gene_counts[chromosome]:
            raise SystemExit(
                f"ERROR: GFF/lens gene-count mismatch for "
                f"{species}:{chromosome}"
            )

    pep_order, pep_ids = read_fasta_ids(pep_path)

    if pep_ids != gff_seen:
        raise SystemExit(
            f"ERROR: WGDI GFF/protein ID mismatch for {species}"
        )

    if pep_order != gff_ids:
        raise SystemExit(
            f"ERROR: WGDI protein order differs from GFF order "
            f"for {species}"
        )

    cds_status = "MISSING"
    cds_count = 0
    exact_cds_match = "NO"

    if cds_path.is_file() and cds_path.stat().st_size > 0:
        cds_order, cds_ids = read_fasta_ids(cds_path)
        cds_count = len(cds_ids)

        if cds_ids != gff_seen:
            raise SystemExit(
                f"ERROR: WGDI GFF/CDS ID mismatch for {species}"
            )

        if cds_order != gff_ids:
            raise SystemExit(
                f"ERROR: WGDI CDS order differs from GFF order "
                f"for {species}"
            )

        cds_status = "PASS"
        exact_cds_match = "YES"

    status = (
        "PASS"
        if cds_status == "PASS"
        else "CDS_REVIEW_REQUIRED"
    )

    manifest_rows.append(
        {
            "species_code": species,
            "gff": str(gff_path),
            "lens": str(lens_path),
            "pep": str(pep_path),
            "cds": (
                str(cds_path)
                if cds_status == "PASS"
                else ""
            ),
            "gene_id_map": str(
                gff_dir.parent / "maps"
                / f"{species}.gene_id_map.tsv"
            ),
            "gene_count": len(gff_ids),
            "chromosome_count": len(lens_counts),
            "status": status,
        }
    )

    summary_rows.append(
        {
            "species_code": species,
            "chromosome_count": len(lens_counts),
            "gene_count": len(gff_ids),
            "protein_count": len(pep_ids),
            "cds_count": cds_count,
            "gff_lens_match": "PASS",
            "gff_pep_exact_id_match": "PASS",
            "gff_cds_exact_id_match": exact_cds_match,
            "gene_order_validation": "PASS",
            "coordinate_validation": "PASS",
            "status": status,
        }
    )

manifest_fields = [
    "species_code",
    "gff",
    "lens",
    "pep",
    "cds",
    "gene_id_map",
    "gene_count",
    "chromosome_count",
    "status",
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
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(manifest_rows)

summary_fields = [
    "species_code",
    "chromosome_count",
    "gene_count",
    "protein_count",
    "cds_count",
    "gff_lens_match",
    "gff_pep_exact_id_match",
    "gff_cds_exact_id_match",
    "gene_order_validation",
    "coordinate_validation",
    "status",
]

with summary_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(summary_rows)

print(
    "Final WGDI species validation completed."
)
PY

###############################################################################
# CREATE WGDI CONFIGURATION TEMPLATES
###############################################################################

{
    echo "# Generated with WGDI ${WGDI_VERSION}"
    echo "# Review each parameter before execution."
    echo
    wgdi -d ? 2>&1 || true
} > "${CONF_DIR}/dotplot.template.conf"

{
    echo "# Generated with WGDI ${WGDI_VERSION}"
    echo "# Review each parameter before execution."
    echo
    wgdi -icl ? 2>&1 || true
} > "${CONF_DIR}/collinearity.template.conf"

{
    echo "# Generated with WGDI ${WGDI_VERSION}"
    echo "# Review each parameter before execution."
    echo
    wgdi -ks ? 2>&1 || true
} > "${CONF_DIR}/ks.template.conf"

###############################################################################
# CHECK COMPLETION
###############################################################################

PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${QC_DIR}/step36A_species_summary.tsv"
)"

REVIEW_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "CDS_REVIEW_REQUIRED" {
            count++
        }
        END {
            print count + 0
        }
    ' "${QC_DIR}/step36A_species_summary.tsv"
)"

TOTAL_COUNT="$(
    awk -F $'\t' '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
    ' "${QC_DIR}/step36A_species_summary.tsv"
)"

printf '%s\t%s\n' \
    "metric" \
    "value" \
    > "${QC_DIR}/step36A_overall_summary.tsv"

printf '%s\t%s\n' \
    "species_expected" \
    "${#SPECIES[@]}" \
    >> "${QC_DIR}/step36A_overall_summary.tsv"

printf '%s\t%s\n' \
    "species_evaluated" \
    "${TOTAL_COUNT}" \
    >> "${QC_DIR}/step36A_overall_summary.tsv"

printf '%s\t%s\n' \
    "species_complete" \
    "${PASS_COUNT}" \
    >> "${QC_DIR}/step36A_overall_summary.tsv"

printf '%s\t%s\n' \
    "species_requiring_cds_review" \
    "${REVIEW_COUNT}" \
    >> "${QC_DIR}/step36A_overall_summary.tsv"

printf '%s\t%s\n' \
    "wgdi_version" \
    "${WGDI_VERSION}" \
    >> "${QC_DIR}/step36A_overall_summary.tsv"

if [[ "${PASS_COUNT}" -eq "${#SPECIES[@]}" ]]; then
    OVERALL_STATUS="PASS"

    cat > "${CHECKPOINT_DIR}/STEP36A_COMPLETE.txt" <<EOF2
checkpoint=step36A_prepare_validate_wgdi_inputs
date=$(date --iso-8601=seconds)
species_expected=${#SPECIES[@]}
species_complete=${PASS_COUNT}
wgdi_version=${WGDI_VERSION}
gff_columns=chromosome,gene_id,start,end,strand,order,original_gene_id
lens_columns=chromosome,length_bp,number_of_genes
identifier_policy=exact_validated_gene_ids
bed_protein_match=exact
protein_cds_match=exact
status=PASS
next_step=step36B_generate_self_and_pairwise_protein_homology
EOF2
else
    OVERALL_STATUS="CDS_REVIEW_REQUIRED"

    cat > "${CHECKPOINT_DIR}/STEP36A_INCOMPLETE.txt" <<EOF2
checkpoint=step36A_prepare_validate_wgdi_inputs
date=$(date --iso-8601=seconds)
species_expected=${#SPECIES[@]}
species_complete=${PASS_COUNT}
species_requiring_cds_review=${REVIEW_COUNT}
wgdi_version=${WGDI_VERSION}
gff_lens_pep_status=PASS
cds_status=REVIEW_REQUIRED
overall_status=INCOMPLETE
next_action=review_species_cds_candidate_scores_and_resolve_exact_CDS_files
EOF2
fi

printf '%s\t%s\n' \
    "overall_status" \
    "${OVERALL_STATUS}" \
    >> "${QC_DIR}/step36A_overall_summary.tsv"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# REPORT
###############################################################################

echo
echo "============================================================"
echo "Step 36A species summary"
echo "============================================================"

column -t -s $'\t' \
    "${QC_DIR}/step36A_species_summary.tsv"

echo
echo "============================================================"
echo "Step 36A overall summary"
echo "============================================================"

column -t -s $'\t' \
    "${QC_DIR}/step36A_overall_summary.tsv"

if [[ -s "${MISSING_CDS}" ]]; then
    echo
    echo "============================================================"
    echo "CDS datasets requiring review"
    echo "============================================================"

    column -t -s $'\t' \
        "${MISSING_CDS}"

    echo
    echo "The BED, lens and protein inputs were prepared successfully."
    echo "Step 36A remains incomplete until all CDS IDs match exactly."
    echo
    echo "Inspect individual candidate tables with:"
    echo "  column -t -s \$'\\t' 11_wgdi/02_qc/input_preparation/<SPECIES>.cds_candidate_scores.tsv | head -n 20"
else
    echo
    echo "Step 36A completed successfully."
fi

echo
echo "End: $(date --iso-8601=seconds)"
