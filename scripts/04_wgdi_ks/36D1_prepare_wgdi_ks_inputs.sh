#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36D1
#SBATCH --output=11_wgdi/logs/step36D1_%j.out
#SBATCH --error=11_wgdi/logs/step36D1_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

COLLINEARITY_MANIFEST="${WGDI_ROOT}/00_admin/wgdi_collinearity_manifest.tsv"

PEP_DIR="${WGDI_ROOT}/01_inputs/pep"
CDS_DIR="${WGDI_ROOT}/01_inputs/cds"

KS_ROOT="${WGDI_ROOT}/06_ks"
PAIR_DIR="${KS_ROOT}/01_pairs"
SEQUENCE_DIR="${KS_ROOT}/02_sequences"
CONFIG_DIR="${KS_ROOT}/03_configs"
RESULT_DIR="${KS_ROOT}/04_results"

QC_DIR="${WGDI_ROOT}/02_qc/ks_preparation"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36D1"
LOG_DIR="${WGDI_ROOT}/logs"

KS_MANIFEST="${ADMIN_DIR}/step36D_ks_manifest.tsv"
SUMMARY="${QC_DIR}/step36D1_ks_input_summary.tsv"

mkdir -p \
    "${PAIR_DIR}" \
    "${SEQUENCE_DIR}" \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

###############################################################################
# ENVIRONMENT
###############################################################################

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

for PROGRAM in python wgdi mafft yn00
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: Required program unavailable: ${PROGRAM}" >&2
        exit 1
    fi
done

if [[ ! -s "${COLLINEARITY_MANIFEST}" ]]; then
    echo "ERROR: Collinearity manifest is missing:" >&2
    echo "${COLLINEARITY_MANIFEST}" >&2
    exit 1
fi

###############################################################################
# TARGETED CLEANUP
###############################################################################

rm -f \
    "${KS_MANIFEST}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36D1_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# PREPARE STRICT COLLINEARITY PAIRS AND CONCATENATED SEQUENCES
###############################################################################

python - \
    "${PROJECT_ROOT}" \
    "${COLLINEARITY_MANIFEST}" \
    "${PEP_DIR}" \
    "${CDS_DIR}" \
    "${PAIR_DIR}" \
    "${SEQUENCE_DIR}" \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${KS_MANIFEST}" \
    "${SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
collinearity_manifest = Path(sys.argv[2])
pep_dir = Path(sys.argv[3])
cds_dir = Path(sys.argv[4])
pair_dir = Path(sys.argv[5])
sequence_dir = Path(sys.argv[6])
config_dir = Path(sys.argv[7])
result_dir = Path(sys.argv[8])
ks_manifest = Path(sys.argv[9])
summary_path = Path(sys.argv[10])


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


def read_fasta(path: Path):
    order = []
    records = {}

    header = None
    sequence_parts = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")

            if line.startswith(">"):
                if header is not None:
                    identifier = header.split()[0]

                    if identifier in records:
                        raise SystemExit(
                            f"ERROR: Duplicate FASTA ID in {path}: "
                            f"{identifier}"
                        )

                    order.append(identifier)
                    records[identifier] = "".join(sequence_parts)

                header = line[1:].strip()
                sequence_parts = []

            else:
                if header is None:
                    raise SystemExit(
                        f"ERROR: Sequence before FASTA header in {path}"
                    )

                sequence_parts.append(
                    re.sub(r"\s+", "", line)
                )

    if header is not None:
        identifier = header.split()[0]

        if identifier in records:
            raise SystemExit(
                f"ERROR: Duplicate FASTA ID in {path}: {identifier}"
            )

        order.append(identifier)
        records[identifier] = "".join(sequence_parts)

    if not records:
        raise SystemExit(f"ERROR: Empty FASTA: {path}")

    return order, records


def write_fasta(path: Path, order, records):
    temporary = Path(str(path) + ".tmp")

    with temporary.open("w", encoding="utf-8") as handle:
        for identifier in order:
            sequence = records[identifier]

            handle.write(f">{identifier}\n")

            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")

    temporary.replace(path)


def merge_fastas(paths, output):
    merged_order = []
    merged_records = {}

    for path in paths:
        order, records = read_fasta(path)

        for identifier in order:
            if identifier in merged_records:
                raise SystemExit(
                    f"ERROR: Duplicate ID while merging FASTAs: "
                    f"{identifier}"
                )

            merged_order.append(identifier)
            merged_records[identifier] = records[identifier]

    write_fasta(
        output,
        merged_order,
        merged_records,
    )

    return merged_order, merged_records


with collinearity_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

strict_rows = [
    row
    for row in rows
    if row["parameter_set"] == "strict"
]

if len(strict_rows) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 strict collinearity runs; "
        f"found {len(strict_rows)}."
    )

strict_rows.sort(
    key=lambda row: int(row["task_id"])
)

summary_rows = []
manifest_rows = []

for ks_task_id, row in enumerate(strict_rows):
    species1 = row["species1"]
    species2 = row["species2"]
    comparison = row["comparison"]
    comparison_type = row["comparison_type"]

    collinearity_file = (
        project_root
        / row["collinearity"]
    )

    if (
        not collinearity_file.is_file()
        or collinearity_file.stat().st_size == 0
    ):
        raise SystemExit(
            f"ERROR: Missing strict collinearity file: "
            f"{collinearity_file}"
        )

    pep_paths = [
        pep_dir / f"{species1}.wgdi.pep.fa"
    ]

    cds_paths = [
        cds_dir / f"{species1}.wgdi.cds.fa"
    ]

    if species2 != species1:
        pep_paths.append(
            pep_dir / f"{species2}.wgdi.pep.fa"
        )

        cds_paths.append(
            cds_dir / f"{species2}.wgdi.cds.fa"
        )

    for path in pep_paths + cds_paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: Missing sequence input: {path}"
            )

    pair_file = (
        pair_dir
        / f"{comparison}.strict.pairs.tsv"
    )

    combined_pep = (
        sequence_dir
        / f"{comparison}.pep.fa"
    )

    combined_cds = (
        sequence_dir
        / f"{comparison}.cds.fa"
    )

    config_file = (
        config_dir
        / f"{comparison}.ks.conf"
    )

    ks_result = (
        result_dir
        / f"{comparison}.ks.tsv"
    )

    raw_assignments = 0
    self_pairs_removed = 0
    repeated_assignments_removed = 0
    reverse_duplicates_removed = 0

    pairs = []
    seen = set()

    with collinearity_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 5:
                raise SystemExit(
                    f"ERROR: Malformed WGDI row at "
                    f"{collinearity_file}:{line_number}"
                )

            gene1 = fields[0]
            gene2 = fields[2]

            raw_assignments += 1

            if gene1 == gene2:
                self_pairs_removed += 1
                continue

            if comparison_type == "self":
                canonical = tuple(
                    sorted((gene1, gene2))
                )

                if canonical in seen:
                    reverse_duplicates_removed += 1
                    continue

                pair = canonical

            else:
                pair = (gene1, gene2)

                if pair in seen:
                    repeated_assignments_removed += 1
                    continue

            seen.add(pair)
            pairs.append(pair)

    pairs.sort()

    if not pairs:
        raise SystemExit(
            f"ERROR: No nonself pairs retained for {comparison}."
        )

    pep_order, pep_records = merge_fastas(
        pep_paths,
        combined_pep,
    )

    cds_order, cds_records = merge_fastas(
        cds_paths,
        combined_cds,
    )

    if pep_order != cds_order:
        raise SystemExit(
            f"ERROR: Protein/CDS order differs for {comparison}."
        )

    missing_pep = set()
    missing_cds = set()

    genes_in_pairs = set()

    for gene1, gene2 in pairs:
        genes_in_pairs.add(gene1)
        genes_in_pairs.add(gene2)

        for gene in (gene1, gene2):
            if gene not in pep_records:
                missing_pep.add(gene)

            if gene not in cds_records:
                missing_cds.add(gene)

    if missing_pep:
        raise SystemExit(
            f"ERROR: {len(missing_pep)} genes lack proteins "
            f"for {comparison}."
        )

    if missing_cds:
        raise SystemExit(
            f"ERROR: {len(missing_cds)} genes lack CDS "
            f"for {comparison}."
        )

    with pair_file.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerows(pairs)

    config_file.write_text(
        "\n".join(
            [
                "[ks]",
                f"cds_file = {combined_cds}",
                f"pep_file = {combined_pep}",
                "align_software = mafft",
                f"pairs_file = {pair_file}",
                f"ks_file = {ks_result}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    summary_rows.append(
        {
            "task_id": str(ks_task_id),
            "comparison": comparison,
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "raw_collinearity_pair_assignments": str(
                raw_assignments
            ),
            "self_pairs_removed": str(
                self_pairs_removed
            ),
            "reverse_or_repeated_pairs_removed": str(
                reverse_duplicates_removed
                + repeated_assignments_removed
            ),
            "final_unique_nonself_pairs": str(
                len(pairs)
            ),
            "unique_genes_in_pairs": str(
                len(genes_in_pairs)
            ),
            "combined_protein_count": str(
                len(pep_records)
            ),
            "combined_cds_count": str(
                len(cds_records)
            ),
            "missing_protein_ids": "0",
            "missing_cds_ids": "0",
            "status": "PASS",
        }
    )

    manifest_rows.append(
        {
            "task_id": str(ks_task_id),
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "comparison": comparison,
            "parameter_set": "strict",
            "pairs": relative(pair_file),
            "pep": relative(combined_pep),
            "cds": relative(combined_cds),
            "ks_config": relative(config_file),
            "ks_result": relative(ks_result),
            "pair_count": str(len(pairs)),
            "status": "READY",
        }
    )

with summary_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(summary_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(summary_rows)

with ks_manifest.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(manifest_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(manifest_rows)

print(
    f"Prepared {len(manifest_rows)} Ks datasets with "
    f"{sum(int(row['pair_count']) for row in manifest_rows):,} "
    "unique nonself pairs."
)
PY

###############################################################################
# FINAL VALIDATION
###############################################################################

READY_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "READY" {
            count++
        }
        END {
            print count + 0
        }
    ' "${KS_MANIFEST}"
)"

TOTAL_PAIRS="$(
    awk -F $'\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "pair_count") {
                    pair_column = i
                }
            }
            next
        }
        {
            sum += $pair_column
        }
        END {
            print sum + 0
        }
    ' "${KS_MANIFEST}"
)"

if [[ "${READY_COUNT}" -ne 12 ]]; then
    echo "ERROR: Only ${READY_COUNT}/12 Ks datasets are ready." >&2
    exit 1
fi

if [[ "${TOTAL_PAIRS}" -le 0 ]]; then
    echo "ERROR: No Ks pairs were prepared." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36D1_COMPLETE.txt" <<EOF2
checkpoint=step36D1_prepare_WGDI_Ks_inputs
date=$(date --iso-8601=seconds)
comparisons_expected=12
comparisons_ready=${READY_COUNT}
parameter_set=strict
mg=25,25
pair_policy=unique_nonself_gene_pairs
self_pair_policy=removed
self_reverse_pair_policy=canonicalized
overlapping_block_pair_policy=deduplicated
alignment_software=mafft
total_pairs=${TOTAL_PAIRS}
manifest=11_wgdi/00_admin/step36D_ks_manifest.tsv
summary=11_wgdi/02_qc/ks_preparation/step36D1_ks_input_summary.tsv
status=PASS
next_step=step36D2_run_WGDI_Ks
EOF2

find \
    "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "============================================================"
echo "Step 36D1 Ks input summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36D Ks manifest"
echo "============================================================"

column -t -s $'\t' \
    "${KS_MANIFEST}"

echo
echo "============================================================"
echo "Step 36D1 checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36D1_COMPLETE.txt"
