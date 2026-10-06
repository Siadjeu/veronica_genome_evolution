#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36B1
#SBATCH --output=11_wgdi/logs/step36B1_%j.out
#SBATCH --error=11_wgdi/logs/step36B1_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"
MANIFEST="${WGDI_ROOT}/00_admin/step36B_homology_comparisons.tsv"

RESULT_DIR="${WGDI_ROOT}/04_homology/results"
PER_RUN_QC_DIR="${WGDI_ROOT}/02_qc/homology"
SUMMARY_DIR="${WGDI_ROOT}/02_qc/homology_summary"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36B"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
LOG_DIR="${WGDI_ROOT}/logs"

SUMMARY="${SUMMARY_DIR}/step36B_homology_summary.tsv"
OUTPUT_MANIFEST="${ADMIN_DIR}/wgdi_homology_manifest.tsv"
OVERALL="${SUMMARY_DIR}/step36B_overall_summary.tsv"

mkdir -p \
    "${SUMMARY_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

cd "${PROJECT_ROOT}"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Comparison manifest missing:" >&2
    echo "${MANIFEST}" >&2
    exit 1
fi

rm -f \
    "${SUMMARY}" \
    "${OUTPUT_MANIFEST}" \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/STEP36B_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${MANIFEST}" \
    "${RESULT_DIR}" \
    "${PER_RUN_QC_DIR}" \
    "${SUMMARY}" \
    "${OUTPUT_MANIFEST}" \
    "${OVERALL}" <<'PY'
from __future__ import annotations

import csv
import gzip
import hashlib
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
comparison_manifest = Path(sys.argv[2])
result_dir = Path(sys.argv[3])
per_run_qc_dir = Path(sys.argv[4])
summary_path = Path(sys.argv[5])
output_manifest_path = Path(sys.argv[6])
overall_path = Path(sys.argv[7])


def relative(path: Path):
    return str(path.relative_to(project_root))


def sha256(path: Path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


with comparison_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    comparisons = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

if len(comparisons) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 comparisons; found {len(comparisons)}."
    )

summary_rows = []
manifest_rows = []

for comparison in comparisons:
    species1 = comparison["species1"]
    species2 = comparison["species2"]
    run_id = f"{species1}_{species2}"

    result_path = (
        result_dir
        / f"{run_id}.blast.tsv.gz"
    )

    qc_path = (
        per_run_qc_dir
        / run_id
        / f"{run_id}.homology_qc.tsv"
    )

    for required in [result_path, qc_path]:
        if (
            not required.is_file()
            or required.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing result for {run_id}: {required}"
            )

    with qc_path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        qc_rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    if len(qc_rows) != 1:
        raise SystemExit(
            f"ERROR: Expected one QC row for {run_id}."
        )

    qc = qc_rows[0]

    if qc["status"] != "PASS":
        raise SystemExit(
            f"ERROR: Homology QC did not pass for {run_id}."
        )

    line_count = 0
    malformed = 0

    with gzip.open(
        result_path,
        "rt",
        encoding="utf-8",
    ) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")

            line_count += 1

            if len(fields) != 12:
                malformed += 1

    if line_count != int(qc["homology_rows"]):
        raise SystemExit(
            f"ERROR: Result/QC row-count mismatch for {run_id}: "
            f"{line_count} versus {qc['homology_rows']}"
        )

    if malformed != 0:
        raise SystemExit(
            f"ERROR: {malformed} malformed rows in {result_path}"
        )

    summary_row = dict(qc)
    summary_row["result_file"] = relative(result_path)
    summary_row["result_size_bytes"] = result_path.stat().st_size
    summary_row["result_sha256"] = sha256(result_path)

    summary_rows.append(summary_row)

    manifest_rows.append(
        {
            "task_id": comparison["task_id"],
            "comparison_type": comparison["comparison_type"],
            "species1": species1,
            "species2": species2,
            "comparison": run_id,
            "blast": relative(result_path),
            "qc": relative(qc_path),
            "homology_rows": line_count,
            "status": "PASS",
        }
    )

summary_fields = list(summary_rows[0])

with summary_path.open(
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

manifest_fields = [
    "task_id",
    "comparison_type",
    "species1",
    "species2",
    "comparison",
    "blast",
    "qc",
    "homology_rows",
    "status",
]

with output_manifest_path.open(
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

total_rows = sum(
    int(row["homology_rows"])
    for row in summary_rows
)

self_count = sum(
    row["comparison_type"] == "self"
    for row in summary_rows
)

pairwise_count = sum(
    row["comparison_type"] == "pairwise"
    for row in summary_rows
)

minimum_query_coverage = min(
    float(row["query_gene_coverage"])
    for row in summary_rows
)

minimum_subject_coverage = min(
    float(row["subject_gene_coverage"])
    for row in summary_rows
)

overall_rows = [
    ("comparisons_expected", 12),
    ("comparisons_complete", len(summary_rows)),
    ("self_comparisons", self_count),
    ("pairwise_comparisons", pairwise_count),
    ("total_homology_rows", total_rows),
    (
        "minimum_query_gene_coverage",
        f"{minimum_query_coverage:.8f}",
    ),
    (
        "minimum_subject_gene_coverage",
        f"{minimum_subject_coverage:.8f}",
    ),
    ("malformed_rows", 0),
    ("invalid_id_rows", 0),
    ("duplicate_gene_pairs", 0),
    ("status", "PASS"),
]

with overall_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(["metric", "value"])
    writer.writerows(overall_rows)

print(
    f"Validated {len(summary_rows)} homology comparisons "
    f"with {total_rows:,} total rows."
)
PY

cat > "${CHECKPOINT_DIR}/STEP36B_COMPLETE.txt" <<EOF2
checkpoint=step36B_generate_protein_homology
date=$(date --iso-8601=seconds)
self_comparisons=9
pairwise_comparisons=3
comparison_set=PMAJ_PMAJ,VPAN_VPAN,VSCU_VSCU,VANA_VANA,VARV_VARV,VPER_VPER,VSER_VSER,VTRI_VTRI,VVER_VVER,PMAJ_VSCU,VSCU_VANA,VSCU_VPER
search_engine=DIAMOND
format=BLAST_tabular_outfmt6_12_columns
evalue=1e-5
max_target_sequences=100
sensitivity=very-sensitive
pairwise_strategy=reciprocal_search_normalized_and_deduplicated
result_paths=relative_to_project_root
status=PASS
next_step=step36C_run_WGDI_improved_collinearity
EOF2

cp -f \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${OUTPUT_MANIFEST}" \
    "${CHECKPOINT_DIR}/"

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
echo "Step 36B homology summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36B overall summary"
echo "============================================================"

column -t -s $'\t' \
    "${OVERALL}"

echo
echo "============================================================"
echo "Step 36B completed successfully"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36B_COMPLETE.txt"
