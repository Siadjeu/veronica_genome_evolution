#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=jcvi34_sum
#SBATCH --output=10_synteny/logs/jcvi34_summary_%j.out
#SBATCH --error=10_synteny/logs/jcvi34_summary_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

STEP_DIR="${SYNTENY_DIR}/jcvi_step34"
COMPARISON_MANIFEST="${STEP_DIR}/manifests/jcvi_comparisons.tsv"
WORK_ROOT="${STEP_DIR}/comparisons"
STATUS_DIR="${STEP_DIR}/status"
TABLE_DIR="${STEP_DIR}/tables"

CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_jcvi_step34"

mkdir -p \
    "${TABLE_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${SYNTENY_DIR}/logs"

cd "${PROJECT_DIR}"

if [[ ! -s "${COMPARISON_MANIFEST}" ]]; then
    echo "ERROR: Missing comparison manifest:" >&2
    echo "${COMPARISON_MANIFEST}" >&2
    exit 1
fi

rm -f \
    "${TABLE_DIR}/jcvi_step34_results.tsv" \
    "${TABLE_DIR}/jcvi_step34_summary.tsv" \
    "${TABLE_DIR}/jcvi_step34_failed_comparisons.tsv" \
    "${TABLE_DIR}/jcvi_step34_anchor_files.tsv"

rm -f "${CHECKPOINT_DIR}"/*

python - \
    "${COMPARISON_MANIFEST}" \
    "${WORK_ROOT}" \
    "${STATUS_DIR}" \
    "${TABLE_DIR}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

comparison_manifest = Path(sys.argv[1])
work_root = Path(sys.argv[2])
status_dir = Path(sys.argv[3])
table_dir = Path(sys.argv[4])

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

if len(comparisons) != 28:
    raise SystemExit(
        f"ERROR: Expected 28 comparisons; found {len(comparisons)}."
    )


def read_key_value_file(path: Path) -> dict[str, str]:
    values = {}

    if not path.is_file():
        return values

    for line in path.read_text(
        encoding="utf-8",
        errors="replace",
    ).splitlines():
        if "=" not in line:
            continue

        key, value = line.split("=", 1)
        values[key.strip()] = value.strip()

    return values


results = []
anchor_file_rows = []

for comparison in comparisons:
    comparison_id = comparison["comparison_id"]
    work_dir = work_root / comparison_id

    complete_file = work_dir / "RUN_COMPLETE.txt"
    failed_file = work_dir / "RUN_FAILED.txt"
    status_file = status_dir / f"{comparison_id}.status.tsv"

    complete = read_key_value_file(
        complete_file
    )

    failed = read_key_value_file(
        failed_file
    )

    if complete.get("status") == "PASS":
        status = "PASS"
        reason = ""
    elif failed:
        status = "FAIL"
        reason = failed.get(
            "reason",
            "UNKNOWN_FAILURE",
        )
    else:
        status = "MISSING"
        reason = "NO_RUN_COMPLETION_RECORD"

    anchor_file = complete.get(
        "anchor_file",
        "",
    )

    anchor_exists = (
        bool(anchor_file)
        and Path(anchor_file).is_file()
        and Path(anchor_file).stat().st_size > 0
    )

    if status == "PASS" and not anchor_exists:
        status = "FAIL"
        reason = "RECORDED_ANCHOR_FILE_MISSING"

    result = {
        "task_id": comparison["task_id"],
        "comparison_id": comparison_id,
        "query_species": comparison[
            "query_species"
        ],
        "subject_species": comparison[
            "subject_species"
        ],
        "comparison_type": comparison[
            "comparison_type"
        ],
        "priority": comparison[
            "priority"
        ],
        "cscore": comparison[
            "cscore"
        ],
        "minspan": comparison[
            "minspan"
        ],
        "align_soft": complete.get(
            "align_soft",
            "",
        ),
        "anchor_file": anchor_file,
        "anchor_pair_count": complete.get(
            "anchor_pair_count",
            "0",
        ),
        "anchor_block_count": complete.get(
            "anchor_block_count",
            "0",
        ),
        "unique_query_anchor_genes": complete.get(
            "unique_query_anchor_genes",
            "0",
        ),
        "unique_subject_anchor_genes": complete.get(
            "unique_subject_anchor_genes",
            "0",
        ),
        "simple_anchor_file": complete.get(
            "simple_anchor_file",
            "",
        ),
        "simple_anchor_pair_count": complete.get(
            "simple_anchor_pair_count",
            "0",
        ),
        "simple_status": complete.get(
            "simple_status",
            "",
        ),
        "status": status,
        "reason": reason,
        "work_directory": str(work_dir),
        "status_file": str(status_file),
    }

    results.append(result)

    anchor_file_rows.append(
        {
            "comparison_id": comparison_id,
            "query_species": comparison[
                "query_species"
            ],
            "subject_species": comparison[
                "subject_species"
            ],
            "comparison_type": comparison[
                "comparison_type"
            ],
            "anchor_file": anchor_file,
            "anchor_file_exists": (
                "YES"
                if anchor_exists
                else "NO"
            ),
            "anchor_pair_count": result[
                "anchor_pair_count"
            ],
            "anchor_block_count": result[
                "anchor_block_count"
            ],
        }
    )

result_fields = [
    "task_id",
    "comparison_id",
    "query_species",
    "subject_species",
    "comparison_type",
    "priority",
    "cscore",
    "minspan",
    "align_soft",
    "anchor_file",
    "anchor_pair_count",
    "anchor_block_count",
    "unique_query_anchor_genes",
    "unique_subject_anchor_genes",
    "simple_anchor_file",
    "simple_anchor_pair_count",
    "simple_status",
    "status",
    "reason",
    "work_directory",
    "status_file",
]

result_file = (
    table_dir
    / "jcvi_step34_results.tsv"
)

with result_file.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=result_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(results)

non_pass = [
    row
    for row in results
    if row["status"] != "PASS"
]

failed_file = (
    table_dir
    / "jcvi_step34_failed_comparisons.tsv"
)

with failed_file.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=result_fields,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(non_pass)

anchor_file_table = (
    table_dir
    / "jcvi_step34_anchor_files.tsv"
)

with anchor_file_table.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "comparison_id",
            "query_species",
            "subject_species",
            "comparison_type",
            "anchor_file",
            "anchor_file_exists",
            "anchor_pair_count",
            "anchor_block_count",
        ],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(anchor_file_rows)

comparison_pass = sum(
    row["status"] == "PASS"
    for row in results
)

comparison_fail = sum(
    row["status"] == "FAIL"
    for row in results
)

comparison_missing = sum(
    row["status"] == "MISSING"
    for row in results
)

self_pass = sum(
    row["status"] == "PASS"
    and row["comparison_type"] == "self"
    for row in results
)

pairwise_pass = sum(
    row["status"] == "PASS"
    and row["comparison_type"] == "pairwise"
    for row in results
)

outgroup_pass = sum(
    row["status"] == "PASS"
    and row["comparison_type"] == "outgroup"
    for row in results
)

simple_pass = sum(
    row["simple_status"] == "PASS"
    for row in results
)

simple_failed = sum(
    row["simple_status"] == "OPTIONAL_SCREEN_FAILED"
    for row in results
)

summary_file = (
    table_dir
    / "jcvi_step34_summary.tsv"
)

with summary_file.open(
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
    writer.writerow(
        ["comparisons_total", len(results)]
    )
    writer.writerow(
        ["comparisons_pass", comparison_pass]
    )
    writer.writerow(
        ["comparisons_fail", comparison_fail]
    )
    writer.writerow(
        ["comparisons_missing", comparison_missing]
    )
    writer.writerow(
        ["self_comparisons_pass", self_pass]
    )
    writer.writerow(
        ["pairwise_comparisons_pass", pairwise_pass]
    )
    writer.writerow(
        ["outgroup_comparisons_pass", outgroup_pass]
    )
    writer.writerow(
        ["simple_screen_pass", simple_pass]
    )
    writer.writerow(
        ["simple_screen_optional_fail", simple_failed]
    )
    writer.writerow(
        [
            "total_anchor_pairs",
            sum(
                int(row["anchor_pair_count"])
                for row in results
            ),
        ]
    )
    writer.writerow(
        [
            "total_anchor_blocks",
            sum(
                int(row["anchor_block_count"])
                for row in results
            ),
        ]
    )
    writer.writerow(
        [
            "total_simple_anchor_pairs",
            sum(
                int(
                    row[
                        "simple_anchor_pair_count"
                    ]
                )
                for row in results
            ),
        ]
    )

print("Step 34 comparisons:", len(results))
print("PASS:", comparison_pass)
print("FAIL:", comparison_fail)
print("MISSING:", comparison_missing)
print("Self PASS:", self_pass)
print("Pairwise PASS:", pairwise_pass)
print("Outgroup PASS:", outgroup_pass)
print("Optional simple screen PASS:", simple_pass)
print("Optional simple screen failures:", simple_failed)

if non_pass:
    print(
        "Non-PASS comparisons:",
        ",".join(
            row["comparison_id"]
            for row in non_pass
        ),
    )

    raise SystemExit(
        "ERROR: Step 34 contains failed or missing comparisons."
    )
PY

echo
echo "============================================================"
echo "Step 34 summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/jcvi_step34_summary.tsv"

echo
echo "============================================================"
echo "Step 34 results"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/jcvi_step34_results.tsv"

COMPARISONS_PASS=$(
    awk -F'\t' '
        $1 == "comparisons_pass" {
            print $2
        }
    ' "${TABLE_DIR}/jcvi_step34_summary.tsv"
)

COMPARISONS_FAIL=$(
    awk -F'\t' '
        $1 == "comparisons_fail" {
            print $2
        }
    ' "${TABLE_DIR}/jcvi_step34_summary.tsv"
)

COMPARISONS_MISSING=$(
    awk -F'\t' '
        $1 == "comparisons_missing" {
            print $2
        }
    ' "${TABLE_DIR}/jcvi_step34_summary.tsv"
)

if [[ "${COMPARISONS_PASS}" -ne 28 ]]; then
    echo "ERROR: Expected 28 successful comparisons; observed ${COMPARISONS_PASS}." >&2
    exit 1
fi

if [[ "${COMPARISONS_FAIL}" -ne 0 ]]; then
    echo "ERROR: Expected 0 failed comparisons; observed ${COMPARISONS_FAIL}." >&2
    exit 1
fi

if [[ "${COMPARISONS_MISSING}" -ne 0 ]]; then
    echo "ERROR: Expected 0 missing comparisons; observed ${COMPARISONS_MISSING}." >&2
    exit 1
fi

cp -f \
    "${TABLE_DIR}/jcvi_step34_summary.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_step34_results.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_step34_failed_comparisons.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/jcvi_step34_anchor_files.tsv" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/JCVI_STEP34_COMPLETE.txt" <<EOF2
checkpoint=jcvi_self_and_pairwise_synteny
date=$(date --iso-8601=seconds)
comparisons_total=28
comparisons_pass=28
comparisons_fail=0
comparisons_missing=0
self_comparisons_pass=9
pairwise_comparisons_pass=16
outgroup_comparisons_pass=3
status=PASS
next_step=quantify_syntenic_depth_and_prepare_dotplots
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
echo "Step 34C completed successfully"
echo "============================================================"
echo "Checkpoint:"
echo "${CHECKPOINT_DIR}/JCVI_STEP34_COMPLETE.txt"
