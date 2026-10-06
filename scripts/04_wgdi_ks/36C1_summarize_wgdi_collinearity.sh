#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36C1
#SBATCH --output=11_wgdi/logs/step36C1_%j.out
#SBATCH --error=11_wgdi/logs/step36C1_%j.err

set -euo pipefail

###############################################################################
# PROJECT CONFIGURATION
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

MANIFEST="${WGDI_ROOT}/00_admin/step36C_collinearity_manifest.tsv"

RESULT_ROOT="${WGDI_ROOT}/05_collinearity/results"
QC_ROOT="${WGDI_ROOT}/02_qc/collinearity"

SUMMARY_DIR="${WGDI_ROOT}/02_qc/collinearity_summary"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36C"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
LOG_DIR="${WGDI_ROOT}/logs"

SUMMARY="${SUMMARY_DIR}/step36C_collinearity_summary.tsv"
COMPARISON_TABLE="${SUMMARY_DIR}/step36C_strict_sensitive_comparison.tsv"
OVERALL="${SUMMARY_DIR}/step36C_overall_summary.tsv"
DUPLICATE_TABLE="${SUMMARY_DIR}/step36C_repeated_pair_summary.tsv"
OUTPUT_MANIFEST="${ADMIN_DIR}/wgdi_collinearity_manifest.tsv"

mkdir -p \
    "${SUMMARY_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Step 36C manifest is missing:" >&2
    echo "${MANIFEST}" >&2
    exit 1
fi

EXPECTED_TASKS="$(
    awk -F $'\t' '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
    ' "${MANIFEST}"
)"

if [[ "${EXPECTED_TASKS}" -ne 24 ]]; then
    echo "ERROR: Expected 24 manifest tasks; found ${EXPECTED_TASKS}." >&2
    exit 1
fi

QC_FILE_COUNT="$(
    find \
        "${QC_ROOT}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name '*.collinearity_qc.tsv' \
        | wc -l
)"

if [[ "${QC_FILE_COUNT}" -ne 24 ]]; then
    echo "ERROR: Expected 24 QC files; found ${QC_FILE_COUNT}." >&2
    exit 1
fi

###############################################################################
# TARGETED CLEANUP
###############################################################################

rm -f \
    "${SUMMARY}" \
    "${COMPARISON_TABLE}" \
    "${OVERALL}" \
    "${DUPLICATE_TABLE}" \
    "${OUTPUT_MANIFEST}" \
    "${CHECKPOINT_DIR}/STEP36C_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# CONSOLIDATE AND VALIDATE
###############################################################################

python - \
    "${PROJECT_ROOT}" \
    "${MANIFEST}" \
    "${RESULT_ROOT}" \
    "${QC_ROOT}" \
    "${SUMMARY}" \
    "${COMPARISON_TABLE}" \
    "${OVERALL}" \
    "${DUPLICATE_TABLE}" \
    "${OUTPUT_MANIFEST}" <<'PY'
from __future__ import annotations

import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

project_root = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
result_root = Path(sys.argv[3])
qc_root = Path(sys.argv[4])
summary_path = Path(sys.argv[5])
comparison_path = Path(sys.argv[6])
overall_path = Path(sys.argv[7])
duplicate_path = Path(sys.argv[8])
output_manifest_path = Path(sys.argv[9])


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def read_single_tsv_row(path: Path) -> dict[str, str]:
    with path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    if len(rows) != 1:
        raise SystemExit(
            f"ERROR: Expected exactly one data row in {path}; "
            f"found {len(rows)}."
        )

    return rows[0]


with manifest_path.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    manifest_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

if len(manifest_rows) != 24:
    raise SystemExit(
        f"ERROR: Expected 24 manifest rows; "
        f"found {len(manifest_rows)}."
    )

observed_task_ids = [
    int(row["task_id"])
    for row in manifest_rows
]

if observed_task_ids != list(range(24)):
    raise SystemExit(
        "ERROR: Task IDs are not exactly 0 through 23."
    )

required_qc_fields = {
    "comparison",
    "comparison_type",
    "parameter_set",
    "mg",
    "block_count",
    "collinear_gene_pairs",
    "unique_collinear_gene_pairs",
    "duplicate_gene_pairs",
    "species1_block_gene_coverage",
    "species2_block_gene_coverage",
    "malformed_header_lines",
    "malformed_data_lines",
    "invalid_gene_pairs",
    "order_mismatches",
    "block_size_mismatches",
    "status",
}

summary_rows: list[dict[str, str]] = []
output_manifest_rows: list[dict[str, str]] = []

for manifest_row in manifest_rows:
    comparison = manifest_row["comparison"]
    parameter_set = manifest_row["parameter_set"]
    run_id = f"{comparison}_{parameter_set}"

    result_file = (
        result_root
        / run_id
        / f"{run_id}.collinearity"
    )

    qc_file = (
        qc_root
        / run_id
        / f"{run_id}.collinearity_qc.tsv"
    )

    for required_file in [result_file, qc_file]:
        if (
            not required_file.is_file()
            or required_file.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing or empty Step 36C output: "
                f"{required_file}"
            )

    qc = read_single_tsv_row(qc_file)

    missing_fields = required_qc_fields.difference(qc)

    if missing_fields:
        raise SystemExit(
            f"ERROR: Missing QC fields for {run_id}: "
            + ",".join(sorted(missing_fields))
        )

    if qc["comparison"] != comparison:
        raise SystemExit(
            f"ERROR: Comparison mismatch for {run_id}."
        )

    if qc["comparison_type"] != manifest_row["comparison_type"]:
        raise SystemExit(
            f"ERROR: Comparison-type mismatch for {run_id}."
        )

    if qc["parameter_set"] != parameter_set:
        raise SystemExit(
            f"ERROR: Parameter-set mismatch for {run_id}."
        )

    if qc["mg"] != manifest_row["mg"]:
        raise SystemExit(
            f"ERROR: mg mismatch for {run_id}: "
            f"{qc['mg']} versus {manifest_row['mg']}"
        )

    if qc["status"] != "PASS":
        raise SystemExit(
            f"ERROR: QC status is not PASS for {run_id}."
        )

    # These indicate malformed or inconsistent output and must remain zero.
    zero_required_fields = [
        "malformed_header_lines",
        "malformed_data_lines",
        "invalid_gene_pairs",
        "order_mismatches",
        "block_size_mismatches",
    ]

    for field in zero_required_fields:
        if int(qc[field]) != 0:
            raise SystemExit(
                f"ERROR: {field} is nonzero for {run_id}: "
                f"{qc[field]}"
            )

    block_count = int(qc["block_count"])
    total_pairs = int(qc["collinear_gene_pairs"])
    unique_pairs = int(qc["unique_collinear_gene_pairs"])
    repeated_assignments = int(qc["duplicate_gene_pairs"])

    if block_count <= 0:
        raise SystemExit(
            f"ERROR: No collinearity blocks for {run_id}."
        )

    if total_pairs <= 0:
        raise SystemExit(
            f"ERROR: No collinear gene pairs for {run_id}."
        )

    if unique_pairs <= 0:
        raise SystemExit(
            f"ERROR: No unique collinear pairs for {run_id}."
        )

    if unique_pairs > total_pairs:
        raise SystemExit(
            f"ERROR: Unique pair count exceeds total pair count "
            f"for {run_id}."
        )

    expected_repeated_assignments = total_pairs - unique_pairs

    if repeated_assignments != expected_repeated_assignments:
        raise SystemExit(
            f"ERROR: Repeated-pair count is inconsistent for "
            f"{run_id}: reported={repeated_assignments}, "
            f"expected={expected_repeated_assignments}."
        )

    repeated_fraction = (
        repeated_assignments / total_pairs
        if total_pairs > 0
        else 0.0
    )

    summary_row = dict(qc)
    summary_row["repeated_pair_fraction"] = (
        f"{repeated_fraction:.8f}"
    )
    summary_row["result_file"] = relative(result_file)
    summary_row["result_size_bytes"] = str(
        result_file.stat().st_size
    )
    summary_row["result_sha256"] = sha256(result_file)
    summary_row["qc_file"] = relative(qc_file)

    summary_rows.append(summary_row)

    output_manifest_rows.append(
        {
            "task_id": manifest_row["task_id"],
            "comparison_type": manifest_row[
                "comparison_type"
            ],
            "species1": manifest_row["species1"],
            "species2": manifest_row["species2"],
            "comparison": comparison,
            "parameter_set": parameter_set,
            "mg": manifest_row["mg"],
            "collinearity": relative(result_file),
            "qc": relative(qc_file),
            "block_count": str(block_count),
            "collinear_gene_pairs": str(total_pairs),
            "unique_collinear_gene_pairs": str(unique_pairs),
            "repeated_pair_assignments": str(
                repeated_assignments
            ),
            "repeated_pair_fraction": (
                f"{repeated_fraction:.8f}"
            ),
            "status": "PASS",
        }
    )

###############################################################################
# WRITE FULL SUMMARY
###############################################################################

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

###############################################################################
# WRITE REPEATED-PAIR SUMMARY
###############################################################################

duplicate_rows = []

for row in summary_rows:
    total_pairs = int(row["collinear_gene_pairs"])
    unique_pairs = int(row["unique_collinear_gene_pairs"])
    repeated_assignments = int(row["duplicate_gene_pairs"])

    duplicate_rows.append(
        {
            "comparison": row["comparison"],
            "comparison_type": row["comparison_type"],
            "parameter_set": row["parameter_set"],
            "mg": row["mg"],
            "total_collinear_pair_assignments": str(
                total_pairs
            ),
            "unique_collinear_gene_pairs": str(
                unique_pairs
            ),
            "repeated_pair_assignments": str(
                repeated_assignments
            ),
            "repeated_pair_fraction": row[
                "repeated_pair_fraction"
            ],
            "interpretation": (
                "overlapping_or_reused_across_WGDI_blocks"
                if repeated_assignments > 0
                else "no_repeated_pair_assignments"
            ),
            "status": "PASS",
        }
    )

with duplicate_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(duplicate_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(duplicate_rows)

###############################################################################
# STRICT VERSUS SENSITIVE COMPARISON
###############################################################################

grouped: dict[str, dict[str, dict[str, str]]] = defaultdict(dict)

for row in summary_rows:
    grouped[row["comparison"]][
        row["parameter_set"]
    ] = row

if len(grouped) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 comparisons; found {len(grouped)}."
    )

comparison_rows: list[dict[str, str]] = []

for comparison in sorted(grouped):
    parameter_rows = grouped[comparison]

    if set(parameter_rows) != {"strict", "sensitive"}:
        raise SystemExit(
            f"ERROR: Missing strict or sensitive result for "
            f"{comparison}."
        )

    strict = parameter_rows["strict"]
    sensitive = parameter_rows["sensitive"]

    strict_blocks = int(strict["block_count"])
    sensitive_blocks = int(sensitive["block_count"])

    strict_pairs = int(strict["collinear_gene_pairs"])
    sensitive_pairs = int(
        sensitive["collinear_gene_pairs"]
    )

    strict_unique_pairs = int(
        strict["unique_collinear_gene_pairs"]
    )
    sensitive_unique_pairs = int(
        sensitive["unique_collinear_gene_pairs"]
    )

    strict_repeated = int(strict["duplicate_gene_pairs"])
    sensitive_repeated = int(
        sensitive["duplicate_gene_pairs"]
    )

    block_ratio = sensitive_blocks / strict_blocks
    pair_ratio = sensitive_pairs / strict_pairs
    unique_pair_ratio = (
        sensitive_unique_pairs / strict_unique_pairs
    )

    block_percent_change = (
        (sensitive_blocks - strict_blocks)
        / strict_blocks
        * 100
    )

    pair_percent_change = (
        (sensitive_pairs - strict_pairs)
        / strict_pairs
        * 100
    )

    unique_pair_percent_change = (
        (sensitive_unique_pairs - strict_unique_pairs)
        / strict_unique_pairs
        * 100
    )

    strict_species1_coverage = float(
        strict["species1_block_gene_coverage"]
    )
    sensitive_species1_coverage = float(
        sensitive["species1_block_gene_coverage"]
    )

    strict_species2_coverage = float(
        strict["species2_block_gene_coverage"]
    )
    sensitive_species2_coverage = float(
        sensitive["species2_block_gene_coverage"]
    )

    maximum_absolute_percent_change = max(
        abs(block_percent_change),
        abs(pair_percent_change),
        abs(unique_pair_percent_change),
    )

    if maximum_absolute_percent_change >= 50:
        review_flag = "HIGH_SENSITIVITY"
    elif maximum_absolute_percent_change >= 20:
        review_flag = "MODERATE_SENSITIVITY"
    else:
        review_flag = "STABLE"

    comparison_rows.append(
        {
            "comparison": comparison,
            "comparison_type": strict[
                "comparison_type"
            ],
            "strict_mg": strict["mg"],
            "sensitive_mg": sensitive["mg"],
            "strict_blocks": str(strict_blocks),
            "sensitive_blocks": str(sensitive_blocks),
            "block_ratio_sensitive_to_strict": (
                f"{block_ratio:.8f}"
            ),
            "block_percent_change": (
                f"{block_percent_change:.4f}"
            ),
            "strict_total_pair_assignments": str(
                strict_pairs
            ),
            "sensitive_total_pair_assignments": str(
                sensitive_pairs
            ),
            "pair_assignment_ratio": (
                f"{pair_ratio:.8f}"
            ),
            "pair_assignment_percent_change": (
                f"{pair_percent_change:.4f}"
            ),
            "strict_unique_gene_pairs": str(
                strict_unique_pairs
            ),
            "sensitive_unique_gene_pairs": str(
                sensitive_unique_pairs
            ),
            "unique_pair_ratio": (
                f"{unique_pair_ratio:.8f}"
            ),
            "unique_pair_percent_change": (
                f"{unique_pair_percent_change:.4f}"
            ),
            "strict_repeated_pair_assignments": str(
                strict_repeated
            ),
            "sensitive_repeated_pair_assignments": str(
                sensitive_repeated
            ),
            "strict_repeated_pair_fraction": strict[
                "repeated_pair_fraction"
            ],
            "sensitive_repeated_pair_fraction": sensitive[
                "repeated_pair_fraction"
            ],
            "strict_species1_coverage": (
                f"{strict_species1_coverage:.8f}"
            ),
            "sensitive_species1_coverage": (
                f"{sensitive_species1_coverage:.8f}"
            ),
            "species1_coverage_change": (
                f"{sensitive_species1_coverage - strict_species1_coverage:.8f}"
            ),
            "strict_species2_coverage": (
                f"{strict_species2_coverage:.8f}"
            ),
            "sensitive_species2_coverage": (
                f"{sensitive_species2_coverage:.8f}"
            ),
            "species2_coverage_change": (
                f"{sensitive_species2_coverage - strict_species2_coverage:.8f}"
            ),
            "review_flag": review_flag,
        }
    )

with comparison_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(comparison_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(comparison_rows)

###############################################################################
# WRITE OUTPUT MANIFEST
###############################################################################

manifest_fields = [
    "task_id",
    "comparison_type",
    "species1",
    "species2",
    "comparison",
    "parameter_set",
    "mg",
    "collinearity",
    "qc",
    "block_count",
    "collinear_gene_pairs",
    "unique_collinear_gene_pairs",
    "repeated_pair_assignments",
    "repeated_pair_fraction",
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
    writer.writerows(output_manifest_rows)

###############################################################################
# OVERALL SUMMARY
###############################################################################

strict_rows = [
    row
    for row in summary_rows
    if row["parameter_set"] == "strict"
]

sensitive_rows = [
    row
    for row in summary_rows
    if row["parameter_set"] == "sensitive"
]

stable_count = sum(
    row["review_flag"] == "STABLE"
    for row in comparison_rows
)

moderate_count = sum(
    row["review_flag"] == "MODERATE_SENSITIVITY"
    for row in comparison_rows
)

high_count = sum(
    row["review_flag"] == "HIGH_SENSITIVITY"
    for row in comparison_rows
)

strict_total_assignments = sum(
    int(row["collinear_gene_pairs"])
    for row in strict_rows
)

sensitive_total_assignments = sum(
    int(row["collinear_gene_pairs"])
    for row in sensitive_rows
)

strict_unique_pairs = sum(
    int(row["unique_collinear_gene_pairs"])
    for row in strict_rows
)

sensitive_unique_pairs = sum(
    int(row["unique_collinear_gene_pairs"])
    for row in sensitive_rows
)

strict_repeated_assignments = sum(
    int(row["duplicate_gene_pairs"])
    for row in strict_rows
)

sensitive_repeated_assignments = sum(
    int(row["duplicate_gene_pairs"])
    for row in sensitive_rows
)

overall_rows = [
    ("tasks_expected", "24"),
    ("tasks_complete", str(len(summary_rows))),
    ("comparisons_expected", "12"),
    ("comparisons_complete", str(len(comparison_rows))),
    ("strict_runs", str(len(strict_rows))),
    ("sensitive_runs", str(len(sensitive_rows))),
    (
        "strict_total_blocks",
        str(
            sum(
                int(row["block_count"])
                for row in strict_rows
            )
        ),
    ),
    (
        "sensitive_total_blocks",
        str(
            sum(
                int(row["block_count"])
                for row in sensitive_rows
            )
        ),
    ),
    (
        "strict_total_pair_assignments",
        str(strict_total_assignments),
    ),
    (
        "sensitive_total_pair_assignments",
        str(sensitive_total_assignments),
    ),
    (
        "strict_unique_collinear_gene_pairs",
        str(strict_unique_pairs),
    ),
    (
        "sensitive_unique_collinear_gene_pairs",
        str(sensitive_unique_pairs),
    ),
    (
        "strict_repeated_pair_assignments",
        str(strict_repeated_assignments),
    ),
    (
        "sensitive_repeated_pair_assignments",
        str(sensitive_repeated_assignments),
    ),
    (
        "strict_repeated_pair_fraction",
        f"{strict_repeated_assignments / strict_total_assignments:.8f}",
    ),
    (
        "sensitive_repeated_pair_fraction",
        f"{sensitive_repeated_assignments / sensitive_total_assignments:.8f}",
    ),
    ("stable_comparisons", str(stable_count)),
    ("moderate_sensitivity_comparisons", str(moderate_count)),
    ("high_sensitivity_comparisons", str(high_count)),
    ("malformed_header_lines", "0"),
    ("malformed_data_lines", "0"),
    ("invalid_gene_pairs", "0"),
    ("order_mismatches", "0"),
    ("block_size_mismatches", "0"),
    (
        "repeated_pair_policy",
        "retained_and_reported_as_overlapping_WGDI_block_assignments",
    ),
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
    "Step 36C consolidation completed: "
    f"{len(summary_rows)} runs and "
    f"{len(comparison_rows)} comparisons validated."
)
PY

###############################################################################
# FINAL CHECKPOINT
###############################################################################

TASKS_COMPLETE="$(
    awk -F $'\t' '
        $1 == "tasks_complete" {
            print $2
        }
    ' "${OVERALL}"
)"

STATUS="$(
    awk -F $'\t' '
        $1 == "status" {
            print $2
        }
    ' "${OVERALL}"
)"

if [[ "${TASKS_COMPLETE}" -ne 24 ]]; then
    echo "ERROR: Only ${TASKS_COMPLETE}/24 tasks were summarized." >&2
    exit 1
fi

if [[ "${STATUS}" != "PASS" ]]; then
    echo "ERROR: Consolidated status is ${STATUS}." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36C_COMPLETE.txt" <<EOF2
checkpoint=step36C_WGDI_improved_collinearity
date=$(date --iso-8601=seconds)
comparisons=12
self_comparisons=9
pairwise_comparisons=3
parameter_sets=strict,sensitive
strict_mg=25,25
sensitive_mg=50,50
evalue=1e-5
score=100
grading=50,30,25
pvalue=1
repeat_number=20
position_parameter_key=positon
position=order
tasks_expected=24
tasks_complete=24
repeated_pair_policy=retained_and_reported_as_overlapping_WGDI_block_assignments
comparison_summary=11_wgdi/02_qc/collinearity_summary/step36C_strict_sensitive_comparison.tsv
repeated_pair_summary=11_wgdi/02_qc/collinearity_summary/step36C_repeated_pair_summary.tsv
overall_summary=11_wgdi/02_qc/collinearity_summary/step36C_overall_summary.tsv
output_manifest=11_wgdi/00_admin/wgdi_collinearity_manifest.tsv
status=PASS
next_step=review_parameter_stability_and_prepare_Ks
EOF2

cp -f \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${COMPARISON_TABLE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${DUPLICATE_TABLE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${OUTPUT_MANIFEST}" \
    "${CHECKPOINT_DIR}/"

find \
    "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# REPORT
###############################################################################

echo
echo "============================================================"
echo "Step 36C overall summary"
echo "============================================================"

column -t -s $'\t' \
    "${OVERALL}"

echo
echo "============================================================"
echo "Strict versus sensitive comparison"
echo "============================================================"

column -t -s $'\t' \
    "${COMPARISON_TABLE}"

echo
echo "============================================================"
echo "Repeated-pair summary"
echo "============================================================"

column -t -s $'\t' \
    "${DUPLICATE_TABLE}"

echo
echo "============================================================"
echo "Step 36C checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36C_COMPLETE.txt"

echo
echo "Completed: $(date --iso-8601=seconds)"
