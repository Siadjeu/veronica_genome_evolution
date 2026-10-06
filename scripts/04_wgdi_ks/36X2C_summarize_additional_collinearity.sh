#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=wgdi36X2C
#SBATCH --output=11_wgdi/logs/step36X2C_%j.out
#SBATCH --error=11_wgdi/logs/step36X2C_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

MANIFEST="${WGDI_ROOT}/00_admin/step36X2_collinearity_manifest.tsv"

RESULT_ROOT="${WGDI_ROOT}/05_collinearity/results"
QC_ROOT="${WGDI_ROOT}/02_qc/collinearity"

SUMMARY_DIR="${WGDI_ROOT}/02_qc/additional_collinearity_summary"

SUMMARY="${SUMMARY_DIR}/step36X2_collinearity_summary.tsv"
COMPARISON_TABLE="${SUMMARY_DIR}/step36X2_strict_sensitive_comparison.tsv"
OVERALL="${SUMMARY_DIR}/step36X2_overall_summary.tsv"
DUPLICATE_TABLE="${SUMMARY_DIR}/step36X2_repeated_pair_summary.tsv"

OUTPUT_MANIFEST="${WGDI_ROOT}/00_admin/step36X2_wgdi_collinearity_manifest.tsv"

CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X2C"

mkdir -p \
    "${SUMMARY_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

###############################################################################
# INPUT VALIDATION
###############################################################################

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Missing manifest: ${MANIFEST}" >&2
    exit 1
fi

###############################################################################
# SUMMARIZE
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


project_root = Path(sys.argv[1]).resolve()
manifest_path = Path(sys.argv[2]).resolve()
result_root = Path(sys.argv[3]).resolve()
qc_root = Path(sys.argv[4]).resolve()

summary_path = Path(sys.argv[5])
comparison_path = Path(sys.argv[6])
overall_path = Path(sys.argv[7])
duplicate_path = Path(sys.argv[8])
output_manifest_path = Path(sys.argv[9])


def relative(path: Path) -> str:
    return str(
        path.resolve().relative_to(
            project_root
        )
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(
                1024 * 1024
            ),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


def read_single_tsv_row(
    path: Path,
) -> dict[str, str]:

    with path.open(
        newline="",
        encoding="utf-8",
    ) as handle:

        rows = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    if len(rows) != 1:
        raise SystemExit(
            f"ERROR: Expected exactly one QC row in {path}; "
            f"found {len(rows)}."
        )

    return rows[0]


###############################################################################
# READ MANIFEST
###############################################################################

with manifest_path.open(
    newline="",
    encoding="utf-8",
) as handle:

    manifest_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )


if len(manifest_rows) != 6:
    raise SystemExit(
        f"ERROR: Expected 6 manifest tasks; "
        f"found {len(manifest_rows)}."
    )


task_ids = [
    int(row["task_id"])
    for row in manifest_rows
]

if task_ids != list(range(6)):
    raise SystemExit(
        f"ERROR: Task IDs must be exactly 0 through 5; "
        f"found {task_ids}."
    )


###############################################################################
# REQUIRED QC FIELDS
###############################################################################

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


summary_rows = []
output_manifest_rows = []


###############################################################################
# VALIDATE EACH OF SIX RUNS
###############################################################################

for manifest_row in manifest_rows:

    comparison = manifest_row[
        "comparison"
    ]

    parameter_set = manifest_row[
        "parameter_set"
    ]

    run_id = (
        f"{comparison}_{parameter_set}"
    )

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

    for required_file in [
        result_file,
        qc_file,
    ]:

        if (
            not required_file.is_file()
            or required_file.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing or empty output: "
                f"{required_file}"
            )

    qc = read_single_tsv_row(
        qc_file
    )

    missing = (
        required_qc_fields
        - set(qc)
    )

    if missing:
        raise SystemExit(
            f"ERROR: Missing QC fields for {run_id}: "
            + ",".join(
                sorted(missing)
            )
        )

    if qc["comparison"] != comparison:
        raise SystemExit(
            f"ERROR: Comparison mismatch for {run_id}."
        )

    if (
        qc["comparison_type"]
        != manifest_row["comparison_type"]
    ):
        raise SystemExit(
            f"ERROR: Comparison type mismatch for {run_id}."
        )

    if (
        qc["parameter_set"]
        != parameter_set
    ):
        raise SystemExit(
            f"ERROR: Parameter-set mismatch for {run_id}."
        )

    if qc["mg"] != manifest_row["mg"]:
        raise SystemExit(
            f"ERROR: mg mismatch for {run_id}."
        )

    if qc["status"] != "PASS":
        raise SystemExit(
            f"ERROR: QC status is not PASS for {run_id}."
        )

    zero_required = [
        "malformed_header_lines",
        "malformed_data_lines",
        "invalid_gene_pairs",
        "order_mismatches",
        "block_size_mismatches",
    ]

    for field in zero_required:

        if int(qc[field]) != 0:
            raise SystemExit(
                f"ERROR: {field}={qc[field]} "
                f"for {run_id}."
            )

    blocks = int(
        qc["block_count"]
    )

    total_pairs = int(
        qc[
            "collinear_gene_pairs"
        ]
    )

    unique_pairs = int(
        qc[
            "unique_collinear_gene_pairs"
        ]
    )

    repeated = int(
        qc[
            "duplicate_gene_pairs"
        ]
    )

    if blocks <= 0:
        raise SystemExit(
            f"ERROR: No blocks for {run_id}."
        )

    if total_pairs <= 0:
        raise SystemExit(
            f"ERROR: No collinear pairs for {run_id}."
        )

    if unique_pairs <= 0:
        raise SystemExit(
            f"ERROR: No unique pairs for {run_id}."
        )

    if unique_pairs > total_pairs:
        raise SystemExit(
            f"ERROR: unique > total pairs for {run_id}."
        )

    expected_repeated = (
        total_pairs
        - unique_pairs
    )

    if repeated != expected_repeated:
        raise SystemExit(
            f"ERROR: repeated pair count mismatch "
            f"for {run_id}: "
            f"{repeated} vs {expected_repeated}."
        )

    repeated_fraction = (
        repeated / total_pairs
    )

    summary_row = dict(qc)

    summary_row[
        "repeated_pair_fraction"
    ] = f"{repeated_fraction:.8f}"

    summary_row[
        "result_file"
    ] = relative(
        result_file
    )

    summary_row[
        "result_size_bytes"
    ] = str(
        result_file.stat().st_size
    )

    summary_row[
        "result_sha256"
    ] = sha256(
        result_file
    )

    summary_row[
        "qc_file"
    ] = relative(
        qc_file
    )

    summary_rows.append(
        summary_row
    )

    output_manifest_rows.append(
        {
            "task_id": manifest_row[
                "task_id"
            ],
            "comparison_type": manifest_row[
                "comparison_type"
            ],
            "species1": manifest_row[
                "species1"
            ],
            "species2": manifest_row[
                "species2"
            ],
            "comparison": comparison,
            "parameter_set": parameter_set,
            "mg": manifest_row["mg"],
            "collinearity": relative(
                result_file
            ),
            "qc": relative(
                qc_file
            ),
            "block_count": str(
                blocks
            ),
            "collinear_gene_pairs": str(
                total_pairs
            ),
            "unique_collinear_gene_pairs": str(
                unique_pairs
            ),
            "repeated_pair_assignments": str(
                repeated
            ),
            "repeated_pair_fraction": (
                f"{repeated_fraction:.8f}"
            ),
            "status": "PASS",
        }
    )


###############################################################################
# FULL SUMMARY
###############################################################################

with summary_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=list(
            summary_rows[0]
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        summary_rows
    )


###############################################################################
# REPEATED-PAIR SUMMARY
###############################################################################

duplicate_rows = []

for row in summary_rows:

    total_pairs = int(
        row[
            "collinear_gene_pairs"
        ]
    )

    unique_pairs = int(
        row[
            "unique_collinear_gene_pairs"
        ]
    )

    repeated = int(
        row[
            "duplicate_gene_pairs"
        ]
    )

    duplicate_rows.append(
        {
            "comparison": row[
                "comparison"
            ],
            "comparison_type": row[
                "comparison_type"
            ],
            "parameter_set": row[
                "parameter_set"
            ],
            "mg": row["mg"],
            "total_collinear_pair_assignments": (
                total_pairs
            ),
            "unique_collinear_gene_pairs": (
                unique_pairs
            ),
            "repeated_pair_assignments": (
                repeated
            ),
            "repeated_pair_fraction": row[
                "repeated_pair_fraction"
            ],
            "interpretation": (
                "overlapping_or_reused_across_WGDI_blocks"
                if repeated > 0
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
        fieldnames=list(
            duplicate_rows[0]
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        duplicate_rows
    )


###############################################################################
# STRICT VS SENSITIVE
###############################################################################

grouped = defaultdict(dict)

for row in summary_rows:

    grouped[
        row["comparison"]
    ][
        row["parameter_set"]
    ] = row


if len(grouped) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 comparisons; "
        f"found {len(grouped)}."
    )


comparison_rows = []


for comparison in sorted(grouped):

    parameter_rows = grouped[
        comparison
    ]

    if set(parameter_rows) != {
        "strict",
        "sensitive",
    }:
        raise SystemExit(
            f"ERROR: {comparison} lacks strict/sensitive pair."
        )

    strict = parameter_rows[
        "strict"
    ]

    sensitive = parameter_rows[
        "sensitive"
    ]

    strict_blocks = int(
        strict[
            "block_count"
        ]
    )

    sensitive_blocks = int(
        sensitive[
            "block_count"
        ]
    )

    strict_pairs = int(
        strict[
            "collinear_gene_pairs"
        ]
    )

    sensitive_pairs = int(
        sensitive[
            "collinear_gene_pairs"
        ]
    )

    strict_unique = int(
        strict[
            "unique_collinear_gene_pairs"
        ]
    )

    sensitive_unique = int(
        sensitive[
            "unique_collinear_gene_pairs"
        ]
    )

    strict_repeated = int(
        strict[
            "duplicate_gene_pairs"
        ]
    )

    sensitive_repeated = int(
        sensitive[
            "duplicate_gene_pairs"
        ]
    )

    block_ratio = (
        sensitive_blocks
        / strict_blocks
    )

    pair_ratio = (
        sensitive_pairs
        / strict_pairs
    )

    unique_ratio = (
        sensitive_unique
        / strict_unique
    )

    block_change = (
        (
            sensitive_blocks
            - strict_blocks
        )
        / strict_blocks
        * 100
    )

    pair_change = (
        (
            sensitive_pairs
            - strict_pairs
        )
        / strict_pairs
        * 100
    )

    unique_change = (
        (
            sensitive_unique
            - strict_unique
        )
        / strict_unique
        * 100
    )

    maximum_change = max(
        abs(block_change),
        abs(pair_change),
        abs(unique_change),
    )

    #
    # Exact original Step 36C1 thresholds:
    #
    # <20%   = STABLE
    # 20-50% = MODERATE_SENSITIVITY
    # >=50%  = HIGH_SENSITIVITY
    #
    if maximum_change >= 50:
        review_flag = (
            "HIGH_SENSITIVITY"
        )

    elif maximum_change >= 20:
        review_flag = (
            "MODERATE_SENSITIVITY"
        )

    else:
        review_flag = (
            "STABLE"
        )

    s1_strict = float(
        strict[
            "species1_block_gene_coverage"
        ]
    )

    s1_sensitive = float(
        sensitive[
            "species1_block_gene_coverage"
        ]
    )

    s2_strict = float(
        strict[
            "species2_block_gene_coverage"
        ]
    )

    s2_sensitive = float(
        sensitive[
            "species2_block_gene_coverage"
        ]
    )

    comparison_rows.append(
        {
            "comparison": comparison,
            "comparison_type": strict[
                "comparison_type"
            ],
            "strict_mg": strict["mg"],
            "sensitive_mg": sensitive[
                "mg"
            ],
            "strict_blocks": strict_blocks,
            "sensitive_blocks": sensitive_blocks,
            "block_ratio_sensitive_to_strict": (
                f"{block_ratio:.8f}"
            ),
            "block_percent_change": (
                f"{block_change:.4f}"
            ),
            "strict_total_pair_assignments": (
                strict_pairs
            ),
            "sensitive_total_pair_assignments": (
                sensitive_pairs
            ),
            "pair_assignment_ratio": (
                f"{pair_ratio:.8f}"
            ),
            "pair_assignment_percent_change": (
                f"{pair_change:.4f}"
            ),
            "strict_unique_gene_pairs": (
                strict_unique
            ),
            "sensitive_unique_gene_pairs": (
                sensitive_unique
            ),
            "unique_pair_ratio": (
                f"{unique_ratio:.8f}"
            ),
            "unique_pair_percent_change": (
                f"{unique_change:.4f}"
            ),
            "strict_repeated_pair_assignments": (
                strict_repeated
            ),
            "sensitive_repeated_pair_assignments": (
                sensitive_repeated
            ),
            "strict_repeated_pair_fraction": strict[
                "repeated_pair_fraction"
            ],
            "sensitive_repeated_pair_fraction": sensitive[
                "repeated_pair_fraction"
            ],
            "strict_species1_coverage": (
                f"{s1_strict:.8f}"
            ),
            "sensitive_species1_coverage": (
                f"{s1_sensitive:.8f}"
            ),
            "species1_coverage_change": (
                f"{s1_sensitive-s1_strict:.8f}"
            ),
            "strict_species2_coverage": (
                f"{s2_strict:.8f}"
            ),
            "sensitive_species2_coverage": (
                f"{s2_sensitive:.8f}"
            ),
            "species2_coverage_change": (
                f"{s2_sensitive-s2_strict:.8f}"
            ),
            "maximum_absolute_percent_change": (
                f"{maximum_change:.4f}"
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
        fieldnames=list(
            comparison_rows[0]
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        comparison_rows
    )


###############################################################################
# OUTPUT MANIFEST
###############################################################################

with output_manifest_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=list(
            output_manifest_rows[0]
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        output_manifest_rows
    )


###############################################################################
# OVERALL SUMMARY
###############################################################################

strict_rows = [
    row
    for row in summary_rows
    if row["parameter_set"]
    == "strict"
]

sensitive_rows = [
    row
    for row in summary_rows
    if row["parameter_set"]
    == "sensitive"
]


stable_count = sum(
    row["review_flag"]
    == "STABLE"
    for row in comparison_rows
)

moderate_count = sum(
    row["review_flag"]
    == "MODERATE_SENSITIVITY"
    for row in comparison_rows
)

high_count = sum(
    row["review_flag"]
    == "HIGH_SENSITIVITY"
    for row in comparison_rows
)


def total(
    rows,
    field,
):
    return sum(
        int(row[field])
        for row in rows
    )


strict_total_assignments = total(
    strict_rows,
    "collinear_gene_pairs",
)

sensitive_total_assignments = total(
    sensitive_rows,
    "collinear_gene_pairs",
)

strict_unique_pairs = total(
    strict_rows,
    "unique_collinear_gene_pairs",
)

sensitive_unique_pairs = total(
    sensitive_rows,
    "unique_collinear_gene_pairs",
)

strict_repeated = total(
    strict_rows,
    "duplicate_gene_pairs",
)

sensitive_repeated = total(
    sensitive_rows,
    "duplicate_gene_pairs",
)


overall_rows = [
    ("tasks_expected", "6"),
    ("tasks_complete", str(len(summary_rows))),
    ("comparisons_expected", "3"),
    ("comparisons_complete", str(len(comparison_rows))),
    ("strict_runs", str(len(strict_rows))),
    ("sensitive_runs", str(len(sensitive_rows))),
    (
        "strict_total_blocks",
        str(
            total(
                strict_rows,
                "block_count",
            )
        ),
    ),
    (
        "sensitive_total_blocks",
        str(
            total(
                sensitive_rows,
                "block_count",
            )
        ),
    ),
    (
        "strict_total_pair_assignments",
        str(
            strict_total_assignments
        ),
    ),
    (
        "sensitive_total_pair_assignments",
        str(
            sensitive_total_assignments
        ),
    ),
    (
        "strict_unique_collinear_gene_pairs",
        str(
            strict_unique_pairs
        ),
    ),
    (
        "sensitive_unique_collinear_gene_pairs",
        str(
            sensitive_unique_pairs
        ),
    ),
    (
        "strict_repeated_pair_assignments",
        str(
            strict_repeated
        ),
    ),
    (
        "sensitive_repeated_pair_assignments",
        str(
            sensitive_repeated
        ),
    ),
    (
        "strict_repeated_pair_fraction",
        f"{strict_repeated / strict_total_assignments:.8f}",
    ),
    (
        "sensitive_repeated_pair_fraction",
        f"{sensitive_repeated / sensitive_total_assignments:.8f}",
    ),
    (
        "stable_comparisons",
        str(stable_count),
    ),
    (
        "moderate_sensitivity_comparisons",
        str(moderate_count),
    ),
    (
        "high_sensitivity_comparisons",
        str(high_count),
    ),
    (
        "malformed_header_lines",
        "0",
    ),
    (
        "malformed_data_lines",
        "0",
    ),
    (
        "invalid_gene_pairs",
        "0",
    ),
    (
        "order_mismatches",
        "0",
    ),
    (
        "block_size_mismatches",
        "0",
    ),
    (
        "repeated_pair_policy",
        "retained_and_reported_as_overlapping_WGDI_block_assignments",
    ),
    (
        "status",
        "PASS",
    ),
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

    writer.writerow(
        [
            "metric",
            "value",
        ]
    )

    writer.writerows(
        overall_rows
    )


print(
    f"PASS: {len(summary_rows)}/6 runs; "
    f"{len(comparison_rows)}/3 comparisons."
)
PY

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36X2C_COMPLETE.txt" <<EOF2
checkpoint=step36X2C_additional_WGDI_collinearity
date=$(date --iso-8601=seconds)
comparisons=3
pairwise_comparisons=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
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
tasks_expected=6
tasks_complete=6
repeated_pair_policy=retained_and_reported_as_overlapping_WGDI_block_assignments
comparison_summary=11_wgdi/02_qc/additional_collinearity_summary/step36X2_strict_sensitive_comparison.tsv
repeated_pair_summary=11_wgdi/02_qc/additional_collinearity_summary/step36X2_repeated_pair_summary.tsv
overall_summary=11_wgdi/02_qc/additional_collinearity_summary/step36X2_overall_summary.tsv
output_manifest=11_wgdi/00_admin/step36X2_wgdi_collinearity_manifest.tsv
status=PASS
next_step=step36X3_prepare_additional_Ks_inputs
EOF2

cp -f \
    "${SUMMARY}" \
    "${COMPARISON_TABLE}" \
    "${DUPLICATE_TABLE}" \
    "${OVERALL}" \
    "${OUTPUT_MANIFEST}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# DISPLAY
###############################################################################

echo
echo "============================================================"
echo "OVERALL"
echo "============================================================"

column -t -s $'\t' \
    "${OVERALL}"

echo
echo "============================================================"
echo "STRICT VS SENSITIVE"
echo "============================================================"

column -t -s $'\t' \
    "${COMPARISON_TABLE}"

echo
echo "============================================================"
echo "REPEATED PAIRS"
echo "============================================================"

column -t -s $'\t' \
    "${DUPLICATE_TABLE}"

echo
echo "============================================================"
echo "CHECKPOINT"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X2C_COMPLETE.txt"

echo
echo "Step 36X2C completed successfully."
