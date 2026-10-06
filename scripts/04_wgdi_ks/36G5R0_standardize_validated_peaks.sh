#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=wgdi36G5R0
#SBATCH --output=11_wgdi/logs/step36G5R0_%j.out
#SBATCH --error=11_wgdi/logs/step36G5R0_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

INPUT="11_wgdi/08_block_ks/08_corrected_peak_validation/tables/step36F4R_validated_peaks.tsv"

OUT_DIR="11_wgdi/11_ksfigure/00_standardized_inputs"
QC_DIR="11_wgdi/02_qc/ksfigure_multipeak"
CHECKPOINT_DIR="11_wgdi/checkpoints/step36G5R0"

OUTPUT="${OUT_DIR}/step36F4R_validated_peaks.standardized.tsv"
COLUMN_REPORT="${QC_DIR}/step36G5R0_peak_column_detection.tsv"
SUMMARY="${QC_DIR}/step36G5R0_standardized_peak_summary.tsv"

mkdir -p \
    "${OUT_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "11_wgdi/logs"

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: Missing validated peak table: ${INPUT}" >&2
    exit 1
fi

rm -f \
    "${OUTPUT}" \
    "${COLUMN_REPORT}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36G5R0_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${INPUT}" \
    "${OUTPUT}" \
    "${COLUMN_REPORT}" \
    "${SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import math
import re
import sys
from collections import Counter
from pathlib import Path

input_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
column_report_path = Path(sys.argv[3])
summary_path = Path(sys.argv[4])


def normalize_header(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


with input_path.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    reader = csv.DictReader(
        handle,
        delimiter="\t",
    )
    rows = list(reader)
    original_columns = reader.fieldnames or []

if not original_columns:
    raise SystemExit(
        "ERROR: Validated peak table has no header."
    )

if not rows:
    raise SystemExit(
        "ERROR: Validated peak table has no data rows."
    )

normalized_to_original: dict[str, str] = {}

for column in original_columns:
    normalized = normalize_header(column)

    if normalized in normalized_to_original:
        raise SystemExit(
            "ERROR: Header normalization created a duplicate: "
            f"{column!r} and "
            f"{normalized_to_original[normalized]!r}"
        )

    normalized_to_original[normalized] = column

print("Original columns:")
for index, column in enumerate(original_columns, start=1):
    print(
        f"  {index:02d}. {column!r} "
        f"-> {normalize_header(column)!r}"
    )


def first_alias(
    aliases: list[str],
) -> str | None:
    for alias in aliases:
        normalized_alias = normalize_header(alias)

        if normalized_alias in normalized_to_original:
            return normalized_to_original[normalized_alias]

    return None


comparison_column = first_alias(
    [
        "comparison",
        "comparison_id",
        "comparison_name",
        "dataset",
        "contrast",
        "analysis",
    ]
)

if comparison_column is None:
    comparison_candidates = [
        original
        for normalized, original
        in normalized_to_original.items()
        if (
            "comparison" in normalized
            or "contrast" in normalized
        )
    ]

    if len(comparison_candidates) == 1:
        comparison_column = comparison_candidates[0]

if comparison_column is None:
    raise SystemExit(
        "ERROR: Could not identify the comparison column. "
        f"Available columns: {original_columns}"
    )

###############################################################################
# DETECT PEAK POSITION COLUMN
###############################################################################

peak_aliases = [
    "peak_ks",
    "component_peak_ks",
    "validated_peak_ks",
    "corrected_peak_ks",
    "final_peak_ks",
    "peak_position",
    "peak_position_ks",
    "peak_center",
    "peak_centre",
    "peak_location",
    "mode_ks",
    "kde_peak",
    "kde_peak_ks",
    "peak",
    "ks_peak",
    "ks_mode",
    "mode",
]

peak_column = first_alias(peak_aliases)

candidate_scores: list[
    tuple[int, int, str]
] = []

if peak_column is None:
    for column in original_columns:
        normalized = normalize_header(column)

        if column == comparison_column:
            continue

        numeric_values = []

        for row in rows:
            value = str(row.get(column, "")).strip()

            if value == "":
                continue

            try:
                number = float(value)
            except ValueError:
                continue

            if math.isfinite(number):
                numeric_values.append(number)

        if not numeric_values:
            continue

        score = 0

        if "peak" in normalized:
            score += 100

        if "ks" in normalized:
            score += 60

        if "centre" in normalized or "center" in normalized:
            score += 50

        if "mode" in normalized:
            score += 35

        if "position" in normalized:
            score += 30

        if "location" in normalized:
            score += 20

        # Peak Ks values should normally be within the analysed interval.
        within_ks_range = sum(
            0 <= value <= 3
            for value in numeric_values
        )

        if within_ks_range == len(numeric_values):
            score += 25

        candidate_scores.append(
            (
                score,
                len(numeric_values),
                column,
            )
        )

    candidate_scores.sort(
        key=lambda item: (
            item[0],
            item[1],
            item[2],
        ),
        reverse=True,
    )

    if candidate_scores:
        best_score, _, best_column = candidate_scores[0]

        if best_score >= 60:
            peak_column = best_column

if peak_column is None:
    print(
        "Numeric candidate columns and scores:",
        file=sys.stderr,
    )

    for score, count, column in candidate_scores:
        print(
            f"  score={score:3d} "
            f"numeric_rows={count:3d} "
            f"column={column!r}",
            file=sys.stderr,
        )

    raise SystemExit(
        "ERROR: Could not identify a Ks peak-position column. "
        f"Available columns: {original_columns}"
    )

###############################################################################
# OPTIONAL COLUMNS
###############################################################################

validation_column = first_alias(
    [
        "validation_class",
        "peak_validation_class",
        "validation_status",
        "robustness_class",
        "support_class",
        "class",
        "status",
    ]
)

component_column = first_alias(
    [
        "component_id",
        "peak_id",
        "component",
        "peak_number",
        "component_number",
    ]
)

lower_column = first_alias(
    [
        "ci_lower",
        "bootstrap_ci_lower",
        "peak_ci_lower",
        "lower_ci",
        "window_lower",
        "lower",
    ]
)

upper_column = first_alias(
    [
        "ci_upper",
        "bootstrap_ci_upper",
        "peak_ci_upper",
        "upper_ci",
        "window_upper",
        "upper",
    ]
)

###############################################################################
# STANDARDIZE AND VALIDATE
###############################################################################

standardized_rows = []

for row_number, row in enumerate(rows, start=2):
    comparison = str(
        row.get(comparison_column, "")
    ).strip()

    raw_peak = str(
        row.get(peak_column, "")
    ).strip()

    if not comparison:
        raise SystemExit(
            f"ERROR: Empty comparison at input line {row_number}."
        )

    try:
        peak_ks = float(raw_peak)
    except ValueError as error:
        raise SystemExit(
            f"ERROR: Non-numeric peak value {raw_peak!r} "
            f"at input line {row_number} in column "
            f"{peak_column!r}."
        ) from error

    if not math.isfinite(peak_ks):
        raise SystemExit(
            f"ERROR: Non-finite peak value at line {row_number}."
        )

    if peak_ks <= 0 or peak_ks > 3:
        raise SystemExit(
            f"ERROR: Peak Ks outside (0,3] at line "
            f"{row_number}: {peak_ks}"
        )

    validation_class = (
        str(row.get(validation_column, "")).strip()
        if validation_column
        else "VALIDATED"
    )

    if not validation_class:
        validation_class = "VALIDATED"

    component_id = (
        str(row.get(component_column, "")).strip()
        if component_column
        else ""
    )

    ci_lower = (
        str(row.get(lower_column, "")).strip()
        if lower_column
        else ""
    )

    ci_upper = (
        str(row.get(upper_column, "")).strip()
        if upper_column
        else ""
    )

    standardized_rows.append(
        {
            "comparison": comparison,
            "component_id": component_id,
            "peak_ks": f"{peak_ks:.12g}",
            "validation_class": validation_class,
            "ci_lower": ci_lower,
            "ci_upper": ci_upper,
            "source_peak_column": peak_column,
            "source_row_number": str(row_number),
            "status": "PASS",
        }
    )

# Ensure comparison/peak combinations are unique.
seen = set()

for row in standardized_rows:
    key = (
        row["comparison"],
        row["peak_ks"],
    )

    if key in seen:
        raise SystemExit(
            "ERROR: Duplicate comparison/peak combination: "
            f"{key}"
        )

    seen.add(key)

standardized_rows.sort(
    key=lambda row: (
        row["comparison"],
        float(row["peak_ks"]),
    )
)

with output_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(standardized_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(standardized_rows)

column_report_rows = [
    {
        "role": "comparison",
        "detected_column": comparison_column,
        "required": "yes",
        "status": "PASS",
    },
    {
        "role": "peak_ks",
        "detected_column": peak_column,
        "required": "yes",
        "status": "PASS",
    },
    {
        "role": "validation_class",
        "detected_column": validation_column or "not_present",
        "required": "no",
        "status": "PASS",
    },
    {
        "role": "component_id",
        "detected_column": component_column or "not_present",
        "required": "no",
        "status": "PASS",
    },
    {
        "role": "ci_lower",
        "detected_column": lower_column or "not_present",
        "required": "no",
        "status": "PASS",
    },
    {
        "role": "ci_upper",
        "detected_column": upper_column or "not_present",
        "required": "no",
        "status": "PASS",
    },
]

with column_report_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(column_report_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(column_report_rows)

counts = Counter(
    row["comparison"]
    for row in standardized_rows
)

summary_rows = []

for comparison in sorted(counts):
    comparison_rows = [
        row
        for row in standardized_rows
        if row["comparison"] == comparison
    ]

    summary_rows.append(
        {
            "comparison": comparison,
            "validated_peak_count": str(
                len(comparison_rows)
            ),
            "validated_peaks": ",".join(
                row["peak_ks"]
                for row in comparison_rows
            ),
            "validation_classes": ",".join(
                row["validation_class"]
                for row in comparison_rows
            ),
            "status": "PASS",
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

print()
print(f"Detected comparison column: {comparison_column}")
print(f"Detected peak Ks column:    {peak_column}")
print(f"Standardized peak rows:     {len(standardized_rows)}")
print(f"Comparisons represented:    {len(summary_rows)}")
PY

PEAK_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${OUTPUT}"
)"

COMPARISON_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${SUMMARY}"
)"

if [[ "${PEAK_COUNT}" -lt 12 ]]; then
    echo "ERROR: Only ${PEAK_COUNT} validated peaks were standardized." >&2
    exit 1
fi

if [[ "${COMPARISON_COUNT}" -ne 12 ]]; then
    echo "ERROR: Expected 12 comparisons; found ${COMPARISON_COUNT}." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36G5R0_COMPLETE.txt" <<EOF2
checkpoint=step36G5R0_standardize_validated_peaks
date=$(date --iso-8601=seconds)
input=${INPUT}
output=${OUTPUT}
validated_peak_rows=${PEAK_COUNT}
comparisons=${COMPARISON_COUNT}
peak_column_detection=automatic_with_alias_and_numeric_validation
ks_range=(0,3]
status=PASS
next_step=step36G5R1_prepare_multipeak_WGDI_KsFigure
EOF2

cp -f \
    "${OUTPUT}" \
    "${COLUMN_REPORT}" \
    "${SUMMARY}" \
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

echo
echo "============================================================"
echo "Detected columns"
echo "============================================================"

column -t -s $'\t' \
    "${COLUMN_REPORT}"

echo
echo "============================================================"
echo "Standardized validated peaks"
echo "============================================================"

column -t -s $'\t' \
    "${OUTPUT}"

echo
echo "============================================================"
echo "Peak summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36G5R0_COMPLETE.txt"
