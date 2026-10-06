#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36F1
#SBATCH --output=11_wgdi/logs/step36F1_%j.out
#SBATCH --error=11_wgdi/logs/step36F1_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

BLOCKINFO_MANIFEST="${WGDI_ROOT}/00_admin/step36E_blockinfo_manifest.tsv"

OUT_ROOT="${WGDI_ROOT}/08_block_ks"
CURATED_DIR="${OUT_ROOT}/01_curated_blockinfo"
PAIR_KS_DIR="${OUT_ROOT}/02_pairwise_ks_tables"

QC_DIR="${WGDI_ROOT}/02_qc/block_ks_curated"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36F1"
LOG_DIR="${WGDI_ROOT}/logs"

CURATED_MANIFEST="${ADMIN_DIR}/step36F_curated_block_ks_manifest.tsv"
SUMMARY="${QC_DIR}/step36F1_block_filter_summary.tsv"
OVERALL="${QC_DIR}/step36F1_overall_summary.tsv"

mkdir -p \
    "${CURATED_DIR}" \
    "${PAIR_KS_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

if [[ ! -s "${BLOCKINFO_MANIFEST}" ]]; then
    echo "ERROR: Missing Step 36E manifest:" >&2
    echo "${BLOCKINFO_MANIFEST}" >&2
    exit 1
fi

rm -f \
    "${CURATED_MANIFEST}" \
    "${SUMMARY}" \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/STEP36F1_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${BLOCKINFO_MANIFEST}" \
    "${CURATED_DIR}" \
    "${PAIR_KS_DIR}" \
    "${CURATED_MANIFEST}" \
    "${SUMMARY}" \
    "${OVERALL}" <<'PY'
from __future__ import annotations

import csv
import math
import statistics
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
curated_dir = Path(sys.argv[3])
pair_ks_dir = Path(sys.argv[4])
output_manifest = Path(sys.argv[5])
summary_path = Path(sys.argv[6])
overall_path = Path(sys.argv[7])

PVALUE_MAX = 0.2
MIN_BLOCK_LENGTH = 10
KS_MIN_EXCLUSIVE = 0.0
KS_MAX_INCLUSIVE = 3.0
TANDEM_RATIO_MAX_SENSITIVITY = 0.5


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


def parse_float(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except Exception as error:
        raise SystemExit(
            f"ERROR: Could not parse {label}: {value!r}"
        ) from error

    if not math.isfinite(parsed):
        raise SystemExit(
            f"ERROR: Nonfinite {label}: {value!r}"
        )

    return parsed


def parse_underscore_floats(value: str) -> list[float]:
    values = []

    for token in value.strip().split("_"):
        if token == "":
            continue

        parsed = float(token)

        if not math.isfinite(parsed):
            continue

        # Normalize negative zero.
        if abs(parsed) <= 1e-12:
            parsed = 0.0

        values.append(parsed)

    return values


with manifest_path.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    manifest_rows = list(
        csv.DictReader(handle, delimiter="\t")
    )

if len(manifest_rows) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 blockinfo datasets; "
        f"found {len(manifest_rows)}."
    )

required_columns = {
    "id",
    "chr1",
    "chr2",
    "start1",
    "end1",
    "start2",
    "end2",
    "pvalue",
    "length",
    "ks_median",
    "ks_average",
    "homo1",
    "homo2",
    "homo3",
    "homo4",
    "homo5",
    "block1",
    "block2",
    "ks",
    "tandem_ratio",
    "density1",
    "density2",
    "class1",
    "class2",
}

summary_rows = []
output_manifest_rows = []

global_raw_blocks = 0
global_primary_blocks = 0
global_sensitivity_blocks = 0
global_positive_pair_values = 0

for manifest_row in manifest_rows:
    comparison = manifest_row["comparison"]
    comparison_type = manifest_row["comparison_type"]
    species1 = manifest_row["species1"]
    species2 = manifest_row["species2"]

    source_path = project_root / manifest_row["result"]

    if not source_path.is_file() or source_path.stat().st_size == 0:
        raise SystemExit(
            f"ERROR: Missing blockinfo file: {source_path}"
        )

    primary_path = (
        curated_dir
        / f"{comparison}.primary.filtered.blockinfo.csv"
    )

    sensitivity_path = (
        curated_dir
        / f"{comparison}.low_tandem.filtered.blockinfo.csv"
    )

    pair_ks_path = (
        pair_ks_dir
        / f"{comparison}.positive_pair_ks.tsv"
    )

    with source_path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []

        missing_columns = required_columns.difference(fieldnames)

        if missing_columns:
            raise SystemExit(
                f"ERROR: Missing columns in {source_path}: "
                + ",".join(sorted(missing_columns))
            )

        rows = list(reader)

    primary_rows = []
    sensitivity_rows = []
    pair_ks_rows = []

    fail_pvalue = 0
    fail_length = 0
    fail_nonpositive_ks = 0
    fail_high_ks = 0
    high_tandem_primary_blocks = 0
    block_ks_list_length_mismatches = 0

    for row in rows:
        pvalue = parse_float(row["pvalue"], "pvalue")
        block_length = int(float(row["length"]))
        ks_median = parse_float(row["ks_median"], "ks_median")
        ks_average = parse_float(row["ks_average"], "ks_average")
        tandem_ratio = parse_float(
            row["tandem_ratio"],
            "tandem_ratio",
        )

        ks_values = parse_underscore_floats(row["ks"])

        if len(ks_values) != block_length:
            block_ks_list_length_mismatches += 1

        for pair_index, ks_value in enumerate(
            ks_values,
            start=1,
        ):
            if 0.0 < ks_value <= KS_MAX_INCLUSIVE:
                pair_ks_rows.append(
                    {
                        "comparison": comparison,
                        "comparison_type": comparison_type,
                        "species1": species1,
                        "species2": species2,
                        "block_id": row["id"],
                        "pair_index": str(pair_index),
                        "block_length": str(block_length),
                        "block_pvalue": f"{pvalue:.10g}",
                        "block_ks_median": f"{ks_median:.10g}",
                        "block_ks_average": f"{ks_average:.10g}",
                        "tandem_ratio": f"{tandem_ratio:.10g}",
                        "pair_ks_YN00": f"{ks_value:.10g}",
                    }
                )

        if pvalue > PVALUE_MAX:
            fail_pvalue += 1
            continue

        if block_length < MIN_BLOCK_LENGTH:
            fail_length += 1
            continue

        if ks_median <= KS_MIN_EXCLUSIVE:
            fail_nonpositive_ks += 1
            continue

        if ks_median > KS_MAX_INCLUSIVE:
            fail_high_ks += 1
            continue

        primary_rows.append(row)

        if tandem_ratio > TANDEM_RATIO_MAX_SENSITIVITY:
            high_tandem_primary_blocks += 1
        else:
            sensitivity_rows.append(row)

    if block_ks_list_length_mismatches != 0:
        raise SystemExit(
            f"ERROR: {block_ks_list_length_mismatches} blocks in "
            f"{comparison} have a Ks-list count different from length."
        )

    for output_path, output_rows in [
        (primary_path, primary_rows),
        (sensitivity_path, sensitivity_rows),
    ]:
        temporary = Path(str(output_path) + ".tmp")

        with temporary.open(
            "w",
            newline="",
            encoding="utf-8",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(output_rows)

        temporary.replace(output_path)

    pair_fields = [
        "comparison",
        "comparison_type",
        "species1",
        "species2",
        "block_id",
        "pair_index",
        "block_length",
        "block_pvalue",
        "block_ks_median",
        "block_ks_average",
        "tandem_ratio",
        "pair_ks_YN00",
    ]

    with pair_ks_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=pair_fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(pair_ks_rows)

    primary_ks = [
        float(row["ks_median"])
        for row in primary_rows
    ]

    sensitivity_ks = [
        float(row["ks_median"])
        for row in sensitivity_rows
    ]

    summary_rows.append(
        {
            "comparison": comparison,
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "raw_blocks": str(len(rows)),
            "failed_pvalue": str(fail_pvalue),
            "failed_minimum_length": str(fail_length),
            "failed_nonpositive_median_ks": str(
                fail_nonpositive_ks
            ),
            "failed_median_ks_above_3": str(fail_high_ks),
            "primary_retained_blocks": str(len(primary_rows)),
            "primary_retained_fraction": (
                f"{len(primary_rows) / len(rows):.8f}"
                if rows else "0"
            ),
            "high_tandem_primary_blocks": str(
                high_tandem_primary_blocks
            ),
            "low_tandem_sensitivity_blocks": str(
                len(sensitivity_rows)
            ),
            "positive_pair_ks_values_0_to_3": str(
                len(pair_ks_rows)
            ),
            "primary_median_of_block_medians": (
                f"{statistics.median(primary_ks):.8f}"
                if primary_ks else "NA"
            ),
            "sensitivity_median_of_block_medians": (
                f"{statistics.median(sensitivity_ks):.8f}"
                if sensitivity_ks else "NA"
            ),
            "ks_list_length_mismatches": "0",
            "status": "PASS",
        }
    )

    output_manifest_rows.append(
        {
            "comparison": comparison,
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "source_blockinfo": relative(source_path),
            "primary_blockinfo": relative(primary_path),
            "low_tandem_blockinfo": relative(sensitivity_path),
            "positive_pair_ks": relative(pair_ks_path),
            "pvalue_max": str(PVALUE_MAX),
            "minimum_block_length": str(MIN_BLOCK_LENGTH),
            "median_ks_range": "(0,3]",
            "low_tandem_ratio_max": str(
                TANDEM_RATIO_MAX_SENSITIVITY
            ),
            "primary_block_count": str(len(primary_rows)),
            "low_tandem_block_count": str(
                len(sensitivity_rows)
            ),
            "status": "PASS",
        }
    )

    global_raw_blocks += len(rows)
    global_primary_blocks += len(primary_rows)
    global_sensitivity_blocks += len(sensitivity_rows)
    global_positive_pair_values += len(pair_ks_rows)

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

with output_manifest.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(output_manifest_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(output_manifest_rows)

overall_rows = [
    ("comparisons", "12"),
    ("raw_blocks", str(global_raw_blocks)),
    ("primary_filtered_blocks", str(global_primary_blocks)),
    (
        "low_tandem_sensitivity_blocks",
        str(global_sensitivity_blocks),
    ),
    (
        "positive_pair_ks_values_0_to_3",
        str(global_positive_pair_values),
    ),
    ("pvalue_max", str(PVALUE_MAX)),
    ("minimum_block_length", str(MIN_BLOCK_LENGTH)),
    ("median_ks_range", "(0,3]"),
    (
        "low_tandem_ratio_max",
        str(TANDEM_RATIO_MAX_SENSITIVITY),
    ),
    ("primary_ks_column", "ks_YN00"),
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
    f"Curated {global_raw_blocks:,} raw blocks into "
    f"{global_primary_blocks:,} primary filtered blocks."
)
PY

PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            n++
        }
        END {
            print n + 0
        }
    ' "${CURATED_MANIFEST}"
)"

if [[ "${PASS_COUNT}" -ne 12 ]]; then
    echo "ERROR: Only ${PASS_COUNT}/12 curated datasets passed." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36F1_COMPLETE.txt" <<EOF2
checkpoint=step36F1_curate_block_Ks
date=$(date --iso-8601=seconds)
comparisons=12
primary_ks_column=ks_YN00
pvalue_max=0.2
minimum_block_length=10
primary_median_ks_range=(0,3]
low_tandem_sensitivity_ratio_max=0.5
manifest=11_wgdi/00_admin/step36F_curated_block_ks_manifest.tsv
summary=11_wgdi/02_qc/block_ks_curated/step36F1_block_filter_summary.tsv
overall=11_wgdi/02_qc/block_ks_curated/step36F1_overall_summary.tsv
status=PASS
next_step=step36F2_prepare_WGDI_block_Ks_plots
EOF2

cp -f \
    "${CURATED_MANIFEST}" \
    "${SUMMARY}" \
    "${OVERALL}" \
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
echo "Step 36F1 filter summary"
echo "============================================================"

column -t -s $'\t' "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36F1 overall"
echo "============================================================"

column -t -s $'\t' "${OVERALL}"

echo
echo "============================================================"
echo "Step 36F1 checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36F1_COMPLETE.txt"
