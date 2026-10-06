#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36X5C
#SBATCH --output=11_wgdi/logs/step36X5C_%j.out
#SBATCH --error=11_wgdi/logs/step36X5C_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

MANIFEST="${WGDI_ROOT}/00_admin/step36X5_blockinfo_manifest.tsv"

COLLINEARITY_SUMMARY="${WGDI_ROOT}/02_qc/additional_collinearity_summary/step36X2_strict_sensitive_comparison.tsv"

SUMMARY_DIR="${WGDI_ROOT}/02_qc/additional_blockinfo_summary"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X5"
LOG_DIR="${WGDI_ROOT}/logs"

SUMMARY="${SUMMARY_DIR}/step36X5_blockinfo_summary.tsv"
OVERALL="${SUMMARY_DIR}/step36X5_blockinfo_overall.tsv"

mkdir -p \
    "${SUMMARY_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

for FILE in \
    "${MANIFEST}" \
    "${COLLINEARITY_SUMMARY}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required file missing: ${FILE}" >&2
        exit 1
    fi
done

rm -f \
    "${SUMMARY}" \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/STEP36X5_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${MANIFEST}" \
    "${COLLINEARITY_SUMMARY}" \
    "${SUMMARY}" \
    "${OVERALL}" <<'PY'
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
collinearity_summary_path = Path(
    sys.argv[3]
)
summary_path = Path(sys.argv[4])
overall_path = Path(sys.argv[5])


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


def read_tsv(path):
    with path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


manifest_rows = read_tsv(
    manifest_path
)

collinearity_rows = read_tsv(
    collinearity_summary_path
)

if len(manifest_rows) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 blockinfo manifest rows; "
        f"found {len(manifest_rows)}."
    )

expected_comparisons = {
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
}

observed = {
    row["comparison"]
    for row in manifest_rows
}

if observed != expected_comparisons:
    raise SystemExit(
        f"ERROR: Unexpected manifest comparison set: "
        f"{sorted(observed)}"
    )

strict_blocks = {}

for row in collinearity_rows:
    comparison = row["comparison"]

    if comparison in expected_comparisons:
        strict_blocks[comparison] = int(
            row["strict_blocks"]
        )

if set(strict_blocks) != expected_comparisons:
    raise SystemExit(
        "ERROR: Could not recover strict block counts "
        "for all 3 comparisons."
    )

summary_rows = []

for row in manifest_rows:

    comparison = row["comparison"]

    result = (
        project_root
        / row["result"]
    )

    qc = (
        project_root
        / row["qc"]
    )

    for required in [
        result,
        qc,
    ]:
        if (
            not required.is_file()
            or required.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing Step 36X5 output: "
                f"{required}"
            )

    qc_rows = read_tsv(
        qc
    )

    if len(qc_rows) != 1:
        raise SystemExit(
            f"ERROR: Expected one QC row "
            f"in {qc}."
        )

    q = qc_rows[0]

    if q.get("status") != "PASS":
        raise SystemExit(
            f"ERROR: QC is not PASS "
            f"for {comparison}."
        )

    blockinfo_rows = int(
        q["blockinfo_rows"]
    )

    rows_with_ks = int(
        q["rows_with_any_ks_value"]
    )

    expected_blocks = strict_blocks[
        comparison
    ]

    if blockinfo_rows != expected_blocks:
        raise SystemExit(
            f"ERROR: Block count mismatch "
            f"for {comparison}: "
            f"strict_collinearity={expected_blocks}, "
            f"blockinfo={blockinfo_rows}"
        )

    if rows_with_ks != blockinfo_rows:
        raise SystemExit(
            f"ERROR: Not every blockinfo row "
            f"has Ks for {comparison}: "
            f"{rows_with_ks}/{blockinfo_rows}"
        )

    if (
        int(
            q[
                "nonfinite_scalar_ks_values"
            ]
        )
        != 0
    ):
        raise SystemExit(
            f"ERROR: Nonfinite scalar Ks "
            f"values detected for {comparison}."
        )

    output = dict(q)

    output[
        "strict_collinearity_blocks"
    ] = str(expected_blocks)

    output[
        "blockinfo_equals_strict_blocks"
    ] = "true"

    output[
        "all_blocks_have_ks"
    ] = "true"

    output["result"] = row[
        "result"
    ]

    output[
        "result_size_bytes"
    ] = str(
        result.stat().st_size
    )

    output[
        "result_sha256"
    ] = sha256(
        result
    )

    summary_rows.append(
        output
    )

summary_rows.sort(
    key=lambda x: x["comparison"]
)

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

total_rows = sum(
    int(row["blockinfo_rows"])
    for row in summary_rows
)

total_rows_with_ks = sum(
    int(
        row["rows_with_any_ks_value"]
    )
    for row in summary_rows
)

total_strict_blocks = sum(
    int(
        row[
            "strict_collinearity_blocks"
        ]
    )
    for row in summary_rows
)

if total_rows != 4080:
    raise SystemExit(
        f"ERROR: Expected 4080 total blockinfo rows; "
        f"found {total_rows}."
    )

if total_strict_blocks != 4080:
    raise SystemExit(
        f"ERROR: Expected 4080 strict blocks; "
        f"found {total_strict_blocks}."
    )

if total_rows_with_ks != 4080:
    raise SystemExit(
        f"ERROR: Expected Ks for all 4080 blocks; "
        f"found {total_rows_with_ks}."
    )

overall_rows = [
    (
        "comparisons_expected",
        "3",
    ),
    (
        "comparisons_complete",
        str(len(summary_rows)),
    ),
    (
        "self_comparisons",
        "0",
    ),
    (
        "pairwise_comparisons",
        "3",
    ),
    (
        "primary_ks_column",
        "ks_YN00",
    ),
    (
        "strict_collinearity_blocks",
        str(total_strict_blocks),
    ),
    (
        "total_blockinfo_rows",
        str(total_rows),
    ),
    (
        "total_rows_with_any_ks_value",
        str(total_rows_with_ks),
    ),
    (
        "blockinfo_equals_strict_blocks",
        "true",
    ),
    (
        "all_blocks_have_ks",
        "true",
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
    f"Validated 3 blockinfo outputs with "
    f"{total_rows:,} total rows; "
    "all strict blocks represented and "
    "all blocks contain Ks."
)
PY

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36X5_COMPLETE.txt" <<EOF2
checkpoint=step36X5_additional_WGDI_blockinfo
date=$(date --iso-8601=seconds)
comparisons_expected=3
comparisons_complete=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
self_comparisons=0
pairwise_comparisons=3
parameter_set=strict
mg=25,25
primary_ks_column=ks_YN00
strict_blocks_expected=4080
blockinfo_rows_expected=4080
all_strict_blocks_represented=true
all_blocks_have_ks=true
summary=11_wgdi/02_qc/additional_blockinfo_summary/step36X5_blockinfo_summary.tsv
overall=11_wgdi/02_qc/additional_blockinfo_summary/step36X5_blockinfo_overall.tsv
results_directory=11_wgdi/07_blockinfo/02_results
status=PASS
next_step=step36X6_additional_block_Ks_filtering
EOF2

cp -f \
    "${SUMMARY}" \
    "${OVERALL}" \
    "${MANIFEST}" \
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
echo "Step 36X5 blockinfo summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36X5 overall"
echo "============================================================"

column -t -s $'\t' \
    "${OVERALL}"

echo
echo "============================================================"
echo "Step 36X5 checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X5_COMPLETE.txt"
