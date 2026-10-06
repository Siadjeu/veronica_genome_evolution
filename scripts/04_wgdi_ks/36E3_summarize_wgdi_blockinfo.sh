#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36E3
#SBATCH --output=11_wgdi/logs/step36E3_%j.out
#SBATCH --error=11_wgdi/logs/step36E3_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

MANIFEST="${WGDI_ROOT}/00_admin/step36E_blockinfo_manifest.tsv"
QC_DIR="${WGDI_ROOT}/02_qc/blockinfo_results"
SUMMARY_DIR="${WGDI_ROOT}/02_qc/blockinfo_summary"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36E"
LOG_DIR="${WGDI_ROOT}/logs"

SUMMARY="${SUMMARY_DIR}/step36E_blockinfo_summary.tsv"
OVERALL="${SUMMARY_DIR}/step36E_blockinfo_overall.tsv"

mkdir -p \
    "${SUMMARY_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Missing Step 36E manifest." >&2
    exit 1
fi

rm -f \
    "${SUMMARY}" \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/STEP36E_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${MANIFEST}" \
    "${SUMMARY}" \
    "${OVERALL}" <<'PY'
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
summary_path = Path(sys.argv[3])
overall_path = Path(sys.argv[4])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        for block in iter(
            lambda: handle.read(1024 * 1024),
            b"",
        ):
            digest.update(block)

    return digest.hexdigest()


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

if len(manifest_rows) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 manifest rows; "
        f"found {len(manifest_rows)}."
    )

summary_rows = []

for row in manifest_rows:
    result = project_root / row["result"]
    qc = project_root / row["qc"]

    for required in [result, qc]:
        if not required.is_file() or required.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: Missing Step 36E output: {required}"
            )

    with qc.open(
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
            f"ERROR: Expected one QC row in {qc}."
        )

    q = qc_rows[0]

    if q.get("status") != "PASS":
        raise SystemExit(
            f"ERROR: QC status is not PASS for "
            f"{row['comparison']}."
        )

    if int(q["blockinfo_rows"]) <= 0:
        raise SystemExit(
            f"ERROR: No blockinfo rows for {row['comparison']}."
        )

    output = dict(q)
    output["result"] = row["result"]
    output["result_size_bytes"] = str(result.stat().st_size)
    output["result_sha256"] = sha256(result)

    summary_rows.append(output)

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

total_rows = sum(
    int(row["blockinfo_rows"])
    for row in summary_rows
)

total_rows_with_ks = sum(
    int(row["rows_with_any_ks_value"])
    for row in summary_rows
)

overall_rows = [
    ("comparisons_expected", "12"),
    ("comparisons_complete", str(len(summary_rows))),
    (
        "self_comparisons",
        str(
            sum(
                row["comparison_type"] == "self"
                for row in summary_rows
            )
        ),
    ),
    (
        "pairwise_comparisons",
        str(
            sum(
                row["comparison_type"] == "pairwise"
                for row in summary_rows
            )
        ),
    ),
    ("primary_ks_column", "ks_YN00"),
    ("total_blockinfo_rows", str(total_rows)),
    (
        "total_rows_with_any_ks_value",
        str(total_rows_with_ks),
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
    f"Validated 12 block-information outputs with "
    f"{total_rows:,} total rows."
)
PY

cat > "${CHECKPOINT_DIR}/STEP36E_COMPLETE.txt" <<EOF2
checkpoint=step36E_WGDI_blockinfo
date=$(date --iso-8601=seconds)
comparisons_expected=12
comparisons_complete=12
self_comparisons=9
pairwise_comparisons=3
parameter_set=strict
mg=25,25
primary_ks_column=ks_YN00
summary=11_wgdi/02_qc/blockinfo_summary/step36E_blockinfo_summary.tsv
overall=11_wgdi/02_qc/blockinfo_summary/step36E_blockinfo_overall.tsv
results_directory=11_wgdi/07_blockinfo/02_results
status=PASS
next_step=step36F_block_Ks_plots_and_peak_analysis
EOF2

cp -f \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/"

cp -f \
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
echo "Step 36E blockinfo summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36E overall"
echo "============================================================"

column -t -s $'\t' \
    "${OVERALL}"

echo
echo "============================================================"
echo "Step 36E checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36E_COMPLETE.txt"
