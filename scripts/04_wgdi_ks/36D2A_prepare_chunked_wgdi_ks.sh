#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36D2A
#SBATCH --output=11_wgdi/logs/step36D2A_%j.out
#SBATCH --error=11_wgdi/logs/step36D2A_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"
SOURCE_MANIFEST="${WGDI_ROOT}/00_admin/step36D_ks_manifest.tsv"

KS_ROOT="${WGDI_ROOT}/06_ks"
CHUNK_PAIR_DIR="${KS_ROOT}/05_chunked/pairs"
CHUNK_CONFIG_DIR="${KS_ROOT}/05_chunked/configs"
CHUNK_RESULT_DIR="${KS_ROOT}/05_chunked/results"

QC_DIR="${WGDI_ROOT}/02_qc/ks_chunks"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36D2A"
LOG_DIR="${WGDI_ROOT}/logs"

CHUNK_MANIFEST="${ADMIN_DIR}/step36D2_chunk_manifest.tsv"
CHUNK_SUMMARY="${QC_DIR}/step36D2A_chunk_preparation_summary.tsv"

CHUNK_SIZE=5000

mkdir -p \
    "${CHUNK_PAIR_DIR}" \
    "${CHUNK_CONFIG_DIR}" \
    "${CHUNK_RESULT_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

if [[ ! -s "${SOURCE_MANIFEST}" ]]; then
    echo "ERROR: Missing Step 36D manifest:" >&2
    echo "${SOURCE_MANIFEST}" >&2
    exit 1
fi

rm -f \
    "${CHUNK_MANIFEST}" \
    "${CHUNK_SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36D2A_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${SOURCE_MANIFEST}" \
    "${CHUNK_PAIR_DIR}" \
    "${CHUNK_CONFIG_DIR}" \
    "${CHUNK_RESULT_DIR}" \
    "${CHUNK_MANIFEST}" \
    "${CHUNK_SUMMARY}" \
    "${CHUNK_SIZE}" <<'PY'
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
source_manifest = Path(sys.argv[2])
chunk_pair_dir = Path(sys.argv[3])
chunk_config_dir = Path(sys.argv[4])
chunk_result_dir = Path(sys.argv[5])
chunk_manifest = Path(sys.argv[6])
chunk_summary = Path(sys.argv[7])
chunk_size = int(sys.argv[8])


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


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

if len(source_rows) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 source comparisons; "
        f"found {len(source_rows)}."
    )

chunk_rows = []
summary_rows = []
global_task_id = 0

for source_row in source_rows:
    comparison = source_row["comparison"]

    source_pairs = (
        project_root
        / source_row["pairs"]
    )

    pep_file = (
        project_root
        / source_row["pep"]
    )

    cds_file = (
        project_root
        / source_row["cds"]
    )

    for required in [
        source_pairs,
        pep_file,
        cds_file,
    ]:
        if (
            not required.is_file()
            or required.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing input for {comparison}: "
                f"{required}"
            )

    with source_pairs.open(
        encoding="utf-8",
    ) as handle:
        pairs = [
            line.rstrip("\n")
            for line in handle
            if line.strip()
        ]

    expected_pairs = int(source_row["pair_count"])

    if len(pairs) != expected_pairs:
        raise SystemExit(
            f"ERROR: Pair-count mismatch for {comparison}: "
            f"manifest={expected_pairs}, file={len(pairs)}"
        )

    number_of_chunks = math.ceil(
        len(pairs) / chunk_size
    )

    comparison_pair_total = 0

    for chunk_index in range(number_of_chunks):
        start = chunk_index * chunk_size
        end = min(
            start + chunk_size,
            len(pairs),
        )

        chunk_pairs = pairs[start:end]
        chunk_number = chunk_index + 1
        chunk_id = (
            f"{comparison}.chunk_{chunk_number:03d}"
        )

        chunk_pair_file = (
            chunk_pair_dir
            / f"{chunk_id}.pairs.tsv"
        )

        chunk_config_file = (
            chunk_config_dir
            / f"{chunk_id}.ks.conf"
        )

        chunk_result_file = (
            chunk_result_dir
            / f"{chunk_id}.ks.tsv"
        )

        chunk_stdout = (
            chunk_result_dir
            / f"{chunk_id}.wgdi.stdout.txt"
        )

        chunk_stderr = (
            chunk_result_dir
            / f"{chunk_id}.wgdi.stderr.txt"
        )

        chunk_qc = (
            project_root
            / "11_wgdi/02_qc/ks_chunks/results"
            / f"{chunk_id}.qc.tsv"
        )

        chunk_pair_file.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        chunk_qc.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

        chunk_pair_file.write_text(
            "\n".join(chunk_pairs) + "\n",
            encoding="utf-8",
        )

        chunk_config_file.write_text(
            "\n".join(
                [
                    "[ks]",
                    f"cds_file = {cds_file}",
                    f"pep_file = {pep_file}",
                    "align_software = mafft",
                    f"pairs_file = {chunk_pair_file}",
                    f"ks_file = {chunk_result_file}",
                    "",
                ]
            ),
            encoding="utf-8",
        )

        comparison_pair_total += len(chunk_pairs)

        chunk_rows.append(
            {
                "task_id": str(global_task_id),
                "comparison": comparison,
                "comparison_type": source_row[
                    "comparison_type"
                ],
                "species1": source_row["species1"],
                "species2": source_row["species2"],
                "chunk_number": str(chunk_number),
                "chunk_count_for_comparison": str(
                    number_of_chunks
                ),
                "start_pair_1based": str(start + 1),
                "end_pair_1based": str(end),
                "pair_count": str(len(chunk_pairs)),
                "pairs": relative(chunk_pair_file),
                "pep": relative(pep_file),
                "cds": relative(cds_file),
                "ks_config": relative(
                    chunk_config_file
                ),
                "ks_result": relative(
                    chunk_result_file
                ),
                "wgdi_stdout": relative(chunk_stdout),
                "wgdi_stderr": relative(chunk_stderr),
                "qc": relative(chunk_qc),
                "status": "READY",
            }
        )

        global_task_id += 1

    if comparison_pair_total != len(pairs):
        raise SystemExit(
            f"ERROR: Chunk total mismatch for {comparison}."
        )

    summary_rows.append(
        {
            "comparison": comparison,
            "source_pair_count": str(len(pairs)),
            "chunk_size": str(chunk_size),
            "number_of_chunks": str(
                number_of_chunks
            ),
            "chunk_pair_total": str(
                comparison_pair_total
            ),
            "status": "PASS",
        }
    )

with chunk_manifest.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(chunk_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(chunk_rows)

with chunk_summary.open(
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

print(
    f"Prepared {len(chunk_rows)} chunks containing "
    f"{sum(int(row['pair_count']) for row in chunk_rows):,} pairs."
)
PY

TASK_COUNT="$(
    awk -F $'\t' '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
    ' "${CHUNK_MANIFEST}"
)"

TOTAL_PAIRS="$(
    awk -F $'\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) {
                if ($i == "pair_count") {
                    pair_col = i
                }
            }
            next
        }
        {
            sum += $pair_col
        }
        END {
            print sum + 0
        }
    ' "${CHUNK_MANIFEST}"
)"

if [[ "${TOTAL_PAIRS}" -ne 339563 ]]; then
    echo "ERROR: Expected 339563 total pairs; found ${TOTAL_PAIRS}." >&2
    exit 1
fi

MAX_TASK_ID="$((TASK_COUNT - 1))"

cat > "${CHECKPOINT_DIR}/STEP36D2A_COMPLETE.txt" <<EOF2
checkpoint=step36D2A_prepare_chunked_WGDI_Ks
date=$(date --iso-8601=seconds)
comparisons=12
chunk_size=${CHUNK_SIZE}
chunk_tasks=${TASK_COUNT}
maximum_task_id=${MAX_TASK_ID}
total_pairs=${TOTAL_PAIRS}
alignment_software=mafft
ks_engine=YN00
manifest=11_wgdi/00_admin/step36D2_chunk_manifest.tsv
summary=11_wgdi/02_qc/ks_chunks/step36D2A_chunk_preparation_summary.tsv
status=PASS
next_step=step36D2B_run_chunked_WGDI_Ks
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
echo "Chunk preparation summary"
echo "============================================================"

column -t -s $'\t' \
    "${CHUNK_SUMMARY}"

echo
echo "============================================================"
echo "Step 36D2A checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36D2A_COMPLETE.txt"
