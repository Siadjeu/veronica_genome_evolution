#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36X4A
#SBATCH --output=11_wgdi/logs/step36X4A_%j.out
#SBATCH --error=11_wgdi/logs/step36X4A_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

SOURCE_MANIFEST="${WGDI_ROOT}/00_admin/step36X3_ks_manifest.tsv"

KS_ROOT="${WGDI_ROOT}/06_ks"
CHUNK_ROOT="${KS_ROOT}/05_chunked_additional"

CHUNK_PAIR_DIR="${CHUNK_ROOT}/pairs"
CHUNK_CONFIG_DIR="${CHUNK_ROOT}/configs"
CHUNK_RESULT_DIR="${CHUNK_ROOT}/results"
CHUNK_WORK_DIR="${CHUNK_ROOT}/work"

QC_DIR="${WGDI_ROOT}/02_qc/additional_ks_chunks"
QC_RESULT_DIR="${QC_DIR}/results"

ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X4A"
LOG_DIR="${WGDI_ROOT}/logs"

CHUNK_MANIFEST="${ADMIN_DIR}/step36X4_chunk_manifest.tsv"
CHUNK_SUMMARY="${QC_DIR}/step36X4A_chunk_preparation_summary.tsv"

CHUNK_SIZE=5000
EXPECTED_TOTAL_PAIRS=104147
EXPECTED_COMPARISONS=3
EXPECTED_CHUNKS=23

mkdir -p \
    "${CHUNK_PAIR_DIR}" \
    "${CHUNK_CONFIG_DIR}" \
    "${CHUNK_RESULT_DIR}" \
    "${CHUNK_WORK_DIR}" \
    "${QC_RESULT_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

if [[ ! -s "${SOURCE_MANIFEST}" ]]; then
    echo "ERROR: Missing Step 36X3 Ks manifest:" >&2
    echo "${SOURCE_MANIFEST}" >&2
    exit 1
fi

rm -f \
    "${CHUNK_MANIFEST}" \
    "${CHUNK_SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36X4A_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${SOURCE_MANIFEST}" \
    "${CHUNK_PAIR_DIR}" \
    "${CHUNK_CONFIG_DIR}" \
    "${CHUNK_RESULT_DIR}" \
    "${QC_RESULT_DIR}" \
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
qc_result_dir = Path(sys.argv[6])
chunk_manifest = Path(sys.argv[7])
chunk_summary = Path(sys.argv[8])
chunk_size = int(sys.argv[9])


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

if len(source_rows) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 source comparisons; "
        f"found {len(source_rows)}."
    )

expected_comparisons = [
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
]

observed_comparisons = [
    row["comparison"]
    for row in source_rows
]

if observed_comparisons != expected_comparisons:
    raise SystemExit(
        "ERROR: Unexpected source comparison order: "
        f"{observed_comparisons}"
    )

chunk_rows = []
summary_rows = []
global_task_id = 0

for source_row in source_rows:

    if source_row["status"] != "READY":
        raise SystemExit(
            f"ERROR: Source dataset not READY: "
            f"{source_row['comparison']}"
        )

    if source_row["parameter_set"] != "strict":
        raise SystemExit(
            f"ERROR: Non-strict source dataset: "
            f"{source_row['comparison']}"
        )

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

    expected_pairs = int(
        source_row["pair_count"]
    )

    if len(pairs) != expected_pairs:
        raise SystemExit(
            f"ERROR: Pair-count mismatch for {comparison}: "
            f"manifest={expected_pairs}, file={len(pairs)}"
        )

    if len(pairs) != len(set(pairs)):
        raise SystemExit(
            f"ERROR: Duplicate source pairs in {comparison}."
        )

    number_of_chunks = math.ceil(
        len(pairs) / chunk_size
    )

    comparison_pair_total = 0

    for chunk_index in range(
        number_of_chunks
    ):
        start = (
            chunk_index
            * chunk_size
        )

        end = min(
            start + chunk_size,
            len(pairs),
        )

        chunk_pairs = pairs[
            start:end
        ]

        chunk_number = (
            chunk_index + 1
        )

        chunk_id = (
            f"{comparison}."
            f"chunk_{chunk_number:03d}"
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
            qc_result_dir
            / f"{chunk_id}.qc.tsv"
        )

        chunk_pair_file.write_text(
            "\n".join(chunk_pairs)
            + "\n",
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

        comparison_pair_total += len(
            chunk_pairs
        )

        chunk_rows.append(
            {
                "task_id": str(
                    global_task_id
                ),
                "comparison": comparison,
                "comparison_type": source_row[
                    "comparison_type"
                ],
                "species1": source_row[
                    "species1"
                ],
                "species2": source_row[
                    "species2"
                ],
                "chunk_number": str(
                    chunk_number
                ),
                "chunk_count_for_comparison": str(
                    number_of_chunks
                ),
                "start_pair_1based": str(
                    start + 1
                ),
                "end_pair_1based": str(
                    end
                ),
                "pair_count": str(
                    len(chunk_pairs)
                ),
                "pairs": relative(
                    chunk_pair_file
                ),
                "pep": relative(
                    pep_file
                ),
                "cds": relative(
                    cds_file
                ),
                "ks_config": relative(
                    chunk_config_file
                ),
                "ks_result": relative(
                    chunk_result_file
                ),
                "wgdi_stdout": relative(
                    chunk_stdout
                ),
                "wgdi_stderr": relative(
                    chunk_stderr
                ),
                "qc": relative(
                    chunk_qc
                ),
                "status": "READY",
            }
        )

        global_task_id += 1

    if comparison_pair_total != len(
        pairs
    ):
        raise SystemExit(
            f"ERROR: Chunk total mismatch "
            f"for {comparison}."
        )

    summary_rows.append(
        {
            "comparison": comparison,
            "source_pair_count": str(
                len(pairs)
            ),
            "chunk_size": str(
                chunk_size
            ),
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
        fieldnames=list(
            chunk_rows[0]
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        chunk_rows
    )

with chunk_summary.open(
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

print(
    f"Prepared {len(chunk_rows)} chunks "
    f"containing "
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

if [[ "${TASK_COUNT}" -ne "${EXPECTED_CHUNKS}" ]]; then
    echo "ERROR: Expected ${EXPECTED_CHUNKS} chunks; found ${TASK_COUNT}." >&2
    exit 1
fi

if [[ "${TOTAL_PAIRS}" -ne "${EXPECTED_TOTAL_PAIRS}" ]]; then
    echo "ERROR: Expected ${EXPECTED_TOTAL_PAIRS} pairs; found ${TOTAL_PAIRS}." >&2
    exit 1
fi

MAX_TASK_ID="$((TASK_COUNT - 1))"

cat > "${CHECKPOINT_DIR}/STEP36X4A_COMPLETE.txt" <<EOF2
checkpoint=step36X4A_prepare_chunked_Ks
date=$(date --iso-8601=seconds)
comparisons=${EXPECTED_COMPARISONS}
chunk_size=${CHUNK_SIZE}
chunk_tasks=${TASK_COUNT}
maximum_task_id=${MAX_TASK_ID}
total_pairs=${TOTAL_PAIRS}
alignment_software=mafft
ks_engine=YN00
primary_ks_column=ks_YN00
diagnostic_ks_column=ks_NG86
manifest=11_wgdi/00_admin/step36X4_chunk_manifest.tsv
summary=11_wgdi/02_qc/additional_ks_chunks/step36X4A_chunk_preparation_summary.tsv
status=PASS
next_step=step36X4B_run_chunked_Ks
EOF2

cp -f \
    "${CHUNK_MANIFEST}" \
    "${CHUNK_SUMMARY}" \
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
echo "Step 36X4A chunk preparation summary"
echo "============================================================"

column -t -s $'\t' \
    "${CHUNK_SUMMARY}"

echo
echo "============================================================"
echo "Checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X4A_COMPLETE.txt"
