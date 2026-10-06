#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36D2C
#SBATCH --output=11_wgdi/logs/step36D2C_%j.out
#SBATCH --error=11_wgdi/logs/step36D2C_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

CHUNK_MANIFEST="${WGDI_ROOT}/00_admin/step36D2_chunk_manifest.tsv"
SOURCE_MANIFEST="${WGDI_ROOT}/00_admin/step36D_ks_manifest.tsv"

MERGED_DIR="${WGDI_ROOT}/06_ks/06_merged"
QC_DIR="${WGDI_ROOT}/02_qc/ks_merged"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36D2C"
LOG_DIR="${WGDI_ROOT}/logs"

MERGED_MANIFEST="${ADMIN_DIR}/step36D2C_merged_ks_manifest.tsv"
SUMMARY="${QC_DIR}/step36D2C_merged_ks_summary.tsv"

mkdir -p \
    "${MERGED_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

for FILE in \
    "${CHUNK_MANIFEST}" \
    "${SOURCE_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required manifest missing: ${FILE}" >&2
        exit 1
    fi
done

rm -f \
    "${MERGED_MANIFEST}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36D2C_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${CHUNK_MANIFEST}" \
    "${SOURCE_MANIFEST}" \
    "${MERGED_DIR}" \
    "${MERGED_MANIFEST}" \
    "${SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path

project_root = Path(sys.argv[1])
chunk_manifest = Path(sys.argv[2])
source_manifest = Path(sys.argv[3])
merged_dir = Path(sys.argv[4])
merged_manifest = Path(sys.argv[5])
summary_path = Path(sys.argv[6])

expected_columns = [
    "id1",
    "id2",
    "ka_NG86",
    "ks_NG86",
    "ka_YN00",
    "ks_YN00",
]


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


def finite_float(value: str) -> float:
    parsed = float(value)

    if not math.isfinite(parsed):
        raise ValueError(value)

    return parsed


with chunk_manifest.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    chunk_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

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

chunks_by_comparison = defaultdict(list)

for row in chunk_rows:
    chunks_by_comparison[row["comparison"]].append(row)

summary_rows = []
manifest_rows = []

grand_expected = 0
grand_observed = 0

for source_row in source_rows:
    comparison = source_row["comparison"]
    expected_pair_count = int(source_row["pair_count"])
    source_pairs_path = project_root / source_row["pairs"]

    chunks = chunks_by_comparison.get(comparison, [])

    if not chunks:
        raise SystemExit(
            f"ERROR: No chunks found for {comparison}."
        )

    chunks.sort(
        key=lambda row: int(row["chunk_number"])
    )

    expected_chunk_numbers = list(
        range(1, len(chunks) + 1)
    )

    observed_chunk_numbers = [
        int(row["chunk_number"])
        for row in chunks
    ]

    if observed_chunk_numbers != expected_chunk_numbers:
        raise SystemExit(
            f"ERROR: Nonconsecutive chunk numbers for "
            f"{comparison}: {observed_chunk_numbers}"
        )

    with source_pairs_path.open(
        encoding="utf-8",
    ) as handle:
        source_pairs = [
            tuple(line.rstrip("\n").split("\t")[:2])
            for line in handle
            if line.strip()
        ]

    if len(source_pairs) != expected_pair_count:
        raise SystemExit(
            f"ERROR: Source pair count mismatch for "
            f"{comparison}."
        )

    all_rows = []
    result_pairs = []
    total_chunk_pairs = 0

    for chunk in chunks:
        result_path = project_root / chunk["ks_result"]
        qc_path = project_root / chunk["qc"]

        for required in [result_path, qc_path]:
            if (
                not required.is_file()
                or required.stat().st_size == 0
            ):
                raise SystemExit(
                    f"ERROR: Missing chunk output: {required}"
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
                f"ERROR: Invalid QC row count in {qc_path}."
            )

        qc = qc_rows[0]

        if qc.get("status") != "PASS":
            raise SystemExit(
                f"ERROR: Chunk QC is not PASS: {qc_path}"
            )

        if qc.get("completion_fraction") != "1.00000000":
            raise SystemExit(
                f"ERROR: Incomplete chunk: {qc_path}"
            )

        with result_path.open(
            newline="",
            encoding="utf-8-sig",
        ) as handle:
            reader = csv.DictReader(
                handle,
                delimiter="\t",
            )

            if reader.fieldnames != expected_columns:
                raise SystemExit(
                    f"ERROR: Unexpected columns in {result_path}: "
                    f"{reader.fieldnames}"
                )

            rows = list(reader)

        expected_chunk_pairs = int(chunk["pair_count"])

        if len(rows) != expected_chunk_pairs:
            raise SystemExit(
                f"ERROR: Chunk row count mismatch for "
                f"{result_path}: {len(rows)} versus "
                f"{expected_chunk_pairs}"
            )

        for row in rows:
            for field in expected_columns[2:]:
                try:
                    finite_float(row[field])
                except Exception:
                    raise SystemExit(
                        f"ERROR: Nonfinite or nonnumeric value "
                        f"in {result_path}: "
                        f"{row['id1']} {row['id2']} {field}="
                        f"{row[field]}"
                    )

        all_rows.extend(rows)
        result_pairs.extend(
            (row["id1"], row["id2"])
            for row in rows
        )
        total_chunk_pairs += len(rows)

    if total_chunk_pairs != expected_pair_count:
        raise SystemExit(
            f"ERROR: Merged pair count mismatch for "
            f"{comparison}: {total_chunk_pairs} versus "
            f"{expected_pair_count}"
        )

    if len(result_pairs) != len(set(result_pairs)):
        raise SystemExit(
            f"ERROR: Duplicate Ks result pairs for {comparison}."
        )

    if result_pairs != source_pairs:
        first_difference = None

        for index, (observed, expected) in enumerate(
            zip(result_pairs, source_pairs),
            start=1,
        ):
            if observed != expected:
                first_difference = (
                    index,
                    observed,
                    expected,
                )
                break

        if first_difference is None:
            first_difference = (
                min(len(result_pairs), len(source_pairs)) + 1,
                None,
                None,
            )

        raise SystemExit(
            f"ERROR: Merged pair order/content differs from "
            f"source pairs for {comparison}. First difference: "
            f"{first_difference}"
        )

    merged_path = (
        merged_dir
        / f"{comparison}.strict.ks.tsv"
    )

    temporary_path = Path(
        str(merged_path) + ".tmp"
    )

    with temporary_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=expected_columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(all_rows)

    temporary_path.replace(merged_path)

    ks_ng86_values = [
        float(row["ks_NG86"])
        for row in all_rows
    ]

    ks_yn00_values = [
        float(row["ks_YN00"])
        for row in all_rows
    ]

    ka_ng86_values = [
        float(row["ka_NG86"])
        for row in all_rows
    ]

    ka_yn00_values = [
        float(row["ka_YN00"])
        for row in all_rows
    ]

    summary_rows.append(
        {
            "comparison": comparison,
            "comparison_type": source_row[
                "comparison_type"
            ],
            "species1": source_row["species1"],
            "species2": source_row["species2"],
            "chunks": str(len(chunks)),
            "expected_pairs": str(expected_pair_count),
            "merged_rows": str(len(all_rows)),
            "unique_pairs": str(len(set(result_pairs))),
            "finite_ka_NG86": str(len(ka_ng86_values)),
            "finite_ks_NG86": str(len(ks_ng86_values)),
            "finite_ka_YN00": str(len(ka_yn00_values)),
            "finite_ks_YN00": str(len(ks_yn00_values)),
            "negative_ks_NG86": str(
                sum(value < -1e-12 for value in ks_ng86_values)
            ),
            "negative_ks_YN00": str(
                sum(value < -1e-12 for value in ks_yn00_values)
            ),
            "zero_ks_NG86": str(
                sum(abs(value) <= 1e-12 for value in ks_ng86_values)
            ),
            "zero_ks_YN00": str(
                sum(abs(value) <= 1e-12 for value in ks_yn00_values)
            ),
            "status": "PASS",
        }
    )

    manifest_rows.append(
        {
            "comparison": comparison,
            "comparison_type": source_row[
                "comparison_type"
            ],
            "species1": source_row["species1"],
            "species2": source_row["species2"],
            "parameter_set": "strict",
            "pair_count": str(expected_pair_count),
            "pairs": source_row["pairs"],
            "pep": source_row["pep"],
            "cds": source_row["cds"],
            "merged_ks": relative(merged_path),
            "merged_ks_sha256": sha256(merged_path),
            "status": "PASS",
        }
    )

    grand_expected += expected_pair_count
    grand_observed += len(all_rows)

if grand_expected != 339563:
    raise SystemExit(
        f"ERROR: Expected global pair count 339563; "
        f"found {grand_expected}."
    )

if grand_observed != grand_expected:
    raise SystemExit(
        f"ERROR: Global merged count mismatch: "
        f"{grand_observed} versus {grand_expected}."
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

with merged_manifest.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(manifest_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(manifest_rows)

print(
    f"Merged {grand_observed:,} Ks rows across "
    f"{len(manifest_rows)} comparisons."
)
PY

COMPARISON_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${MERGED_MANIFEST}"
)"

MERGED_PAIR_COUNT="$(
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
    ' "${MERGED_MANIFEST}"
)"

if [[ "${COMPARISON_COUNT}" -ne 12 ]]; then
    echo "ERROR: Only ${COMPARISON_COUNT}/12 comparisons merged." >&2
    exit 1
fi

if [[ "${MERGED_PAIR_COUNT}" -ne 339563 ]]; then
    echo "ERROR: Merged pair total is ${MERGED_PAIR_COUNT}, expected 339563." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36D2C_COMPLETE.txt" <<EOF2
checkpoint=step36D2C_merge_WGDI_Ks_chunks
date=$(date --iso-8601=seconds)
comparisons_expected=12
comparisons_merged=${COMPARISON_COUNT}
total_pairs_expected=339563
total_pairs_merged=${MERGED_PAIR_COUNT}
columns=id1,id2,ka_NG86,ks_NG86,ka_YN00,ks_YN00
primary_ks_column=ks_NG86
robustness_ks_column=ks_YN00
merged_directory=11_wgdi/06_ks/06_merged
manifest=11_wgdi/00_admin/step36D2C_merged_ks_manifest.tsv
summary=11_wgdi/02_qc/ks_merged/step36D2C_merged_ks_summary.tsv
status=PASS
next_step=step36E_prepare_and_run_WGDI_blockinfo
EOF2

cp -f \
    "${MERGED_MANIFEST}" \
    "${CHECKPOINT_DIR}/"

cp -f \
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
echo "Merged Ks summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Merged Ks manifest"
echo "============================================================"

column -t -s $'\t' \
    "${MERGED_MANIFEST}"

echo
echo "============================================================"
echo "Step 36D2C checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36D2C_COMPLETE.txt"
