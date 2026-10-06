#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36E1
#SBATCH --output=11_wgdi/logs/step36E1_%j.out
#SBATCH --error=11_wgdi/logs/step36E1_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

COLLINEARITY_MANIFEST="${WGDI_ROOT}/00_admin/wgdi_collinearity_manifest.tsv"
HOMOLOGY_MANIFEST="${WGDI_ROOT}/00_admin/wgdi_homology_manifest.tsv"
KS_MANIFEST="${WGDI_ROOT}/00_admin/step36D2C_merged_ks_manifest.tsv"

GFF_DIR="${WGDI_ROOT}/01_inputs/gff"
LENS_DIR="${WGDI_ROOT}/01_inputs/lens"

BLOCKINFO_ROOT="${WGDI_ROOT}/07_blockinfo"
CONFIG_DIR="${BLOCKINFO_ROOT}/01_configs"
RESULT_DIR="${BLOCKINFO_ROOT}/02_results"
WORK_DIR="${BLOCKINFO_ROOT}/03_work"

QC_DIR="${WGDI_ROOT}/02_qc/blockinfo_preparation"
ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36E1"
LOG_DIR="${WGDI_ROOT}/logs"

BLOCKINFO_MANIFEST="${ADMIN_DIR}/step36E_blockinfo_manifest.tsv"
SUMMARY="${QC_DIR}/step36E1_blockinfo_input_summary.tsv"

mkdir -p \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${WORK_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

for FILE in \
    "${COLLINEARITY_MANIFEST}" \
    "${HOMOLOGY_MANIFEST}" \
    "${KS_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required manifest missing: ${FILE}" >&2
        exit 1
    fi
done

rm -f \
    "${BLOCKINFO_MANIFEST}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36E1_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

python - \
    "${PROJECT_ROOT}" \
    "${COLLINEARITY_MANIFEST}" \
    "${HOMOLOGY_MANIFEST}" \
    "${KS_MANIFEST}" \
    "${GFF_DIR}" \
    "${LENS_DIR}" \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${WORK_DIR}" \
    "${BLOCKINFO_MANIFEST}" \
    "${SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
collinearity_manifest = Path(sys.argv[2])
homology_manifest = Path(sys.argv[3])
ks_manifest = Path(sys.argv[4])
gff_dir = Path(sys.argv[5])
lens_dir = Path(sys.argv[6])
config_dir = Path(sys.argv[7])
result_dir = Path(sys.argv[8])
work_dir = Path(sys.argv[9])
output_manifest = Path(sys.argv[10])
summary_path = Path(sys.argv[11])


def read_tsv(path: Path) -> list[dict[str, str]]:
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


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


def first_existing(candidates: list[Path], label: str) -> Path:
    for candidate in candidates:
        if candidate.is_file() and candidate.stat().st_size > 0:
            return candidate

    candidate_text = "\n".join(
        f"  {candidate}"
        for candidate in candidates
    )

    raise SystemExit(
        f"ERROR: Could not locate {label}. Tried:\n"
        f"{candidate_text}"
    )


def find_manifest_path(
    row: dict[str, str],
    candidate_columns: list[str],
    label: str,
) -> Path:
    for column in candidate_columns:
        value = row.get(column, "").strip()

        if not value:
            continue

        path = Path(value)

        if not path.is_absolute():
            path = project_root / path

        if path.is_file() and path.stat().st_size > 0:
            return path

    raise SystemExit(
        f"ERROR: Could not resolve {label} from manifest row. "
        f"Available columns: {sorted(row)}"
    )


collinearity_rows = read_tsv(collinearity_manifest)
homology_rows = read_tsv(homology_manifest)
ks_rows = read_tsv(ks_manifest)

strict_collinearity = {
    row["comparison"]: row
    for row in collinearity_rows
    if row.get("parameter_set") == "strict"
}

if len(strict_collinearity) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 strict collinearity datasets; "
        f"found {len(strict_collinearity)}."
    )

homology_by_comparison = {
    row["comparison"]: row
    for row in homology_rows
}

ks_by_comparison = {
    row["comparison"]: row
    for row in ks_rows
}

if len(ks_by_comparison) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 merged Ks datasets; "
        f"found {len(ks_by_comparison)}."
    )

manifest_rows = []
summary_rows = []

for task_id, comparison in enumerate(
    sorted(strict_collinearity)
):
    col_row = strict_collinearity[comparison]

    if comparison not in homology_by_comparison:
        raise SystemExit(
            f"ERROR: Homology manifest lacks {comparison}."
        )

    if comparison not in ks_by_comparison:
        raise SystemExit(
            f"ERROR: Merged Ks manifest lacks {comparison}."
        )

    hom_row = homology_by_comparison[comparison]
    ks_row = ks_by_comparison[comparison]

    species1 = col_row["species1"]
    species2 = col_row["species2"]
    comparison_type = col_row["comparison_type"]

    collinearity = find_manifest_path(
        col_row,
        [
            "collinearity",
            "result",
            "result_file",
            "output",
        ],
        f"strict collinearity file for {comparison}",
    )

    blast = find_manifest_path(
        hom_row,
        [
            "blast",
            "homology",
            "homology_file",
            "blast_file",
            "result",
            "result_file",
            "output",
        ],
        f"homology file for {comparison}",
    )

    merged_ks = find_manifest_path(
        ks_row,
        [
            "merged_ks",
            "ks",
            "ks_file",
            "result",
        ],
        f"merged Ks file for {comparison}",
    )

    gff1 = first_existing(
        [
            gff_dir / f"{species1}.wgdi.gff",
            gff_dir / f"{species1}.gff",
            gff_dir / f"{species1}.wgdi.gff.tsv",
        ],
        f"WGDI GFF for {species1}",
    )

    gff2 = first_existing(
        [
            gff_dir / f"{species2}.wgdi.gff",
            gff_dir / f"{species2}.gff",
            gff_dir / f"{species2}.wgdi.gff.tsv",
        ],
        f"WGDI GFF for {species2}",
    )

    lens1 = first_existing(
        [
            lens_dir / f"{species1}.wgdi.lens",
            lens_dir / f"{species1}.lens",
            lens_dir / f"{species1}.wgdi.lens.tsv",
        ],
        f"WGDI lens file for {species1}",
    )

    lens2 = first_existing(
        [
            lens_dir / f"{species2}.wgdi.lens",
            lens_dir / f"{species2}.lens",
            lens_dir / f"{species2}.wgdi.lens.tsv",
        ],
        f"WGDI lens file for {species2}",
    )

    config = (
        config_dir
        / f"{comparison}.blockinfo.conf"
    )

    result = (
        result_dir
        / f"{comparison}.blockinfo.csv"
    )

    qc = (
        project_root
        / "11_wgdi/02_qc/blockinfo_results"
        / f"{comparison}.blockinfo.qc.tsv"
    )

    stdout = (
        result_dir
        / f"{comparison}.blockinfo.wgdi.stdout.txt"
    )

    stderr = (
        result_dir
        / f"{comparison}.blockinfo.wgdi.stderr.txt"
    )

    comparison_work = (
        work_dir
        / comparison
    )

    config.write_text(
        "\n".join(
            [
                "[blockinfo]",
                f"blast = {blast}",
                f"gff1 = {gff1}",
                f"gff2 = {gff2}",
                f"lens1 = {lens1}",
                f"lens2 = {lens2}",
                f"collinearity = {collinearity}",
                "score = 100",
                "evalue = 1e-5",
                "repeat_number = 20",
                "position = order",
                f"ks = {merged_ks}",
                "ks_col = ks_YN00",
                f"savefile = {result}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    input_files = [
        blast,
        gff1,
        gff2,
        lens1,
        lens2,
        collinearity,
        merged_ks,
    ]

    for path in input_files:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: Missing or empty input for "
                f"{comparison}: {path}"
            )

    manifest_rows.append(
        {
            "task_id": str(task_id),
            "comparison": comparison,
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "blast": relative(blast),
            "gff1": relative(gff1),
            "gff2": relative(gff2),
            "lens1": relative(lens1),
            "lens2": relative(lens2),
            "collinearity": relative(collinearity),
            "ks": relative(merged_ks),
            "ks_col": "ks_YN00",
            "config": relative(config),
            "result": relative(result),
            "stdout": relative(stdout),
            "stderr": relative(stderr),
            "qc": relative(qc),
            "work_dir": relative(comparison_work),
            "status": "READY",
        }
    )

    summary_rows.append(
        {
            "task_id": str(task_id),
            "comparison": comparison,
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "blast_size_bytes": str(blast.stat().st_size),
            "collinearity_size_bytes": str(
                collinearity.stat().st_size
            ),
            "ks_size_bytes": str(merged_ks.stat().st_size),
            "gff1_size_bytes": str(gff1.stat().st_size),
            "gff2_size_bytes": str(gff2.stat().st_size),
            "lens1_size_bytes": str(lens1.stat().st_size),
            "lens2_size_bytes": str(lens2.stat().st_size),
            "ks_column": "ks_YN00",
            "status": "PASS",
        }
    )

with output_manifest.open(
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

print(
    f"Prepared {len(manifest_rows)} WGDI block-information "
    "configurations using ks_YN00."
)
PY

READY_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "READY" {
            count++
        }
        END {
            print count + 0
        }
    ' "${BLOCKINFO_MANIFEST}"
)"

if [[ "${READY_COUNT}" -ne 12 ]]; then
    echo "ERROR: Only ${READY_COUNT}/12 configurations are ready." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36E1_COMPLETE.txt" <<EOF2
checkpoint=step36E1_prepare_WGDI_blockinfo
date=$(date --iso-8601=seconds)
comparisons_expected=12
comparisons_ready=${READY_COUNT}
parameter_set=strict
mg=25,25
primary_ks_column=ks_YN00
score=100
evalue=1e-5
repeat_number=20
position=order
manifest=11_wgdi/00_admin/step36E_blockinfo_manifest.tsv
summary=11_wgdi/02_qc/blockinfo_preparation/step36E1_blockinfo_input_summary.tsv
status=PASS
next_step=step36E2_run_WGDI_blockinfo
EOF2

cp -f \
    "${BLOCKINFO_MANIFEST}" \
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
echo "Step 36E1 preparation summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36E blockinfo manifest"
echo "============================================================"

column -t -s $'\t' \
    "${BLOCKINFO_MANIFEST}"

echo
echo "============================================================"
echo "Step 36E1 checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36E1_COMPLETE.txt"
