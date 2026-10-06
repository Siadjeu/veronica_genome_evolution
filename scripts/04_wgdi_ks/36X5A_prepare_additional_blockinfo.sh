#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36X5A
#SBATCH --output=11_wgdi/logs/step36X5A_%j.out
#SBATCH --error=11_wgdi/logs/step36X5A_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

COLLINEARITY_MANIFEST="${WGDI_ROOT}/00_admin/step36X2_wgdi_collinearity_manifest.tsv"
KS_MANIFEST="${WGDI_ROOT}/00_admin/step36X4C_merged_ks_manifest.tsv"

GFF_DIR="${WGDI_ROOT}/01_inputs/gff"
LENS_DIR="${WGDI_ROOT}/01_inputs/lens"
HOMOLOGY_DIR="${WGDI_ROOT}/04_homology/results"

BLOCKINFO_ROOT="${WGDI_ROOT}/07_blockinfo"
CONFIG_DIR="${BLOCKINFO_ROOT}/01_configs"
RESULT_DIR="${BLOCKINFO_ROOT}/02_results"
WORK_DIR="${BLOCKINFO_ROOT}/03_work"

QC_DIR="${WGDI_ROOT}/02_qc/additional_blockinfo_preparation"
QC_RESULT_DIR="${WGDI_ROOT}/02_qc/additional_blockinfo_results"

ADMIN_DIR="${WGDI_ROOT}/00_admin"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X5A"
LOG_DIR="${WGDI_ROOT}/logs"

BLOCKINFO_MANIFEST="${ADMIN_DIR}/step36X5_blockinfo_manifest.tsv"
SUMMARY="${QC_DIR}/step36X5A_blockinfo_input_summary.tsv"

mkdir -p \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${WORK_DIR}" \
    "${QC_DIR}" \
    "${QC_RESULT_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

###############################################################################
# REQUIRE PREVIOUS VALIDATED STEPS
###############################################################################

for FILE in \
    "${COLLINEARITY_MANIFEST}" \
    "${KS_MANIFEST}" \
    "${WGDI_ROOT}/checkpoints/step36X2C/STEP36X2C_COMPLETE.txt" \
    "${WGDI_ROOT}/checkpoints/step36X4C/STEP36X4C_COMPLETE.txt"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required previous-step file missing: ${FILE}" >&2
        exit 1
    fi
done

for CHECKPOINT in \
    "${WGDI_ROOT}/checkpoints/step36X2C/STEP36X2C_COMPLETE.txt" \
    "${WGDI_ROOT}/checkpoints/step36X4C/STEP36X4C_COMPLETE.txt"
do
    if ! grep -q '^status=PASS$' "${CHECKPOINT}"; then
        echo "ERROR: Previous checkpoint is not PASS: ${CHECKPOINT}" >&2
        exit 1
    fi
done

###############################################################################
# TARGETED CLEANUP
###############################################################################

rm -f \
    "${BLOCKINFO_MANIFEST}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36X5A_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# PREPARE CONFIGURATIONS
###############################################################################

python - \
    "${PROJECT_ROOT}" \
    "${COLLINEARITY_MANIFEST}" \
    "${KS_MANIFEST}" \
    "${GFF_DIR}" \
    "${LENS_DIR}" \
    "${HOMOLOGY_DIR}" \
    "${CONFIG_DIR}" \
    "${RESULT_DIR}" \
    "${WORK_DIR}" \
    "${QC_RESULT_DIR}" \
    "${BLOCKINFO_MANIFEST}" \
    "${SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
collinearity_manifest = Path(sys.argv[2])
ks_manifest = Path(sys.argv[3])
gff_dir = Path(sys.argv[4])
lens_dir = Path(sys.argv[5])
homology_dir = Path(sys.argv[6])
config_dir = Path(sys.argv[7])
result_dir = Path(sys.argv[8])
work_dir = Path(sys.argv[9])
qc_result_dir = Path(sys.argv[10])
output_manifest = Path(sys.argv[11])
summary_path = Path(sys.argv[12])


def read_tsv(path: Path):
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


def resolve_manifest_path(row, columns, label):
    for column in columns:
        value = row.get(column, "").strip()

        if not value:
            continue

        path = Path(value)

        if not path.is_absolute():
            path = project_root / path

        if path.is_file() and path.stat().st_size > 0:
            return path

    raise SystemExit(
        f"ERROR: Could not resolve {label}. "
        f"Available columns: {sorted(row)}"
    )


def first_existing(candidates, label):
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            return path

    raise SystemExit(
        f"ERROR: Could not locate {label}. Tried:\n"
        + "\n".join(str(x) for x in candidates)
    )


collinearity_rows = read_tsv(
    collinearity_manifest
)

ks_rows = read_tsv(
    ks_manifest
)

strict_rows = [
    row
    for row in collinearity_rows
    if row.get("parameter_set") == "strict"
]

if len(strict_rows) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 strict collinearity datasets; "
        f"found {len(strict_rows)}."
    )

expected = {
    "VSCU_VSER": ("VSCU", "VSER"),
    "VSCU_VPAN": ("VSCU", "VPAN"),
    "VPAN_VPER": ("VPAN", "VPER"),
}

strict_by_comparison = {
    row["comparison"]: row
    for row in strict_rows
}

if set(strict_by_comparison) != set(expected):
    raise SystemExit(
        "ERROR: Unexpected strict comparison set: "
        f"{sorted(strict_by_comparison)}"
    )

ks_by_comparison = {
    row["comparison"]: row
    for row in ks_rows
}

if set(ks_by_comparison) != set(expected):
    raise SystemExit(
        "ERROR: Unexpected merged-Ks comparison set: "
        f"{sorted(ks_by_comparison)}"
    )

manifest_rows = []
summary_rows = []

for task_id, comparison in enumerate(
    [
        "VSCU_VSER",
        "VSCU_VPAN",
        "VPAN_VPER",
    ]
):
    col_row = strict_by_comparison[
        comparison
    ]

    ks_row = ks_by_comparison[
        comparison
    ]

    species1 = col_row["species1"]
    species2 = col_row["species2"]
    comparison_type = col_row[
        "comparison_type"
    ]

    if (
        species1,
        species2,
    ) != expected[comparison]:
        raise SystemExit(
            f"ERROR: Species orientation mismatch "
            f"for {comparison}: "
            f"{species1},{species2}"
        )

    if comparison_type != "pairwise":
        raise SystemExit(
            f"ERROR: Expected pairwise comparison "
            f"for {comparison}."
        )

    if ks_row.get("status") != "PASS":
        raise SystemExit(
            f"ERROR: Merged Ks is not PASS "
            f"for {comparison}."
        )

    if ks_row.get(
        "parameter_set"
    ) != "strict":
        raise SystemExit(
            f"ERROR: Merged Ks is not strict "
            f"for {comparison}."
        )

    collinearity = resolve_manifest_path(
        col_row,
        [
            "collinearity",
            "result",
            "result_file",
            "output",
        ],
        f"strict collinearity for {comparison}",
    )

    merged_ks = resolve_manifest_path(
        ks_row,
        [
            "merged_ks",
            "ks",
            "ks_file",
            "result",
        ],
        f"merged Ks for {comparison}",
    )

    blast = first_existing(
        [
            homology_dir
            / f"{comparison}.blast.tsv.gz",
            homology_dir
            / f"{comparison}.blast.tsv",
        ],
        f"DIAMOND homology for {comparison}",
    )

    gff1 = first_existing(
        [
            gff_dir
            / f"{species1}.wgdi.gff",
            gff_dir
            / f"{species1}.gff",
            gff_dir
            / f"{species1}.wgdi.gff.tsv",
        ],
        f"WGDI GFF for {species1}",
    )

    gff2 = first_existing(
        [
            gff_dir
            / f"{species2}.wgdi.gff",
            gff_dir
            / f"{species2}.gff",
            gff_dir
            / f"{species2}.wgdi.gff.tsv",
        ],
        f"WGDI GFF for {species2}",
    )

    lens1 = first_existing(
        [
            lens_dir
            / f"{species1}.wgdi.lens",
            lens_dir
            / f"{species1}.lens",
            lens_dir
            / f"{species1}.wgdi.lens.tsv",
        ],
        f"WGDI lens for {species1}",
    )

    lens2 = first_existing(
        [
            lens_dir
            / f"{species2}.wgdi.lens",
            lens_dir
            / f"{species2}.lens",
            lens_dir
            / f"{species2}.wgdi.lens.tsv",
        ],
        f"WGDI lens for {species2}",
    )

    config = (
        config_dir
        / f"{comparison}.blockinfo.conf"
    )

    result = (
        result_dir
        / f"{comparison}.blockinfo.csv"
    )

    stdout = (
        result_dir
        / f"{comparison}.blockinfo.wgdi.stdout.txt"
    )

    stderr = (
        result_dir
        / f"{comparison}.blockinfo.wgdi.stderr.txt"
    )

    qc = (
        qc_result_dir
        / f"{comparison}.blockinfo.qc.tsv"
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

    for path in [
        blast,
        gff1,
        gff2,
        lens1,
        lens2,
        collinearity,
        merged_ks,
    ]:
        if (
            not path.is_file()
            or path.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: Missing/empty input "
                f"for {comparison}: {path}"
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
            "collinearity": relative(
                collinearity
            ),
            "ks": relative(
                merged_ks
            ),
            "ks_col": "ks_YN00",
            "config": relative(
                config
            ),
            "result": relative(
                result
            ),
            "stdout": relative(
                stdout
            ),
            "stderr": relative(
                stderr
            ),
            "qc": relative(
                qc
            ),
            "work_dir": relative(
                comparison_work
            ),
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
            "blast_size_bytes": str(
                blast.stat().st_size
            ),
            "collinearity_size_bytes": str(
                collinearity.stat().st_size
            ),
            "ks_size_bytes": str(
                merged_ks.stat().st_size
            ),
            "gff1_size_bytes": str(
                gff1.stat().st_size
            ),
            "gff2_size_bytes": str(
                gff2.stat().st_size
            ),
            "lens1_size_bytes": str(
                lens1.stat().st_size
            ),
            "lens2_size_bytes": str(
                lens2.stat().st_size
            ),
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
        fieldnames=list(
            manifest_rows[0]
        ),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        manifest_rows
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

print(
    "Prepared 3 additional WGDI blockinfo "
    "configurations using ks_YN00."
)
PY

###############################################################################
# FINAL PREPARATION QC
###############################################################################

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

if [[ "${READY_COUNT}" -ne 3 ]]; then
    echo "ERROR: Only ${READY_COUNT}/3 configurations READY." >&2
    exit 1
fi

for COMPARISON in \
    VSCU_VSER \
    VSCU_VPAN \
    VPAN_VPER
do
    CONFIG="${CONFIG_DIR}/${COMPARISON}.blockinfo.conf"

    if [[ ! -s "${CONFIG}" ]]; then
        echo "ERROR: Missing config: ${CONFIG}" >&2
        exit 1
    fi

    for REQUIRED in \
        "[blockinfo]" \
        "score = 100" \
        "evalue = 1e-5" \
        "repeat_number = 20" \
        "position = order" \
        "ks_col = ks_YN00"
    do
        if ! grep -Fq \
            "${REQUIRED}" \
            "${CONFIG}"
        then
            echo "ERROR: Missing config setting '${REQUIRED}' in ${CONFIG}" >&2
            exit 1
        fi
    done
done

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36X5A_COMPLETE.txt" <<EOF2
checkpoint=step36X5A_prepare_additional_WGDI_blockinfo
date=$(date --iso-8601=seconds)
comparisons_expected=3
comparisons_ready=${READY_COUNT}
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
comparison_type=pairwise
parameter_set=strict
mg=25,25
primary_ks_column=ks_YN00
score=100
evalue=1e-5
repeat_number=20
position=order
manifest=11_wgdi/00_admin/step36X5_blockinfo_manifest.tsv
summary=11_wgdi/02_qc/additional_blockinfo_preparation/step36X5A_blockinfo_input_summary.tsv
original_step36E_outputs_modified=false
status=PASS
next_step=step36X5B_run_additional_WGDI_blockinfo
EOF2

cp -f \
    "${BLOCKINFO_MANIFEST}" \
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
echo "Step 36X5A preparation summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36X5 blockinfo manifest"
echo "============================================================"

column -t -s $'\t' \
    "${BLOCKINFO_MANIFEST}"

echo
echo "============================================================"
echo "Step 36X5A checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X5A_COMPLETE.txt"
