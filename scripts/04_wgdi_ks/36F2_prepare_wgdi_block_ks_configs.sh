#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=3000
#SBATCH --job-name=wgdi36F2
#SBATCH --output=11_wgdi/logs/step36F2_%j.out
#SBATCH --error=11_wgdi/logs/step36F2_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

CURATED_MANIFEST="${WGDI_ROOT}/00_admin/step36F_curated_block_ks_manifest.tsv"
BLOCKINFO_MANIFEST="${WGDI_ROOT}/00_admin/step36E_blockinfo_manifest.tsv"

OUT_ROOT="${WGDI_ROOT}/08_block_ks"
CONFIG_DIR="${OUT_ROOT}/03_configs"
PLOT_DIR="${OUT_ROOT}/04_plots"
DATA_DIR="${OUT_ROOT}/05_peak_tables"
WORK_DIR="${OUT_ROOT}/06_work"

ADMIN_DIR="${WGDI_ROOT}/00_admin"
QC_DIR="${WGDI_ROOT}/02_qc/block_ks_config"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36F2"

PLOT_MANIFEST="${ADMIN_DIR}/step36F_plot_manifest.tsv"
SUMMARY="${QC_DIR}/step36F2_config_summary.tsv"

mkdir -p \
    "${CONFIG_DIR}" \
    "${PLOT_DIR}" \
    "${DATA_DIR}" \
    "${WORK_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

for FILE in \
    "${CURATED_MANIFEST}" \
    "${BLOCKINFO_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing manifest: ${FILE}" >&2
        exit 1
    fi
done

rm -f \
    "${PLOT_MANIFEST}" \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36F2_COMPLETE.txt"

python - \
    "${PROJECT_ROOT}" \
    "${CURATED_MANIFEST}" \
    "${BLOCKINFO_MANIFEST}" \
    "${CONFIG_DIR}" \
    "${PLOT_DIR}" \
    "${DATA_DIR}" \
    "${WORK_DIR}" \
    "${PLOT_MANIFEST}" \
    "${SUMMARY}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

project_root = Path(sys.argv[1])
curated_manifest = Path(sys.argv[2])
blockinfo_manifest = Path(sys.argv[3])
config_dir = Path(sys.argv[4])
plot_dir = Path(sys.argv[5])
data_dir = Path(sys.argv[6])
work_dir = Path(sys.argv[7])
plot_manifest = Path(sys.argv[8])
summary_path = Path(sys.argv[9])


def read_tsv(path: Path):
    with path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


curated_rows = {
    row["comparison"]: row
    for row in read_tsv(curated_manifest)
}

input_rows = {
    row["comparison"]: row
    for row in read_tsv(blockinfo_manifest)
}

if set(curated_rows) != set(input_rows):
    raise SystemExit(
        "ERROR: Step 36E and Step 36F1 comparison sets differ."
    )

manifest_rows = []
summary_rows = []

for task_id, comparison in enumerate(sorted(curated_rows)):
    curated = curated_rows[comparison]
    source = input_rows[comparison]

    primary_blockinfo = (
        project_root / curated["primary_blockinfo"]
    )

    lens1 = project_root / source["lens1"]
    lens2 = project_root / source["lens2"]

    for required in [primary_blockinfo, lens1, lens2]:
        if not required.is_file() or required.stat().st_size == 0:
            raise SystemExit(
                f"ERROR: Missing input for {comparison}: {required}"
            )

    species1 = curated["species1"]
    species2 = curated["species2"]

    bk_config = config_dir / f"{comparison}.blockks.conf"
    kp_config = config_dir / f"{comparison}.kspeaks.conf"
    pf_config = config_dir / f"{comparison}.peaksfit.conf"

    bk_plot = plot_dir / f"{comparison}.blockks.pdf"
    kp_plot = plot_dir / f"{comparison}.kspeaks.pdf"
    pf_plot = plot_dir / f"{comparison}.peaksfit.pdf"

    kp_table = data_dir / f"{comparison}.kspeaks.csv"

    bk_config.write_text(
        "\n".join(
            [
                "[blockks]",
                f"lens1 = {lens1}",
                f"lens2 = {lens2}",
                f"genome1_name = {species1}",
                f"genome2_name = {species2}",
                f"blockinfo = {primary_blockinfo}",
                "pvalue = 0.2",
                "tandem = true",
                "tandem_length = 200",
                "markersize = 1",
                "area = 0,3",
                "block_length = 10",
                "figsize = 10,10",
                f"savefig = {bk_plot}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    kp_config.write_text(
        "\n".join(
            [
                "[kspeaks]",
                f"blockinfo = {primary_blockinfo}",
                "pvalue = 0.2",
                "tandem = true",
                "block_length = 10",
                "ks_area = 0,3",
                "multiple = 1",
                "homo = 0,1",
                "fontsize = 10",
                "area = 0,3",
                "figsize = 10,6.18",
                f"savefig = {kp_plot}",
                f"savefile = {kp_table}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    pf_config.write_text(
        "\n".join(
            [
                "[peaksfit]",
                f"blockinfo = {primary_blockinfo}",
                "mode = median",
                "bins_number = 200",
                "ks_area = 0,3",
                "fontsize = 10",
                "area = 0,3",
                "figsize = 10,6.18",
                "shadow = true",
                f"savefig = {pf_plot}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    manifest_rows.append(
        {
            "task_id": str(task_id),
            "comparison": comparison,
            "comparison_type": curated["comparison_type"],
            "species1": species1,
            "species2": species2,
            "blockinfo": relative(primary_blockinfo),
            "bk_config": relative(bk_config),
            "kp_config": relative(kp_config),
            "pf_config": relative(pf_config),
            "bk_plot": relative(bk_plot),
            "kp_plot": relative(kp_plot),
            "pf_plot": relative(pf_plot),
            "kp_table": relative(kp_table),
            "work_dir": relative(work_dir / comparison),
            "status": "READY",
        }
    )

    summary_rows.append(
        {
            "comparison": comparison,
            "primary_blocks": curated["primary_block_count"],
            "ks_area": "0,3",
            "minimum_block_length": "10",
            "pvalue_max": "0.2",
            "configs_created": "3",
            "status": "PASS",
        }
    )

with plot_manifest.open(
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

print("Prepared 36 WGDI plotting configurations.")
PY

cat > "${CHECKPOINT_DIR}/STEP36F2_COMPLETE.txt" <<EOF2
checkpoint=step36F2_prepare_WGDI_block_Ks_configs
date=$(date --iso-8601=seconds)
comparisons=12
configs_per_comparison=3
total_configs=36
modules=blockks,kspeaks,peaksfit
ks_area=0,3
minimum_block_length=10
manifest=11_wgdi/00_admin/step36F_plot_manifest.tsv
status=PASS
next_step=step36F3_run_WGDI_block_Ks_modules
EOF2

echo
column -t -s $'\t' "${SUMMARY}"

echo
cat "${CHECKPOINT_DIR}/STEP36F2_COMPLETE.txt"
