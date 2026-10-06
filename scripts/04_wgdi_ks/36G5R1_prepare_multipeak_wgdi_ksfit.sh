#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36G5R1
#SBATCH --output=11_wgdi/logs/step36G5R1_%j.out
#SBATCH --error=11_wgdi/logs/step36G5R1_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

BLOCKINFO_DIR="${WGDI_ROOT}/08_block_ks/01_curated_blockinfo"

VALIDATED_PEAKS="${WGDI_ROOT}/11_ksfigure/00_standardized_inputs/step36F4R_validated_peaks.standardized.tsv"

OUT_ROOT="${WGDI_ROOT}/11_ksfigure"
KSFIT_DIR="${OUT_ROOT}/06_multipeak_ksfit"
CONFIG_DIR="${OUT_ROOT}/07_multipeak_configs"
PLOT_DIR="${OUT_ROOT}/08_multipeak_plots"
DIAGNOSTIC_DIR="${OUT_ROOT}/09_multipeak_diagnostics"

QC_DIR="${WGDI_ROOT}/02_qc/ksfigure_multipeak"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36G5R1"

PARAMETER_TABLE="${QC_DIR}/step36G5R1_multipeak_parameters.tsv"
FIT_SUMMARY="${QC_DIR}/step36G5R1_multipeak_fit_summary.tsv"
RUN_MANIFEST="${WGDI_ROOT}/00_admin/step36G5R_multipeak_manifest.tsv"

mkdir -p \
    "${KSFIT_DIR}" \
    "${CONFIG_DIR}" \
    "${PLOT_DIR}" \
    "${DIAGNOSTIC_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

if [[ ! -s "${VALIDATED_PEAKS}" ]]; then
    echo "ERROR: Missing validated peak table: ${VALIDATED_PEAKS}" >&2
    exit 1
fi

rm -f \
    "${PARAMETER_TABLE}" \
    "${FIT_SUMMARY}" \
    "${RUN_MANIFEST}" \
    "${CHECKPOINT_DIR}/STEP36G5R1_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

find "${KSFIT_DIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.csv' \
    -delete

find "${CONFIG_DIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.conf' \
    -delete

find "${DIAGNOSTIC_DIR}" \
    -maxdepth 1 \
    -type f \
    \( -name '*.pdf' -o -name '*.svg' \) \
    -delete

module purge
module load hpc-env/13.1 2>/dev/null || true

CONDA_EXE_PATH=""

if [[ -n "${CONDA_EXE:-}" && -x "${CONDA_EXE}" ]]; then
    CONDA_EXE_PATH="${CONDA_EXE}"
elif command -v conda >/dev/null 2>&1; then
    CONDA_EXE_PATH="$(command -v conda)"
elif [[ -x "${HOME}/miniforge3/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/miniforge3/bin/conda"
elif [[ -x "${HOME}/mambaforge/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/mambaforge/bin/conda"
elif [[ -x "${HOME}/miniconda3/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/miniconda3/bin/conda"
else
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate wgdi_env

export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

python - \
    "${PROJECT_ROOT}" \
    "${BLOCKINFO_DIR}" \
    "${VALIDATED_PEAKS}" \
    "${KSFIT_DIR}" \
    "${CONFIG_DIR}" \
    "${PLOT_DIR}" \
    "${DIAGNOSTIC_DIR}" \
    "${PARAMETER_TABLE}" \
    "${FIT_SUMMARY}" \
    "${RUN_MANIFEST}" <<'PY'
from __future__ import annotations

import csv
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde

project_root = Path(sys.argv[1]).resolve()
blockinfo_dir = Path(sys.argv[2]).resolve()
validated_peaks_path = Path(sys.argv[3]).resolve()
ksfit_dir = Path(sys.argv[4]).resolve()
config_dir = Path(sys.argv[5]).resolve()
plot_dir = Path(sys.argv[6]).resolve()
diagnostic_dir = Path(sys.argv[7]).resolve()
parameter_table_path = Path(sys.argv[8])
fit_summary_path = Path(sys.argv[9])
manifest_path = Path(sys.argv[10])

comparisons = [
    "PMAJ_PMAJ",
    "VPAN_VPAN",
    "VSCU_VSCU",
    "VANA_VANA",
    "VARV_VARV",
    "VPER_VPER",
    "VSER_VSER",
    "VTRI_VTRI",
    "VVER_VVER",
    "PMAJ_VSCU",
    "VSCU_VANA",
    "VSCU_VPER",
]

self_comparisons = [
    comparison
    for comparison in comparisons
    if comparison.split("_")[0] == comparison.split("_")[1]
]

pairwise_comparisons = [
    comparison
    for comparison in comparisons
    if comparison not in self_comparisons
]

colors = {
    "PMAJ_PMAJ": "black",
    "VPAN_VPAN": "#1f77b4",
    "VSCU_VSCU": "#17becf",
    "VANA_VANA": "#d62728",
    "VARV_VARV": "#2ca02c",
    "VPER_VPER": "#9467bd",
    "VSER_VSER": "#8c564b",
    "VTRI_VTRI": "#e377c2",
    "VVER_VVER": "#7f7f7f",
    "PMAJ_VSCU": "#bcbd22",
    "VSCU_VANA": "#ff7f0e",
    "VSCU_VPER": "#4b0082",
}

linestyles = {
    comparison: (
        "--"
        if comparison in pairwise_comparisons
        or comparison == "PMAJ_PMAJ"
        else "-"
    )
    for comparison in comparisons
}


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
    return str(path.resolve().relative_to(project_root))


def first_present(
    row: dict[str, str],
    candidates: list[str],
) -> str:
    for candidate in candidates:
        value = row.get(candidate, "").strip()

        if value:
            return value

    raise KeyError(
        f"None of the expected columns was found: {candidates}"
    )


def multigaussian(
    x: np.ndarray,
    *parameters: float,
) -> np.ndarray:
    """
    Same functional form used by WGDI peaksfit:

        amplitude * exp(-((x - centre) / width)^2)
    """
    y = np.zeros_like(
        x,
        dtype=float,
    )

    for index in range(
        0,
        len(parameters),
        3,
    ):
        amplitude = parameters[index]
        centre = parameters[index + 1]
        width = parameters[index + 2]

        y += amplitude * np.exp(
            -((x - centre) / width) ** 2
        )

    return y


peak_rows = read_tsv(validated_peaks_path)

peaks_by_comparison: dict[
    str,
    list[dict[str, str]],
] = {
    comparison: []
    for comparison in comparisons
}

for row in peak_rows:
    comparison = first_present(
        row,
        ["comparison", "comparison_id"],
    )

    if comparison not in peaks_by_comparison:
        continue

    peak_value = float(
        first_present(
            row,
            [
                "peak_ks",
                "component_peak_ks",
                "validated_peak_ks",
                "peak",
            ],
        )
    )

    validation_class = row.get(
        "validation_class",
        row.get("status", "VALIDATED"),
    ).strip()

    row["_peak_value"] = str(peak_value)
    row["_validation_class"] = validation_class

    peaks_by_comparison[comparison].append(row)

for comparison in comparisons:
    if not peaks_by_comparison[comparison]:
        raise SystemExit(
            f"ERROR: No validated peaks found for {comparison}."
        )

    peaks_by_comparison[comparison].sort(
        key=lambda row: float(row["_peak_value"])
    )

parameter_rows: list[dict[str, str]] = []
fit_summary_rows: list[dict[str, str]] = []

fitted_parameters: dict[
    str,
    list[float],
] = {}

for comparison in comparisons:
    blockinfo_path = (
        blockinfo_dir
        / f"{comparison}.primary.filtered.blockinfo.csv"
    )

    if (
        not blockinfo_path.is_file()
        or blockinfo_path.stat().st_size == 0
    ):
        raise SystemExit(
            f"ERROR: Missing curated blockinfo: {blockinfo_path}"
        )

    blockinfo = pd.read_csv(blockinfo_path)

    if "ks_median" not in blockinfo.columns:
        raise SystemExit(
            f"ERROR: ks_median absent from {blockinfo_path}."
        )

    ks_values = pd.to_numeric(
        blockinfo["ks_median"],
        errors="coerce",
    ).to_numpy(dtype=float)

    ks_values = ks_values[
        np.isfinite(ks_values)
        & (ks_values > 0)
        & (ks_values <= 3)
    ]

    if len(ks_values) < 10:
        raise SystemExit(
            f"ERROR: Too few positive Ks values for {comparison}: "
            f"{len(ks_values)}"
        )

    x = np.linspace(
        0.0,
        3.0,
        2000,
    )

    kde = gaussian_kde(ks_values)

    # Same bandwidth reduction used in WGDI peaksfit.
    kde.set_bandwidth(
        bw_method=kde.factor / 3.0
    )

    target_density = kde(x)

    validated_centres = np.asarray(
        [
            float(row["_peak_value"])
            for row in peaks_by_comparison[comparison]
        ],
        dtype=float,
    )

    number_of_components = len(validated_centres)

    initial_parameters: list[float] = []
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []

    for peak_index, centre in enumerate(
        validated_centres
    ):
        peak_density = float(
            np.interp(
                centre,
                x,
                target_density,
            )
        )

        initial_amplitude = max(
            peak_density,
            0.01,
        )

        if number_of_components == 1:
            centre_half_window = 0.25
        else:
            neighbouring_distances = []

            if peak_index > 0:
                neighbouring_distances.append(
                    centre
                    - validated_centres[peak_index - 1]
                )

            if peak_index < number_of_components - 1:
                neighbouring_distances.append(
                    validated_centres[peak_index + 1]
                    - centre
                )

            nearest_distance = min(
                neighbouring_distances
            )

            centre_half_window = min(
                0.20,
                max(
                    0.03,
                    nearest_distance * 0.30,
                ),
            )

        initial_width = min(
            0.25,
            max(
                0.04,
                0.10,
            ),
        )

        initial_parameters.extend(
            [
                initial_amplitude,
                centre,
                initial_width,
            ]
        )

        lower_bounds.extend(
            [
                0.0,
                max(
                    0.001,
                    centre - centre_half_window,
                ),
                0.01,
            ]
        )

        upper_bounds.extend(
            [
                max(
                    20.0,
                    float(target_density.max()) * 10.0,
                ),
                min(
                    2.999,
                    centre + centre_half_window,
                ),
                1.0,
            ]
        )

    try:
        optimized_parameters, covariance = curve_fit(
            multigaussian,
            x,
            target_density,
            p0=np.asarray(
                initial_parameters,
                dtype=float,
            ),
            bounds=(
                np.asarray(
                    lower_bounds,
                    dtype=float,
                ),
                np.asarray(
                    upper_bounds,
                    dtype=float,
                ),
            ),
            maxfev=300000,
        )

    except Exception as error:
        raise SystemExit(
            f"ERROR: Multi-Gaussian fit failed for "
            f"{comparison}: {error}"
        ) from error

    fitted_curve = multigaussian(
        x,
        *optimized_parameters,
    )

    residual_sum_squares = float(
        np.sum(
            (target_density - fitted_curve) ** 2
        )
    )

    total_sum_squares = float(
        np.sum(
            (
                target_density
                - np.mean(target_density)
            )
            ** 2
        )
    )

    r_squared = (
        1.0
        - residual_sum_squares
        / total_sum_squares
        if total_sum_squares > 0
        else float("nan")
    )

    fitted_parameters[comparison] = (
        optimized_parameters.tolist()
    )

    for component_index in range(
        number_of_components
    ):
        offset = component_index * 3

        amplitude = float(
            optimized_parameters[offset]
        )

        centre = float(
            optimized_parameters[offset + 1]
        )

        width = float(
            optimized_parameters[offset + 2]
        )

        source_peak = peaks_by_comparison[
            comparison
        ][component_index]

        parameter_rows.append(
            {
                "comparison": comparison,
                "component_number": str(
                    component_index + 1
                ),
                "validation_class": (
                    source_peak["_validation_class"]
                ),
                "initial_validated_peak_ks": (
                    source_peak["_peak_value"]
                ),
                "fitted_amplitude": (
                    f"{amplitude:.12g}"
                ),
                "fitted_centre_ks": (
                    f"{centre:.12g}"
                ),
                "fitted_width": (
                    f"{width:.12g}"
                ),
                "kde_bandwidth_factor": (
                    f"{kde.factor:.12g}"
                ),
                "comparison_r_squared": (
                    f"{r_squared:.12g}"
                ),
                "status": "PASS",
            }
        )

    fit_summary_rows.append(
        {
            "comparison": comparison,
            "comparison_type": (
                "self"
                if comparison in self_comparisons
                else "pairwise"
            ),
            "syntenic_blocks_used": str(
                len(ks_values)
            ),
            "validated_components": str(
                number_of_components
            ),
            "validated_peak_centres": ",".join(
                f"{value:.8f}"
                for value in validated_centres
            ),
            "fitted_peak_centres": ",".join(
                f"{optimized_parameters[index]:.8f}"
                for index in range(
                    1,
                    len(optimized_parameters),
                    3,
                )
            ),
            "r_squared": f"{r_squared:.10f}",
            "kde_maximum": (
                f"{target_density.max():.10f}"
            ),
            "fitted_maximum": (
                f"{fitted_curve.max():.10f}"
            ),
            "blockinfo": relative(
                blockinfo_path
            ),
            "status": (
                "PASS"
                if math.isfinite(r_squared)
                else "FAIL"
            ),
        }
    )

    ###########################################################################
    # DIAGNOSTIC: KDE VERSUS MULTI-GAUSSIAN FIT
    ###########################################################################

    figure, axis = plt.subplots(
        figsize=(10, 6.18)
    )

    axis.plot(
        x,
        target_density,
        linewidth=2.0,
        label="Syntenic-block KDE",
    )

    axis.plot(
        x,
        fitted_curve,
        linewidth=2.0,
        linestyle="--",
        label="Multi-Gaussian fit",
    )

    for component_index in range(
        number_of_components
    ):
        offset = component_index * 3

        component_curve = multigaussian(
            x,
            *optimized_parameters[
                offset:offset + 3
            ],
        )

        axis.plot(
            x,
            component_curve,
            linewidth=1.0,
            alpha=0.7,
            label=(
                f"Component {component_index + 1}: "
                f"Ks={optimized_parameters[offset + 1]:.3f}"
            ),
        )

    axis.set_xlim(0, 3)
    axis.set_ylim(bottom=0)

    axis.set_xlabel(
        "Synonymous substitutions per synonymous site (Ks)"
    )

    axis.set_ylabel(
        "Kernel density of syntenic blocks"
    )

    axis.set_title(
        f"{comparison}: validated multi-Gaussian fit\n"
        f"blocks={len(ks_values)}; "
        f"components={number_of_components}; "
        f"R²={r_squared:.3f}"
    )

    axis.legend(
        fontsize=8,
        frameon=False,
    )

    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)

    figure.tight_layout()

    diagnostic_pdf = (
        diagnostic_dir
        / f"{comparison}.multipeak_fit_diagnostic.pdf"
    )

    diagnostic_svg = (
        diagnostic_dir
        / f"{comparison}.multipeak_fit_diagnostic.svg"
    )

    figure.savefig(
        diagnostic_pdf,
        bbox_inches="tight",
    )

    figure.savefig(
        diagnostic_svg,
        bbox_inches="tight",
    )

    plt.close(figure)

###############################################################################
# WRITE EXACT WGDI KSFIT TABLES
###############################################################################

maximum_parameter_count = max(
    len(values)
    for values in fitted_parameters.values()
)

if maximum_parameter_count % 3 != 0:
    raise SystemExit(
        "ERROR: Parameter count is not divisible by three."
    )

datasets = {
    "all_comparisons": comparisons,
    "self_comparisons": self_comparisons,
    "pairwise_comparisons": pairwise_comparisons,
}

manifest_rows: list[dict[str, str]] = []

for task_id, (
    dataset,
    selected_comparisons,
) in enumerate(datasets.items()):

    ksfit_path = (
        ksfit_dir
        / f"{dataset}.multipeak.ksfit.csv"
    )

    # Match the WGDI example exactly:
    # index,color,linewidth,linestyle,<unnamed Gaussian parameters>
    header = [
        "",
        "color",
        "linewidth",
        "linestyle",
    ] + [
        ""
        for _ in range(
            maximum_parameter_count
        )
    ]

    with ksfit_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.writer(
            handle,
            lineterminator="\n",
        )

        writer.writerow(header)

        for comparison in selected_comparisons:
            parameters = fitted_parameters[
                comparison
            ]

            padded_parameters: list[
                float | str
            ] = (
                parameters
                + [""]
                * (
                    maximum_parameter_count
                    - len(parameters)
                )
            )

            writer.writerow(
                [
                    comparison,
                    colors[comparison],
                    2.2,
                    linestyles[comparison],
                ]
                + padded_parameters
            )

    # Confirm that pandas will expose the blank headers as Unnamed:*.
    validation = pd.read_csv(
        ksfit_path,
        index_col=0,
    )

    unnamed_columns = [
        column
        for column in validation.columns
        if re.match(r"Unnamed:", str(column))
    ]

    if len(unnamed_columns) != maximum_parameter_count:
        raise SystemExit(
            f"ERROR: {dataset} has "
            f"{len(unnamed_columns)} Gaussian columns; "
            f"expected {maximum_parameter_count}."
        )

    for comparison, row in validation.iterrows():
        finite_parameters = pd.to_numeric(
            row[unnamed_columns],
            errors="coerce",
        ).dropna()

        if len(finite_parameters) < 3:
            raise SystemExit(
                f"ERROR: Too few parameters for {comparison}."
            )

        if len(finite_parameters) % 3 != 0:
            raise SystemExit(
                f"ERROR: Incomplete Gaussian triplet for "
                f"{comparison}."
            )

    config_path = (
        config_dir
        / f"{dataset}.total.conf"
    )

    pdf_output = (
        plot_dir
        / f"{dataset}.multipeak.ksfigure.pdf"
    )

    if dataset == "all_comparisons":
        title = "Ks distributions of syntenic blocks"
        figsize = "13,8"
        legendfontsize = "10"

    elif dataset == "self_comparisons":
        title = "Within-genome Ks distributions of syntenic blocks"
        figsize = "12,7.5"
        legendfontsize = "11"

    else:
        title = "Between-genome Ks distributions of syntenic blocks"
        figsize = "11,7"
        legendfontsize = "12"

    config_path.write_text(
        f"""[ksfigure]
ksfit = {ksfit_path}
labelfontsize = 16
legendfontsize = {legendfontsize}
xlabel = Synonymous substitutions per synonymous site (Ks)
ylabel = Kernel density of syntenic blocks
title = {title}
area = 0,3
shadow = true
figsize = {figsize}
savefig = {pdf_output}
""",
        encoding="utf-8",
    )

    manifest_rows.append(
        {
            "task_id": str(task_id),
            "dataset": dataset,
            "comparison_count": str(
                len(selected_comparisons)
            ),
            "ksfit": relative(ksfit_path),
            "config": relative(config_path),
            "pdf_output": relative(pdf_output),
            "svg_output": relative(
                plot_dir
                / f"{dataset}.multipeak.ksfigure.svg"
            ),
            "png_output": relative(
                plot_dir
                / f"{dataset}.multipeak.ksfigure.png"
            ),
            "status": "READY",
        }
    )

with parameter_table_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(parameter_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(parameter_rows)

with fit_summary_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(fit_summary_rows[0]),
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(fit_summary_rows)

with manifest_path.open(
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
    f"Comparisons fitted: {len(fit_summary_rows)}"
)

print(
    f"Gaussian components fitted: {len(parameter_rows)}"
)

print(
    f"WGDI KsFigure datasets prepared: {len(manifest_rows)}"
)
PY

PARAMETER_PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${PARAMETER_TABLE}"
)"

FIT_PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${FIT_SUMMARY}"
)"

READY_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "READY" {
            count++
        }
        END {
            print count + 0
        }
    ' "${RUN_MANIFEST}"
)"

if [[ "${PARAMETER_PASS_COUNT}" -lt 12 ]]; then
    echo "ERROR: Only ${PARAMETER_PASS_COUNT} Gaussian components passed." >&2
    exit 1
fi

if [[ "${FIT_PASS_COUNT}" -ne 12 ]]; then
    echo "ERROR: Only ${FIT_PASS_COUNT}/12 comparisons passed fitting." >&2
    exit 1
fi

if [[ "${READY_COUNT}" -ne 3 ]]; then
    echo "ERROR: Only ${READY_COUNT}/3 KsFigure datasets are ready." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36G5R1_COMPLETE.txt" <<EOF2
checkpoint=step36G5R1_prepare_multipeak_WGDI_KsFigure
date=$(date --iso-8601=seconds)
comparisons=12
source_values=primary_curated_blockinfo_ks_median
ks_range=0,3
density_method=scipy_gaussian_kde
bandwidth_method=WGDI_factor_divided_by_3
peak_centres_source=step36F4R_validated_peaks
gaussian_function=amplitude*exp(-((x-centre)/width)^2)
all_validated_components_retained=true
curve_normalization=false
y_axis=kernel_density_of_syntenic_blocks
ksfit_format=WGDI_example_blank_parameter_headers
prepared_datasets=3
parameter_table=11_wgdi/02_qc/ksfigure_multipeak/step36G5R1_multipeak_parameters.tsv
fit_summary=11_wgdi/02_qc/ksfigure_multipeak/step36G5R1_multipeak_fit_summary.tsv
run_manifest=11_wgdi/00_admin/step36G5R_multipeak_manifest.tsv
status=PASS
next_step=step36G5R2_run_authentic_WGDI_ksfigure
EOF2

cp -f \
    "${PARAMETER_TABLE}" \
    "${FIT_SUMMARY}" \
    "${RUN_MANIFEST}" \
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
echo "Multi-peak Gaussian parameters"
echo "============================================================"

column -t -s $'\t' \
    "${PARAMETER_TABLE}"

echo
echo "============================================================"
echo "Fit summary"
echo "============================================================"

column -t -s $'\t' \
    "${FIT_SUMMARY}"

echo
echo "============================================================"
echo "WGDI KsFigure manifest"
echo "============================================================"

column -t -s $'\t' \
    "${RUN_MANIFEST}"

echo
echo "============================================================"
echo "Checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36G5R1_COMPLETE.txt"
