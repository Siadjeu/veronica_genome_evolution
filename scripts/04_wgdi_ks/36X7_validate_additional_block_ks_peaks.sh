#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=5000
#SBATCH --job-name=wgdi36X7
#SBATCH --output=11_wgdi/logs/step36X7_%j.out
#SBATCH --error=11_wgdi/logs/step36X7_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

INPUT_MANIFEST="11_wgdi/00_admin/step36X6_curated_block_ks_manifest.tsv"

OUTPUT_ROOT="11_wgdi/08_block_ks/09_additional_corrected_peak_validation"
TABLE_DIR="${OUTPUT_ROOT}/tables"
PLOT_DIR="${OUTPUT_ROOT}/plots"

QC_DIR="11_wgdi/02_qc/additional_block_ks_peak_validation_corrected"
CHECKPOINT_DIR="11_wgdi/checkpoints/step36X7"

SUMMARY="${QC_DIR}/step36X7_peak_validation_summary.tsv"
PEAK_TABLE="${TABLE_DIR}/step36X7_validated_peaks.tsv"
BANDWIDTH_TABLE="${TABLE_DIR}/step36X7_bandwidth_peak_results.tsv"
BOOTSTRAP_TABLE="${TABLE_DIR}/step36X7_bootstrap_peak_support.tsv"
LOW_KS_TABLE="${TABLE_DIR}/step36X7_low_ks_component_summary.tsv"

mkdir -p \
    "${TABLE_DIR}" \
    "${PLOT_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "11_wgdi/logs"

if [[ ! -s "${INPUT_MANIFEST}" ]]; then
    echo "ERROR: Missing manifest: ${INPUT_MANIFEST}" >&2
    exit 1
fi

rm -f \
    "${SUMMARY}" \
    "${PEAK_TABLE}" \
    "${BANDWIDTH_TABLE}" \
    "${BOOTSTRAP_TABLE}" \
    "${LOW_KS_TABLE}" \
    "${CHECKPOINT_DIR}/STEP36X7_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

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
    "${INPUT_MANIFEST}" \
    "${PLOT_DIR}" \
    "${SUMMARY}" \
    "${PEAK_TABLE}" \
    "${BANDWIDTH_TABLE}" \
    "${BOOTSTRAP_TABLE}" \
    "${LOW_KS_TABLE}" <<'PY'
from __future__ import annotations

import csv
import math
import statistics
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.stats import gaussian_kde

project_root = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
plot_dir = Path(sys.argv[3])
summary_path = Path(sys.argv[4])
peak_table_path = Path(sys.argv[5])
bandwidth_table_path = Path(sys.argv[6])
bootstrap_table_path = Path(sys.argv[7])
low_ks_table_path = Path(sys.argv[8])

###############################################################################
# ANALYSIS SETTINGS
###############################################################################

WGD_KS_MIN = 0.10
WGD_KS_MAX = 3.00

BANDWIDTH_FACTORS = [0.6, 0.8, 1.0, 1.2, 1.5]

BOOTSTRAP_REPLICATES = 500
BOOTSTRAP_BANDWIDTH_FACTOR = 1.0
RANDOM_SEED = 20260805

MIN_BANDWIDTH_SUPPORT = 3
MIN_BOOTSTRAP_SUPPORT = 0.60

# Peak matching is performed in log-Ks space.
LOG_MATCH_TOLERANCE = 0.12

# Peaks must be separated by approximately 12% in log space.
MIN_LOG_PEAK_DISTANCE = 0.12

MIN_PROMINENCE_FRACTION = 0.05

LOG_GRID = np.linspace(
    math.log(WGD_KS_MIN),
    math.log(WGD_KS_MAX),
    1601,
)

KS_GRID = np.exp(LOG_GRID)
LOG_GRID_STEP = float(LOG_GRID[1] - LOG_GRID[0])

rng = np.random.default_rng(RANDOM_SEED)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_all_positive_block_medians(path: Path) -> np.ndarray:
    with path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        reader = csv.DictReader(handle)

        values = []

        for row in reader:
            raw = row.get("ks_median", "").strip()

            if not raw:
                continue

            value = float(raw)

            if math.isfinite(value) and value > 0:
                values.append(value)

    array = np.asarray(values, dtype=float)

    if array.size == 0:
        raise SystemExit(
            f"ERROR: No positive block-median Ks values in {path}."
        )

    return array


def log_kde(
    values: np.ndarray,
    bandwidth_factor: float,
) -> np.ndarray:
    analysis_values = values[
        (values >= WGD_KS_MIN)
        & (values <= WGD_KS_MAX)
    ]

    if analysis_values.size < 20:
        raise ValueError(
            "Fewer than 20 values in the WGD-analysis range."
        )

    log_values = np.log(analysis_values)

    if np.std(log_values, ddof=1) <= 0:
        raise ValueError(
            "Log-transformed values have zero variance."
        )

    kde = gaussian_kde(
        log_values,
        bw_method=lambda estimator: (
            estimator.scotts_factor()
            * bandwidth_factor
        ),
    )

    log_density = kde(LOG_GRID)

    # Convert density in log(Ks) back to density in Ks.
    ks_density = log_density / KS_GRID

    ks_density[~np.isfinite(ks_density)] = 0
    ks_density[ks_density < 0] = 0

    return ks_density


def detect_peaks(
    density: np.ndarray,
    observed_min: float,
    observed_max: float,
) -> list[dict[str, float]]:
    maximum = float(np.max(density))

    if maximum <= 0:
        return []

    minimum_distance = max(
        1,
        int(MIN_LOG_PEAK_DISTANCE / LOG_GRID_STEP),
    )

    indices, properties = find_peaks(
        density,
        prominence=maximum * MIN_PROMINENCE_FRACTION,
        distance=minimum_distance,
    )

    peaks = []

    for index, prominence in zip(
        indices,
        properties["prominences"],
    ):
        location = float(KS_GRID[index])

        # Do not accept extrapolated or boundary peaks.
        if location < observed_min or location > observed_max:
            continue

        # Require a peak to be genuinely internal rather than exactly
        # at the fixed analysis cutoffs.
        if location <= WGD_KS_MIN * 1.01:
            continue

        if location >= WGD_KS_MAX / 1.01:
            continue

        peaks.append(
            {
                "location": location,
                "log_location": math.log(location),
                "height": float(density[index]),
                "prominence": float(prominence),
            }
        )

    peaks.sort(
        key=lambda item: item["prominence"],
        reverse=True,
    )

    return peaks


def cluster_log_locations(
    log_locations: list[float],
) -> list[list[float]]:
    if not log_locations:
        return []

    values = sorted(log_locations)
    clusters = [[values[0]]]

    for value in values[1:]:
        center = statistics.median(clusters[-1])

        if abs(value - center) <= LOG_MATCH_TOLERANCE:
            clusters[-1].append(value)
        else:
            clusters.append([value])

    return clusters


manifest_rows = read_tsv(manifest_path)

if len(manifest_rows) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 comparisons; found {len(manifest_rows)}."
    )

summary_rows = []
validated_peak_rows = []
bandwidth_rows = []
bootstrap_rows = []
low_ks_rows = []

for manifest_row in manifest_rows:
    comparison = manifest_row["comparison"]
    comparison_type = manifest_row["comparison_type"]
    species1 = manifest_row["species1"]
    species2 = manifest_row["species2"]

    dataset_paths = {
        "primary": (
            project_root
            / manifest_row["primary_blockinfo"]
        ),
        "low_tandem": (
            project_root
            / manifest_row["low_tandem_blockinfo"]
        ),
    }

    dataset_results = {}

    for dataset_name, dataset_path in dataset_paths.items():
        all_values = load_all_positive_block_medians(
            dataset_path
        )

        analysis_values = all_values[
            (all_values >= WGD_KS_MIN)
            & (all_values <= WGD_KS_MAX)
        ]

        low_values = all_values[
            all_values < WGD_KS_MIN
        ]

        above_values = all_values[
            all_values > WGD_KS_MAX
        ]

        if analysis_values.size < 20:
            raise SystemExit(
                f"ERROR: Only {analysis_values.size} usable values "
                f"for {comparison}, {dataset_name}."
            )

        observed_min = float(analysis_values.min())
        observed_max = float(analysis_values.max())

        low_ks_rows.append(
            {
                "comparison": comparison,
                "dataset": dataset_name,
                "all_positive_blocks": str(all_values.size),
                "blocks_ks_below_0.1": str(low_values.size),
                "fraction_ks_below_0.1": (
                    f"{low_values.size / all_values.size:.8f}"
                ),
                "blocks_ks_0.1_to_3": str(
                    analysis_values.size
                ),
                "fraction_ks_0.1_to_3": (
                    f"{analysis_values.size / all_values.size:.8f}"
                ),
                "blocks_ks_above_3": str(above_values.size),
                "minimum_positive_ks": (
                    f"{all_values.min():.8f}"
                ),
                "minimum_analysis_ks": (
                    f"{observed_min:.8f}"
                ),
                "maximum_analysis_ks": (
                    f"{observed_max:.8f}"
                ),
            }
        )

        bandwidth_peak_sets = {}
        densities = {}
        all_log_locations = []

        for bandwidth_factor in BANDWIDTH_FACTORS:
            density = log_kde(
                analysis_values,
                bandwidth_factor,
            )

            peaks = detect_peaks(
                density,
                observed_min,
                observed_max,
            )

            densities[bandwidth_factor] = density
            bandwidth_peak_sets[bandwidth_factor] = peaks

            for rank, peak in enumerate(peaks, start=1):
                all_log_locations.append(
                    peak["log_location"]
                )

                bandwidth_rows.append(
                    {
                        "comparison": comparison,
                        "dataset": dataset_name,
                        "bandwidth_factor": (
                            f"{bandwidth_factor:.2f}"
                        ),
                        "peak_rank": str(rank),
                        "peak_ks": (
                            f"{peak['location']:.8f}"
                        ),
                        "peak_log_ks": (
                            f"{peak['log_location']:.8f}"
                        ),
                        "peak_height": (
                            f"{peak['height']:.10f}"
                        ),
                        "prominence": (
                            f"{peak['prominence']:.10f}"
                        ),
                    }
                )

        clusters = cluster_log_locations(
            all_log_locations
        )

        references = []

        for cluster in clusters:
            log_center = float(statistics.median(cluster))
            center = math.exp(log_center)

            bandwidth_support = 0

            for bandwidth_factor in BANDWIDTH_FACTORS:
                matched = any(
                    abs(
                        peak["log_location"]
                        - log_center
                    )
                    <= LOG_MATCH_TOLERANCE
                    for peak in bandwidth_peak_sets[
                        bandwidth_factor
                    ]
                )

                if matched:
                    bandwidth_support += 1

            if bandwidth_support >= MIN_BANDWIDTH_SUPPORT:
                references.append(
                    {
                        "location": center,
                        "log_location": log_center,
                        "bandwidth_support": bandwidth_support,
                    }
                )

        references.sort(
            key=lambda item: item["location"]
        )

        recoveries = [
            []
            for _ in references
        ]

        bootstrap_peak_counts = []

        for _ in range(BOOTSTRAP_REPLICATES):
            sample = rng.choice(
                analysis_values,
                size=analysis_values.size,
                replace=True,
            )

            sample_min = float(sample.min())
            sample_max = float(sample.max())

            density = log_kde(
                sample,
                BOOTSTRAP_BANDWIDTH_FACTOR,
            )

            sample_peaks = detect_peaks(
                density,
                sample_min,
                sample_max,
            )

            bootstrap_peak_counts.append(
                len(sample_peaks)
            )

            for reference_index, reference in enumerate(
                references
            ):
                matches = [
                    peak
                    for peak in sample_peaks
                    if abs(
                        peak["log_location"]
                        - reference["log_location"]
                    )
                    <= LOG_MATCH_TOLERANCE
                ]

                if matches:
                    closest = min(
                        matches,
                        key=lambda peak: abs(
                            peak["log_location"]
                            - reference["log_location"]
                        ),
                    )

                    recoveries[reference_index].append(
                        closest["location"]
                    )

        validated = []

        for rank, reference in enumerate(
            references,
            start=1,
        ):
            observed_recoveries = recoveries[rank - 1]

            support = (
                len(observed_recoveries)
                / BOOTSTRAP_REPLICATES
            )

            if observed_recoveries:
                bootstrap_median = float(
                    np.median(observed_recoveries)
                )

                lower = float(
                    np.quantile(
                        observed_recoveries,
                        0.025,
                    )
                )

                upper = float(
                    np.quantile(
                        observed_recoveries,
                        0.975,
                    )
                )
            else:
                bootstrap_median = math.nan
                lower = math.nan
                upper = math.nan

            is_valid = (
                support >= MIN_BOOTSTRAP_SUPPORT
                and math.isfinite(bootstrap_median)
                and observed_min
                <= bootstrap_median
                <= observed_max
            )

            bootstrap_rows.append(
                {
                    "comparison": comparison,
                    "dataset": dataset_name,
                    "reference_peak_rank": str(rank),
                    "reference_peak_ks": (
                        f"{reference['location']:.8f}"
                    ),
                    "bandwidth_support_count": str(
                        reference["bandwidth_support"]
                    ),
                    "bandwidths_tested": str(
                        len(BANDWIDTH_FACTORS)
                    ),
                    "bootstrap_recoveries": str(
                        len(observed_recoveries)
                    ),
                    "bootstrap_replicates": str(
                        BOOTSTRAP_REPLICATES
                    ),
                    "bootstrap_support_fraction": (
                        f"{support:.8f}"
                    ),
                    "bootstrap_peak_median": (
                        f"{bootstrap_median:.8f}"
                        if math.isfinite(
                            bootstrap_median
                        )
                        else "NA"
                    ),
                    "bootstrap_95pct_lower": (
                        f"{lower:.8f}"
                        if math.isfinite(lower)
                        else "NA"
                    ),
                    "bootstrap_95pct_upper": (
                        f"{upper:.8f}"
                        if math.isfinite(upper)
                        else "NA"
                    ),
                    "observed_minimum_ks": (
                        f"{observed_min:.8f}"
                    ),
                    "observed_maximum_ks": (
                        f"{observed_max:.8f}"
                    ),
                    "validated": (
                        "YES" if is_valid else "NO"
                    ),
                }
            )

            if is_valid:
                validated.append(
                    {
                        "location": bootstrap_median,
                        "lower": lower,
                        "upper": upper,
                        "support": support,
                        "bandwidth_support": (
                            reference["bandwidth_support"]
                        ),
                    }
                )

        dataset_results[dataset_name] = {
            "all_values": all_values,
            "analysis_values": analysis_values,
            "low_values": low_values,
            "observed_min": observed_min,
            "observed_max": observed_max,
            "densities": densities,
            "validated": validated,
            "bootstrap_peak_mean": float(
                np.mean(bootstrap_peak_counts)
            ),
        }

        #######################################################################
        # PLOT
        #######################################################################

        figure, axis = plt.subplots(
            figsize=(9, 5.5)
        )

        bins = np.geomspace(
            WGD_KS_MIN,
            WGD_KS_MAX,
            46,
        )

        axis.hist(
            analysis_values,
            bins=bins,
            density=True,
            alpha=0.25,
            label="Block-median Ks",
        )

        for bandwidth_factor in BANDWIDTH_FACTORS:
            axis.plot(
                KS_GRID,
                densities[bandwidth_factor],
                linewidth=1.2,
                label=(
                    f"KDE bandwidth × "
                    f"{bandwidth_factor:.1f}"
                ),
            )

        for peak in validated:
            axis.axvline(
                peak["location"],
                linestyle="--",
                linewidth=1.2,
            )

            axis.axvspan(
                peak["lower"],
                peak["upper"],
                alpha=0.12,
            )

        axis.set_xscale("log")
        axis.set_xlim(WGD_KS_MIN, WGD_KS_MAX)

        axis.set_xlabel(
            "Block-median Ks (YN00; logarithmic scale)"
        )
        axis.set_ylabel("Density")

        axis.set_title(
            f"{comparison}: "
            f"{dataset_name.replace('_', ' ')}"
        )

        axis.legend(
            frameon=False,
            fontsize=8,
        )

        figure.tight_layout()

        pdf_path = (
            plot_dir
            / f"{comparison}.{dataset_name}.corrected_kde.pdf"
        )

        svg_path = (
            plot_dir
            / f"{comparison}.{dataset_name}.corrected_kde.svg"
        )

        figure.savefig(
            pdf_path,
            bbox_inches="tight",
        )

        figure.savefig(
            svg_path,
            bbox_inches="tight",
        )

        plt.close(figure)

    ###########################################################################
    # CROSS-DATASET ROBUSTNESS
    ###########################################################################

    primary_peaks = dataset_results[
        "primary"
    ]["validated"]

    low_tandem_peaks = dataset_results[
        "low_tandem"
    ]["validated"]

    robust_count = 0

    for rank, primary_peak in enumerate(
        primary_peaks,
        start=1,
    ):
        primary_log = math.log(
            primary_peak["location"]
        )

        matches = [
            peak
            for peak in low_tandem_peaks
            if abs(
                math.log(peak["location"])
                - primary_log
            )
            <= LOG_MATCH_TOLERANCE
        ]

        reproduced = bool(matches)

        if reproduced:
            closest = min(
                matches,
                key=lambda peak: abs(
                    math.log(peak["location"])
                    - primary_log
                ),
            )

            low_tandem_location = closest["location"]
            validation_class = "ROBUST"
            robust_count += 1
        else:
            low_tandem_location = math.nan
            validation_class = "PRIMARY_ONLY"

        # Final impossibility check.
        primary_observed_min = dataset_results[
            "primary"
        ]["observed_min"]

        primary_observed_max = dataset_results[
            "primary"
        ]["observed_max"]

        if not (
            primary_observed_min
            <= primary_peak["location"]
            <= primary_observed_max
        ):
            raise SystemExit(
                f"ERROR: Peak outside observed range for "
                f"{comparison}: {primary_peak['location']}"
            )

        validated_peak_rows.append(
            {
                "comparison": comparison,
                "comparison_type": comparison_type,
                "species1": species1,
                "species2": species2,
                "peak_rank": str(rank),
                "primary_peak_ks": (
                    f"{primary_peak['location']:.8f}"
                ),
                "primary_bootstrap_95pct_lower": (
                    f"{primary_peak['lower']:.8f}"
                ),
                "primary_bootstrap_95pct_upper": (
                    f"{primary_peak['upper']:.8f}"
                ),
                "primary_bootstrap_support": (
                    f"{primary_peak['support']:.8f}"
                ),
                "primary_bandwidth_support": str(
                    primary_peak["bandwidth_support"]
                ),
                "primary_observed_minimum_ks": (
                    f"{primary_observed_min:.8f}"
                ),
                "primary_observed_maximum_ks": (
                    f"{primary_observed_max:.8f}"
                ),
                "recovered_in_low_tandem": (
                    "YES" if reproduced else "NO"
                ),
                "low_tandem_peak_ks": (
                    f"{low_tandem_location:.8f}"
                    if math.isfinite(
                        low_tandem_location
                    )
                    else "NA"
                ),
                "validation_class": validation_class,
            }
        )

    summary_rows.append(
        {
            "comparison": comparison,
            "comparison_type": comparison_type,
            "species1": species1,
            "species2": species2,
            "primary_all_positive_blocks": str(
                dataset_results[
                    "primary"
                ]["all_values"].size
            ),
            "primary_blocks_below_0.1": str(
                dataset_results[
                    "primary"
                ]["low_values"].size
            ),
            "primary_blocks_used_0.1_to_3": str(
                dataset_results[
                    "primary"
                ]["analysis_values"].size
            ),
            "low_tandem_blocks_used_0.1_to_3": str(
                dataset_results[
                    "low_tandem"
                ]["analysis_values"].size
            ),
            "primary_validated_peak_count": str(
                len(primary_peaks)
            ),
            "low_tandem_validated_peak_count": str(
                len(low_tandem_peaks)
            ),
            "robust_peak_count": str(robust_count),
            "primary_mean_bootstrap_peak_count": (
                f"{dataset_results['primary']['bootstrap_peak_mean']:.6f}"
            ),
            "low_tandem_mean_bootstrap_peak_count": (
                f"{dataset_results['low_tandem']['bootstrap_peak_mean']:.6f}"
            ),
            "artificial_boundary_reflection": "NO",
            "single_gaussian_recommended": "NO",
            "status": "PASS",
        }
    )

###############################################################################
# WRITE TABLES
###############################################################################

summary_fields = [
    "comparison",
    "comparison_type",
    "species1",
    "species2",
    "primary_all_positive_blocks",
    "primary_blocks_below_0.1",
    "primary_blocks_used_0.1_to_3",
    "low_tandem_blocks_used_0.1_to_3",
    "primary_validated_peak_count",
    "low_tandem_validated_peak_count",
    "robust_peak_count",
    "primary_mean_bootstrap_peak_count",
    "low_tandem_mean_bootstrap_peak_count",
    "artificial_boundary_reflection",
    "single_gaussian_recommended",
    "status",
]

with summary_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=summary_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(summary_rows)

peak_fields = [
    "comparison",
    "comparison_type",
    "species1",
    "species2",
    "peak_rank",
    "primary_peak_ks",
    "primary_bootstrap_95pct_lower",
    "primary_bootstrap_95pct_upper",
    "primary_bootstrap_support",
    "primary_bandwidth_support",
    "primary_observed_minimum_ks",
    "primary_observed_maximum_ks",
    "recovered_in_low_tandem",
    "low_tandem_peak_ks",
    "validation_class",
]

with peak_table_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=peak_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(validated_peak_rows)

bandwidth_fields = [
    "comparison",
    "dataset",
    "bandwidth_factor",
    "peak_rank",
    "peak_ks",
    "peak_log_ks",
    "peak_height",
    "prominence",
]

with bandwidth_table_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=bandwidth_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(bandwidth_rows)

bootstrap_fields = [
    "comparison",
    "dataset",
    "reference_peak_rank",
    "reference_peak_ks",
    "bandwidth_support_count",
    "bandwidths_tested",
    "bootstrap_recoveries",
    "bootstrap_replicates",
    "bootstrap_support_fraction",
    "bootstrap_peak_median",
    "bootstrap_95pct_lower",
    "bootstrap_95pct_upper",
    "observed_minimum_ks",
    "observed_maximum_ks",
    "validated",
]

with bootstrap_table_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=bootstrap_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(bootstrap_rows)

low_ks_fields = [
    "comparison",
    "dataset",
    "all_positive_blocks",
    "blocks_ks_below_0.1",
    "fraction_ks_below_0.1",
    "blocks_ks_0.1_to_3",
    "fraction_ks_0.1_to_3",
    "blocks_ks_above_3",
    "minimum_positive_ks",
    "minimum_analysis_ks",
    "maximum_analysis_ks",
]

with low_ks_table_path.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=low_ks_fields,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(low_ks_rows)

print(
    f"Corrected peak validation completed for "
    f"{len(summary_rows)} comparisons."
)
print(
    f"Validated primary peaks: "
    f"{len(validated_peak_rows)}"
)
PY

PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            count++
        }
        END {
            print count + 0
        }
    ' "${SUMMARY}"
)"

if [[ "${PASS_COUNT}" -ne 3 ]]; then
    echo "ERROR: Only ${PASS_COUNT}/3 comparisons passed." >&2
    exit 1
fi

# Preserve the original analysis but explicitly mark it superseded.
cat > \
  "11_wgdi/checkpoints/step36F4/STEP36F4_SUPERSEDED.txt" <<EOF2
original_checkpoint=step36F4_robust_block_Ks_peak_validation
status=SUPERSEDED
reason=lower-boundary reflection generated peaks below the observed data range
examples=PMAJ_VSCU_peak_0.1_below_observed_minimum_0.6476;VSCU_VANA_peak_0.0_below_observed_minimum_0.2019
replacement=step36F4R_corrected_log_Ks_peak_validation
do_not_use_original_peak_table_for_inference=true
EOF2

cat > "${CHECKPOINT_DIR}/STEP36X7_COMPLETE.txt" <<EOF2
checkpoint=step36X7_additional_corrected_block_Ks_peak_validation
date=$(date --iso-8601=seconds)
comparisons=12
primary_ks_column=ks_YN00
wgd_candidate_ks_range=[0.1,3]
values_below_0.1=summarized_separately
transformation=natural_log_Ks
artificial_boundary_reflection=false
kde_bandwidth_factors=0.6,0.8,1.0,1.2,1.5
bootstrap_replicates=500
bootstrap_peak_support_threshold=0.60
bandwidth_support_threshold=3_of_5
log_peak_matching_tolerance=0.12
single_gaussian_assumed=false
supersedes=step36F4
summary=11_wgdi/02_qc/additional_block_ks_peak_validation_corrected/step36X7_peak_validation_summary.tsv
validated_peaks=11_wgdi/08_block_ks/09_additional_corrected_peak_validation/tables/step36X7_validated_peaks.tsv
low_ks_summary=11_wgdi/08_block_ks/09_additional_corrected_peak_validation/tables/step36X7_low_ks_component_summary.tsv
plots_directory=11_wgdi/08_block_ks/09_additional_corrected_peak_validation/plots
status=PASS
next_step=interpret_corrected_robust_peaks_and_compare_self_pairwise_signals
EOF2

cp -f \
    "${SUMMARY}" \
    "${PEAK_TABLE}" \
    "${LOW_KS_TABLE}" \
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
echo "Corrected peak-validation summary"
echo "============================================================"

column -t -s $'\t' "${SUMMARY}"

echo
echo "============================================================"
echo "Corrected validated peaks"
echo "============================================================"

if [[ "$(wc -l < "${PEAK_TABLE}")" -gt 1 ]]; then
    column -t -s $'\t' "${PEAK_TABLE}"
else
    cat "${PEAK_TABLE}"
    echo "No primary peaks passed all validation criteria."
fi

echo
echo "============================================================"
echo "Low-Ks component summary"
echo "============================================================"

column -t -s $'\t' "${LOW_KS_TABLE}"

echo
echo "============================================================"
echo "Checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36X7_COMPLETE.txt"
