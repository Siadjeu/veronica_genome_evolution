#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36F3Q
#SBATCH --output=11_wgdi/logs/step36F3Q_%j.out
#SBATCH --error=11_wgdi/logs/step36F3Q_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

INPUT_DIR="11_wgdi/08_block_ks/01_curated_blockinfo"
OUTPUT_DIR="11_wgdi/02_qc/block_ks_fit_audit"
CHECKPOINT_DIR="11_wgdi/checkpoints/step36F3Q"

SUMMARY="${OUTPUT_DIR}/step36F3Q_distribution_audit.tsv"
PEAKS="${OUTPUT_DIR}/step36F3Q_candidate_peaks.tsv"

mkdir -p \
  "${OUTPUT_DIR}" \
  "${CHECKPOINT_DIR}" \
  "11_wgdi/logs"

rm -f \
  "${SUMMARY}" \
  "${PEAKS}" \
  "${CHECKPOINT_DIR}/STEP36F3Q_COMPLETE.txt"

python - \
  "${INPUT_DIR}" \
  "${SUMMARY}" \
  "${PEAKS}" <<'PY'
from __future__ import annotations

import csv
import math
import statistics
import sys
from pathlib import Path

import numpy as np

try:
    from scipy.signal import find_peaks
    from scipy.stats import skew, kurtosis
except ImportError as error:
    raise SystemExit(
        "ERROR: scipy is required in the active environment."
    ) from error

input_dir = Path(sys.argv[1])
summary_path = Path(sys.argv[2])
peaks_path = Path(sys.argv[3])

files = sorted(
    input_dir.glob("*.primary.filtered.blockinfo.csv")
)

if len(files) != 12:
    raise SystemExit(
        f"ERROR: Expected 12 curated blockinfo files; found {len(files)}."
    )

summary_rows = []
peak_rows = []

BIN_WIDTH = 0.05
BIN_EDGES = np.arange(0.0, 3.0 + BIN_WIDTH, BIN_WIDTH)
BIN_CENTERS = (BIN_EDGES[:-1] + BIN_EDGES[1:]) / 2

for path in files:
    comparison = path.name.split(
        ".primary.filtered.blockinfo.csv"
    )[0]

    with path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        rows = list(csv.DictReader(handle))

    values = np.array(
        [
            float(row["ks_median"])
            for row in rows
            if 0.0 < float(row["ks_median"]) <= 3.0
        ],
        dtype=float,
    )

    if values.size == 0:
        raise SystemExit(
            f"ERROR: No positive Ks values for {comparison}."
        )

    histogram, _ = np.histogram(
        values,
        bins=BIN_EDGES,
    )

    smoothed = np.convolve(
        histogram.astype(float),
        np.array([1, 2, 3, 2, 1], dtype=float) / 9,
        mode="same",
    )

    minimum_prominence = max(
        2.0,
        0.05 * float(smoothed.max()),
    )

    peak_indices, properties = find_peaks(
        smoothed,
        prominence=minimum_prominence,
        distance=3,
    )

    ranked_peaks = sorted(
        zip(
            peak_indices,
            properties.get(
                "prominences",
                np.zeros(len(peak_indices)),
            ),
        ),
        key=lambda item: item[1],
        reverse=True,
    )

    histogram_mode_index = int(np.argmax(histogram))
    histogram_mode = float(BIN_CENTERS[histogram_mode_index])

    near_zero_count = int(np.sum(values <= 0.1))
    high_boundary_count = int(np.sum(values >= 2.9))

    q1, median, q3 = np.quantile(
        values,
        [0.25, 0.5, 0.75],
    )

    skewness = float(skew(values, bias=False))
    excess_kurtosis = float(
        kurtosis(values, fisher=True, bias=False)
    )

    if len(ranked_peaks) == 0:
        distribution_flag = "NO_CLEAR_INTERNAL_PEAK"
    elif len(ranked_peaks) == 1:
        distribution_flag = "ONE_CANDIDATE_PEAK"
    else:
        distribution_flag = "MULTIMODAL_CANDIDATE"

    if near_zero_count / values.size >= 0.20:
        boundary_flag = "STRONG_LOW_KS_BOUNDARY"
    elif high_boundary_count / values.size >= 0.10:
        boundary_flag = "STRONG_HIGH_KS_BOUNDARY"
    else:
        boundary_flag = "NO_STRONG_BOUNDARY"

    summary_rows.append(
        {
            "comparison": comparison,
            "block_count": int(values.size),
            "minimum_ks": f"{values.min():.6f}",
            "q1_ks": f"{q1:.6f}",
            "median_ks": f"{median:.6f}",
            "q3_ks": f"{q3:.6f}",
            "maximum_ks": f"{values.max():.6f}",
            "mean_ks": f"{values.mean():.6f}",
            "standard_deviation": (
                f"{values.std(ddof=1):.6f}"
                if values.size > 1
                else "NA"
            ),
            "skewness": f"{skewness:.6f}",
            "excess_kurtosis": f"{excess_kurtosis:.6f}",
            "histogram_mode_0.05_bins": f"{histogram_mode:.6f}",
            "values_ks_le_0.1": near_zero_count,
            "fraction_ks_le_0.1": (
                f"{near_zero_count / values.size:.8f}"
            ),
            "values_ks_ge_2.9": high_boundary_count,
            "fraction_ks_ge_2.9": (
                f"{high_boundary_count / values.size:.8f}"
            ),
            "candidate_peak_count": len(ranked_peaks),
            "distribution_flag": distribution_flag,
            "boundary_flag": boundary_flag,
            "single_gaussian_supported": (
                "NO"
                if (
                    len(ranked_peaks) != 1
                    or abs(skewness) > 1
                    or boundary_flag != "NO_STRONG_BOUNDARY"
                )
                else "POSSIBLY"
            ),
            "status": "PASS",
        }
    )

    for rank, (peak_index, prominence) in enumerate(
        ranked_peaks,
        start=1,
    ):
        peak_rows.append(
            {
                "comparison": comparison,
                "peak_rank": rank,
                "candidate_peak_ks": (
                    f"{BIN_CENTERS[peak_index]:.6f}"
                ),
                "smoothed_bin_height": (
                    f"{smoothed[peak_index]:.6f}"
                ),
                "prominence": f"{prominence:.6f}",
                "bin_width": f"{BIN_WIDTH:.6f}",
                "method": (
                    "smoothed_histogram_peak_detection"
                ),
            }
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

peak_fields = [
    "comparison",
    "peak_rank",
    "candidate_peak_ks",
    "smoothed_bin_height",
    "prominence",
    "bin_width",
    "method",
]

with peaks_path.open(
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
    writer.writerows(peak_rows)

print(
    f"Audited {len(summary_rows)} Ks distributions and "
    f"identified {len(peak_rows)} candidate peaks."
)
PY

cat > "${CHECKPOINT_DIR}/STEP36F3Q_COMPLETE.txt" <<EOF2
checkpoint=step36F3Q_audit_block_Ks_distributions
date=$(date --iso-8601=seconds)
comparisons=12
primary_data=filtered_block_median_ks_YN00
ks_range=(0,3]
histogram_bin_width=0.05
single_gaussian_not_assumed=true
summary=11_wgdi/02_qc/block_ks_fit_audit/step36F3Q_distribution_audit.tsv
candidate_peaks=11_wgdi/02_qc/block_ks_fit_audit/step36F3Q_candidate_peaks.tsv
status=PASS
next_step=review_distributions_before_gaussian_or_mixture_fitting
EOF2

echo
echo "============================================================"
echo "Distribution audit"
echo "============================================================"

column -t -s $'\t' "${SUMMARY}"

echo
echo "============================================================"
echo "Candidate peaks"
echo "============================================================"

column -t -s $'\t' "${PEAKS}"

echo
cat "${CHECKPOINT_DIR}/STEP36F3Q_COMPLETE.txt"
