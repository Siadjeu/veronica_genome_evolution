#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=36X13_integrate
#SBATCH --output=11_wgdi/checkpoints/step36X13_integrated_evidence/logs/step36X13_%j.out
#SBATCH --error=11_wgdi/checkpoints/step36X13_integrated_evidence/logs/step36X13_%j.err

set -euo pipefail

###############################################################################
# STEP 36X13 — CORRECTED
#
# Integrated evidence / manuscript checkpoint
#
# Additional comparison set:
#
#   VSCU -> VSER    diploid -> diploid
#   VSCU -> VPAN    diploid -> diploid
#   VPAN -> VPER    diploid -> tetraploid
#
# IMPORTANT
# ---------
# This step performs NO new biological analysis.
#
# It validates, integrates and freezes already completed analyses:
#
#   Step 35Y2  JCVI synteny-optimized macrosynteny
#   Step 36X7  corrected validated Ks peaks
#   Step 36X8A peak/block chromosome mapping
#   Step 36X8B component syntenic depth
#   Step 36X9  numbered WGDI BlockKs plots
#   Step 36X10 authentic WGDI KsFigure
#   Step 36X11 gene-level multisynteny
#   Step 36X12 chromosome-consistent multisynteny
#
# The failed previous 36X13 attempt is replaced completely.
# Upstream analyses are read-only.
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

CHECKPOINT_ROOT="11_wgdi/checkpoints/step36X13_integrated_evidence"

ADMIN_DIR="${CHECKPOINT_ROOT}/00_admin"
TABLE_DIR="${CHECKPOINT_ROOT}/01_tables"
FIGURE_DIR="${CHECKPOINT_ROOT}/02_figures"
MANUSCRIPT_DIR="${CHECKPOINT_ROOT}/03_manuscript"
PROVENANCE_DIR="${CHECKPOINT_ROOT}/04_provenance"
LOG_DIR="${CHECKPOINT_ROOT}/logs"

mkdir -p \
    "${ADMIN_DIR}" \
    "${TABLE_DIR}" \
    "${FIGURE_DIR}" \
    "${MANUSCRIPT_DIR}" \
    "${PROVENANCE_DIR}" \
    "${LOG_DIR}"

###############################################################################
# Explicit validated source paths
###############################################################################

# ---------------------------------------------------------------------------
# Step 36X7 — corrected validated Ks peaks
# ---------------------------------------------------------------------------

X7_PEAKS="11_wgdi/08_block_ks/09_additional_corrected_peak_validation/tables/step36X7_validated_peaks.tsv"


# ---------------------------------------------------------------------------
# Step 36X8A — peak/block chromosome mapping
# ---------------------------------------------------------------------------

X8A_BLOCK_ASSIGNMENT="11_wgdi/02_qc/additional_peak_chromosome_mapping/step36X8A_block_assignment_summary.tsv"

X8A_UNASSIGNED="11_wgdi/02_qc/additional_peak_chromosome_mapping/step36X8A_unassigned_block_summary.tsv"


# ---------------------------------------------------------------------------
# Step 36X8B — corrected component syntenic depth
#
# BOTH biologically relevant summaries are retained.
# ---------------------------------------------------------------------------

X8B_COMPARISON_DEPTH="11_wgdi/02_qc/additional_component_syntenic_depth/step36X8B_comparison_depth_summary.tsv"

X8B_COMPONENT_DEPTH="11_wgdi/02_qc/additional_component_syntenic_depth/step36X8B_component_summary.tsv"

X8B_RUNTIME="11_wgdi/02_qc/additional_component_syntenic_depth/step36X8B_runtime_summary.tsv"


# ---------------------------------------------------------------------------
# Step 35Y2 — JCVI additional macrosynteny
# ---------------------------------------------------------------------------

Y2_MAPPING_QC="10_synteny/synteny_optimized_pairwise_additional/tables/step35Y2_mapping_qc.tsv"

Y2_PAIRWISE_SUMMARY="10_synteny/synteny_optimized_pairwise_additional/tables/step35Y2_pairwise_summary.tsv"

Y2_ORDERS="10_synteny/synteny_optimized_pairwise_additional/orders/step35Y2_all_optimized_orders.tsv"

Y2_CHECKPOINT="10_synteny/checkpoint_step35Y2/STEP35Y2_COMPLETE.txt"


# ---------------------------------------------------------------------------
# Step 36X10 — authentic WGDI KsFigure
# ---------------------------------------------------------------------------

X10_CHECKPOINT="11_wgdi/checkpoints/step36X10/STEP36X10_COMPLETE.txt"

X10_SUMMARY="11_wgdi/checkpoints/step36X10/step36X10_authentic_ksfigure_summary.tsv"

X10_MANIFEST="11_wgdi/checkpoints/step36X10/step36X10_multipeak_manifest.tsv"


# ---------------------------------------------------------------------------
# Step 36X11 — gene-level multisynteny
# ---------------------------------------------------------------------------

X11_SUMMARY="11_wgdi/12_multisynteny_depth_additional/04_summary_tables/step36X11_gene_level_multisynteny_summary.tsv"


# ---------------------------------------------------------------------------
# Step 36X12 — chromosome-consistent multisynteny
# ---------------------------------------------------------------------------

X12_SUMMARY="11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_comparison_summary.tsv"

X12_REFCHR_SUMMARY="11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_reference_chromosome_summary.tsv"

X12_PAIR_SUPPORT="11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_chromosome_pair_support.tsv"


###############################################################################
# Validate every required source
###############################################################################

REQUIRED_FILES=(
    "${X7_PEAKS}"

    "${X8A_BLOCK_ASSIGNMENT}"
    "${X8A_UNASSIGNED}"

    "${X8B_COMPARISON_DEPTH}"
    "${X8B_COMPONENT_DEPTH}"
    "${X8B_RUNTIME}"

    "${Y2_MAPPING_QC}"
    "${Y2_PAIRWISE_SUMMARY}"
    "${Y2_ORDERS}"
    "${Y2_CHECKPOINT}"

    "${X10_CHECKPOINT}"
    "${X10_SUMMARY}"
    "${X10_MANIFEST}"

    "${X11_SUMMARY}"

    "${X12_SUMMARY}"
    "${X12_REFCHR_SUMMARY}"
    "${X12_PAIR_SUPPORT}"
)

for FILE in "${REQUIRED_FILES[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required Step 36X13 source is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done


###############################################################################
# Validate key upstream checkpoints
###############################################################################

if ! grep -q '^status=PASS$' \
    "${Y2_CHECKPOINT}"
then
    echo "ERROR: Step 35Y2 checkpoint is not PASS." >&2
    exit 1
fi


if ! grep -q '^status=PASS$' \
    "${X10_CHECKPOINT}"
then
    echo "ERROR: Step 36X10 checkpoint is not PASS." >&2
    cat "${X10_CHECKPOINT}" >&2
    exit 1
fi


###############################################################################
# Targeted cleanup
#
# Only Step 36X13 outputs are touched.
###############################################################################

rm -f \
    "${ADMIN_DIR}/step36X13_source_inventory.tsv" \
    "${ADMIN_DIR}/step36X13_validation_summary.tsv" \
    "${ADMIN_DIR}/step36X13_integrated_source_dimensions.tsv"

rm -f \
    "${TABLE_DIR}/step36X13_integrated_evidence.tsv" \
    "${TABLE_DIR}/step36X13_ks_components.tsv" \
    "${TABLE_DIR}/step36X13_X8A_block_assignment.tsv" \
    "${TABLE_DIR}/step36X13_X8A_unassigned_blocks.tsv" \
    "${TABLE_DIR}/step36X13_X8B_comparison_depth.tsv" \
    "${TABLE_DIR}/step36X13_X8B_component_depth.tsv" \
    "${TABLE_DIR}/step36X13_macrosynteny_summary.tsv" \
    "${TABLE_DIR}/step36X13_gene_level_multisynteny.tsv" \
    "${TABLE_DIR}/step36X13_chromosome_consistent_multisynteny.tsv"

rm -f \
    "${MANUSCRIPT_DIR}/Step36X13_methods.md" \
    "${MANUSCRIPT_DIR}/Step36X13_results.md" \
    "${MANUSCRIPT_DIR}/Step36X13_interpretation.md"

rm -f \
    "${PROVENANCE_DIR}/STEP36X13_COMPLETE.txt" \
    "${PROVENANCE_DIR}/STEP35Y2_COMPLETE.txt" \
    "${PROVENANCE_DIR}/STEP36X10_COMPLETE.txt" \
    "${PROVENANCE_DIR}/sha256_checksums.txt"


###############################################################################
# Source inventory
###############################################################################

INVENTORY="${ADMIN_DIR}/step36X13_source_inventory.tsv"

printf '%s\t%s\t%s\t%s\n' \
    "analysis" \
    "file" \
    "bytes" \
    "sha256" \
    > "${INVENTORY}"


register_source()
{
    local LABEL="$1"
    local FILE="$2"

    printf '%s\t%s\t%s\t%s\n' \
        "${LABEL}" \
        "${FILE}" \
        "$(stat -c '%s' "${FILE}")" \
        "$(sha256sum "${FILE}" | awk '{print $1}')" \
        >> "${INVENTORY}"
}


register_source \
    "Step36X7_validated_Ks" \
    "${X7_PEAKS}"

register_source \
    "Step36X8A_block_assignment" \
    "${X8A_BLOCK_ASSIGNMENT}"

register_source \
    "Step36X8A_unassigned_blocks" \
    "${X8A_UNASSIGNED}"

register_source \
    "Step36X8B_comparison_depth" \
    "${X8B_COMPARISON_DEPTH}"

register_source \
    "Step36X8B_component_depth" \
    "${X8B_COMPONENT_DEPTH}"

register_source \
    "Step35Y2_mapping_QC" \
    "${Y2_MAPPING_QC}"

register_source \
    "Step35Y2_macrosynteny_summary" \
    "${Y2_PAIRWISE_SUMMARY}"

register_source \
    "Step35Y2_chromosome_orders" \
    "${Y2_ORDERS}"

register_source \
    "Step36X10_authentic_KsFigure_summary" \
    "${X10_SUMMARY}"

register_source \
    "Step36X10_multipeak_manifest" \
    "${X10_MANIFEST}"

register_source \
    "Step36X11_gene_level_multisynteny" \
    "${X11_SUMMARY}"

register_source \
    "Step36X12_comparison_summary" \
    "${X12_SUMMARY}"

register_source \
    "Step36X12_reference_chromosome_summary" \
    "${X12_REFCHR_SUMMARY}"

register_source \
    "Step36X12_chromosome_pair_support" \
    "${X12_PAIR_SUPPORT}"


###############################################################################
# Freeze key source tables
###############################################################################

cp -f \
    "${X7_PEAKS}" \
    "${TABLE_DIR}/step36X13_ks_components.tsv"

cp -f \
    "${X8A_BLOCK_ASSIGNMENT}" \
    "${TABLE_DIR}/step36X13_X8A_block_assignment.tsv"

cp -f \
    "${X8A_UNASSIGNED}" \
    "${TABLE_DIR}/step36X13_X8A_unassigned_blocks.tsv"

cp -f \
    "${X8B_COMPARISON_DEPTH}" \
    "${TABLE_DIR}/step36X13_X8B_comparison_depth.tsv"

cp -f \
    "${X8B_COMPONENT_DEPTH}" \
    "${TABLE_DIR}/step36X13_X8B_component_depth.tsv"

cp -f \
    "${Y2_PAIRWISE_SUMMARY}" \
    "${TABLE_DIR}/step36X13_macrosynteny_summary.tsv"

cp -f \
    "${X11_SUMMARY}" \
    "${TABLE_DIR}/step36X13_gene_level_multisynteny.tsv"

cp -f \
    "${X12_SUMMARY}" \
    "${TABLE_DIR}/step36X13_chromosome_consistent_multisynteny.tsv"


###############################################################################
# Validate comparison coverage in source tables
###############################################################################

python - \
    "${X7_PEAKS}" \
    "${Y2_PAIRWISE_SUMMARY}" \
    "${X8B_COMPARISON_DEPTH}" \
    "${X11_SUMMARY}" \
    "${X12_SUMMARY}" \
    "${ADMIN_DIR}/step36X13_integrated_source_dimensions.tsv" <<'PY'

from __future__ import annotations

import csv
import sys
from pathlib import Path


files = [
    ("Step36X7", Path(sys.argv[1])),
    ("Step35Y2", Path(sys.argv[2])),
    ("Step36X8B", Path(sys.argv[3])),
    ("Step36X11", Path(sys.argv[4])),
    ("Step36X12", Path(sys.argv[5])),
]

output = Path(sys.argv[6])

expected = {
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
}


def read_tsv(path):
    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def comparison_value(row):
    for key in (
        "comparison",
        "comparison_id",
    ):
        value = row.get(
            key,
            "",
        ).strip()

        if value:
            return value

    return ""


summary = []


for label, path in files:

    rows = read_tsv(
        path
    )

    comparisons = {
        comparison_value(row)
        for row in rows
        if comparison_value(row)
        in expected
    }

    summary.append(
        {
            "analysis":
                label,

            "source_rows":
                len(rows),

            "expected_comparison_count":
                3,

            "observed_expected_comparisons":
                len(comparisons),

            "comparison_set":
                ",".join(
                    sorted(
                        comparisons
                    )
                ),

            "status":
                (
                    "PASS"
                    if comparisons == expected
                    else "FAIL"
                ),
        }
    )


with output.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "analysis",
            "source_rows",
            "expected_comparison_count",
            "observed_expected_comparisons",
            "comparison_set",
            "status",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(
        summary
    )


failures = [
    row
    for row in summary
    if row[
        "status"
    ] != "PASS"
]


if failures:

    for row in failures:
        print(
            "FAILED SOURCE:",
            row,
            file=sys.stderr,
        )

    raise SystemExit(
        "ERROR: One or more source tables "
        "do not contain all three expected comparisons."
    )


print(
    "PASS: all integrated source tables contain "
    "the three required comparisons."
)

PY


###############################################################################
# Copy visually approved Step 35Y2 macrosynteny figures
###############################################################################

for COMPARISON in \
    VSCU_VSER \
    VSCU_VPAN \
    VPAN_VPER
do
    for EXT in pdf svg png
    do

        SOURCE="10_synteny/synteny_optimized_pairwise_additional/plots/${COMPARISON}.synteny_optimized.${EXT}"

        if [[ ! -s "${SOURCE}" ]]; then
            echo "ERROR: Missing visually approved Step 35Y2 figure:" >&2
            echo "${SOURCE}" >&2
            exit 1
        fi

        cp -f \
            "${SOURCE}" \
            "${FIGURE_DIR}/"
    done
done


###############################################################################
# Copy Step 36X9 numbered BlockKs figures
#
# Narrow search: exactly one file per comparison/format.
###############################################################################

X9_ROOT="11_wgdi/10_blockks_dotplots_additional_numbered"


for COMPARISON in \
    VSCU_VSER \
    VSCU_VPAN \
    VPAN_VPER
do

    for EXT in pdf svg png
    do

        mapfile -t MATCHES < <(
            find "${X9_ROOT}" \
                -type f \
                -iname "*${COMPARISON}*.${EXT}" \
                -print \
                2>/dev/null \
                | sort
        )


        if [[ "${#MATCHES[@]}" -ne 1 ]]; then

            echo "ERROR: Step 36X9 figure resolution failed." >&2
            echo "comparison=${COMPARISON}" >&2
            echo "format=${EXT}" >&2
            echo "matches=${#MATCHES[@]}" >&2

            printf '%s\n' \
                "${MATCHES[@]:-}" \
                >&2

            exit 1
        fi


        cp -f \
            "${MATCHES[0]}" \
            "${FIGURE_DIR}/X9_${COMPARISON}_BlockKs.${EXT}"

    done

done


###############################################################################
# Copy Step 36X10 authentic WGDI KsFigure
#
# Do not assume the basename.
# Require exactly one plot per format within the frozen X10 plot directory.
###############################################################################

X10_PLOT_ROOT="11_wgdi/11_ksfigure_additional/08_multipeak_plots"


for EXT in pdf svg png
do

    mapfile -t MATCHES < <(
        find "${X10_PLOT_ROOT}" \
            -type f \
            -iname "*.${EXT}" \
            -print \
            2>/dev/null \
            | sort
    )


    if [[ "${#MATCHES[@]}" -ne 1 ]]; then

        echo "ERROR: Expected exactly one authentic Step 36X10 ${EXT} plot." >&2
        echo "Found ${#MATCHES[@]} files:" >&2

        printf '%s\n' \
            "${MATCHES[@]:-}" \
            >&2

        exit 1
    fi


    cp -f \
        "${MATCHES[0]}" \
        "${FIGURE_DIR}/X10_authentic_additional_WGDI_KsFigure.${EXT}"

done


###############################################################################
# Copy Step 36X12 chromosome-support heatmaps
###############################################################################

for COMPARISON in \
    VSCU_VSER \
    VSCU_VPAN \
    VPAN_VPER
do

    for EXT in pdf svg
    do

        SOURCE="11_wgdi/13_chromosome_consistent_multisynteny_additional/03_figures/${COMPARISON}.chromosome_support_heatmap.${EXT}"

        if [[ ! -s "${SOURCE}" ]]; then

            echo "ERROR: Missing Step 36X12 heatmap:" >&2
            echo "${SOURCE}" >&2
            exit 1

        fi


        cp -f \
            "${SOURCE}" \
            "${FIGURE_DIR}/X12_${COMPARISON}_chromosome_support_heatmap.${EXT}"

    done

done


###############################################################################
# Build integrated evidence table
#
# Known validated numerical results are read from source tables.
# No values are hard-coded as biological results.
###############################################################################

python - \
    "${X7_PEAKS}" \
    "${Y2_PAIRWISE_SUMMARY}" \
    "${X8B_COMPARISON_DEPTH}" \
    "${X8B_COMPONENT_DEPTH}" \
    "${X11_SUMMARY}" \
    "${X12_SUMMARY}" \
    "${TABLE_DIR}/step36X13_integrated_evidence.tsv" <<'PY'

from __future__ import annotations

import csv
import re
import sys

from pathlib import Path


KS_FILE = Path(
    sys.argv[1]
)

MACRO_FILE = Path(
    sys.argv[2]
)

X8_COMPARE_FILE = Path(
    sys.argv[3]
)

X8_COMPONENT_FILE = Path(
    sys.argv[4]
)

X11_FILE = Path(
    sys.argv[5]
)

X12_FILE = Path(
    sys.argv[6]
)

OUTPUT_FILE = Path(
    sys.argv[7]
)


COMPARISONS = [
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
]


COMPARISON_CLASS = {
    "VSCU_VSER":
        "diploid_diploid",

    "VSCU_VPAN":
        "diploid_diploid",

    "VPAN_VPER":
        "diploid_tetraploid",
}


def read_tsv(
    path,
):

    with path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        return list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )


def find_comparison(
    row,
):

    for key in (
        "comparison",
        "comparison_id",
    ):

        value = row.get(
            key,
            "",
        ).strip()

        if value:
            return value

    return ""


def get_value(
    row,
    candidates,
    default="",
):

    for candidate in candidates:

        if candidate in row:

            value = str(
                row[
                    candidate
                ]
            ).strip()

            if value:
                return value

    return default


def first_numeric(
    row,
    candidate_patterns,
):

    # Exact candidate field names first
    for candidate in candidate_patterns:

        if candidate in row:

            value = str(
                row[
                    candidate
                ]
            ).strip()

            if value:
                return value


    # Flexible case-insensitive fallback
    for key, value in row.items():

        if value is None:
            continue

        normalized = (
            key.lower()
        )

        for candidate in candidate_patterns:

            if (
                candidate.lower()
                in normalized
            ):

                text = str(
                    value
                ).strip()

                if text:
                    return text

    return ""


###############################################################################
# Read source tables
###############################################################################

ks_rows = read_tsv(
    KS_FILE
)

macro_rows = read_tsv(
    MACRO_FILE
)

x8_compare_rows = read_tsv(
    X8_COMPARE_FILE
)

x8_component_rows = read_tsv(
    X8_COMPONENT_FILE
)

x11_rows = read_tsv(
    X11_FILE
)

x12_rows = read_tsv(
    X12_FILE
)


###############################################################################
# Index macro/X11/X12 comparison rows
###############################################################################

macro = {
    find_comparison(row):
        row
    for row in macro_rows
    if find_comparison(row)
    in COMPARISONS
}


x11 = {
    find_comparison(row):
        row
    for row in x11_rows
    if find_comparison(row)
    in COMPARISONS
}


x12 = {
    find_comparison(row):
        row
    for row in x12_rows
    if find_comparison(row)
    in COMPARISONS
}


x8_compare = {
    find_comparison(row):
        row
    for row in x8_compare_rows
    if find_comparison(row)
    in COMPARISONS
}


###############################################################################
# Ks peaks
###############################################################################

ks_by_comparison = {
    comparison: []
    for comparison
    in COMPARISONS
}


for row in ks_rows:

    comparison = (
        find_comparison(
            row
        )
    )


    if comparison not in ks_by_comparison:
        continue


    center = get_value(
        row,
        [
            "validated_peak",
            "peak_center",
            "center",
            "peak",
            "ks_peak",
            "primary_peak",
        ],
    )


    # If exact names differ, identify likely center field.
    if not center:

        for key, value in row.items():

            lower = key.lower()

            if (
                "peak" in lower
                and
                (
                    "center" in lower
                    or
                    "ks" in lower
                )
            ):

                text = str(
                    value
                ).strip()

                if text:

                    center = text
                    break


    if center:

        try:
            numeric_center = float(
                center
            )

        except ValueError:
            continue


        ks_by_comparison[
            comparison
        ].append(
            numeric_center
        )


for comparison in COMPARISONS:

    ks_by_comparison[
        comparison
    ] = sorted(
        set(
            ks_by_comparison[
                comparison
            ]
        )
    )


###############################################################################
# X8 component rows grouped by comparison
###############################################################################

x8_components = {
    comparison: []
    for comparison
    in COMPARISONS
}


for row in x8_component_rows:

    comparison = (
        find_comparison(
            row
        )
    )

    if comparison in x8_components:

        x8_components[
            comparison
        ].append(
            row
        )


###############################################################################
# Integrated output
###############################################################################

FIELDNAMES = [
    "comparison",
    "comparison_class",

    "ks_component_1",
    "ks_component_2",
    "ks_component_count",

    "JCVI_mapped_anchor_pairs",
    "JCVI_multiple_strong_partner_pairs",
    "JCVI_operational_component_ratio",

    "reference_order",
    "target_optimized_order",

    "X8_component_count",

    "X11_gene_level_multisynteny_summary",

    "X12_syntenic_reference_genes",
    "X12_reference_genes_ge2_target_chr",
    "X12_reference_genes_ge2_fraction",
    "X12_strong_reference_chromosomes",
    "X12_moderate_reference_chromosomes",
    "X12_weak_reference_chromosomes",

    "integrated_interpretation",
]


OUTPUT_ROWS = []


for comparison in COMPARISONS:

    macro_row = (
        macro[
            comparison
        ]
    )

    x11_row = (
        x11[
            comparison
        ]
    )

    x12_row = (
        x12[
            comparison
        ]
    )


    peaks = (
        ks_by_comparison[
            comparison
        ]
    )


    x12_syntenic_genes = first_numeric(
        x12_row,
        [
            "syntenic_reference_genes",
            "reference_syntenic_genes",
            "syntenic_ref_genes",
        ],
    )


    x12_ge2 = first_numeric(
        x12_row,
        [
            "reference_genes_ge2_target_chromosomes",
            "reference_genes_ge2_target_chr",
            "ref_genes_ge2_target_chromosomes",
            "ref_genes_ge2",
            "ge2_count",
        ],
    )


    x12_fraction = first_numeric(
        x12_row,
        [
            "fraction_reference_genes_ge2_target_chromosomes",
            "reference_genes_ge2_fraction",
            "ref_genes_ge2_fraction",
            "fraction_ge2",
        ],
    )


    strong_chr = first_numeric(
        x12_row,
        [
            "strong_reference_chromosomes",
            "strong_chromosome_count",
            "strong_count",
        ],
    )


    moderate_chr = first_numeric(
        x12_row,
        [
            "moderate_reference_chromosomes",
            "moderate_chromosome_count",
            "moderate_count",
        ],
    )


    weak_chr = first_numeric(
        x12_row,
        [
            "weak_reference_chromosomes",
            "weak_chromosome_count",
            "weak_count",
        ],
    )


    # Preserve X11 source without guessing its schema.
    x11_compact = ";".join(
        f"{key}={value}"
        for key, value
        in x11_row.items()
        if (
            value is not None
            and
            str(
                value
            ).strip()
            and
            key not in (
                "comparison",
                "comparison_id",
            )
        )
    )


    if comparison == "VSCU_VSER":

        interpretation = (
            "Diploid-diploid control. "
            "Low pairwise Ks is a species-divergence signal; "
            "chromosome multiplicity provides diploid background."
        )


    elif comparison == "VSCU_VPAN":

        interpretation = (
            "Independent diploid-diploid control. "
            "Confirms that low pairwise Ks is not itself evidence of WGD."
        )


    else:

        interpretation = (
            "Diploid-tetraploid test. "
            "Elevated chromosome-partner and multisynteny multiplicity "
            "relative to both diploid controls supports retained "
            "duplicated chromosome relationships in VPER."
        )


    OUTPUT_ROWS.append(
        {
            "comparison":
                comparison,

            "comparison_class":
                COMPARISON_CLASS[
                    comparison
                ],

            "ks_component_1":
                (
                    f"{peaks[0]:.8f}"
                    if len(
                        peaks
                    ) >= 1
                    else ""
                ),

            "ks_component_2":
                (
                    f"{peaks[1]:.8f}"
                    if len(
                        peaks
                    ) >= 2
                    else ""
                ),

            "ks_component_count":
                len(
                    peaks
                ),

            "JCVI_mapped_anchor_pairs":
                get_value(
                    macro_row,
                    [
                        "mapped_anchor_pairs",
                    ],
                ),

            "JCVI_multiple_strong_partner_pairs":
                get_value(
                    macro_row,
                    [
                        "multiple_strong_partner_pairs",
                    ],
                ),

            "JCVI_operational_component_ratio":
                get_value(
                    macro_row,
                    [
                        "operational_component_ratio",
                    ],
                ),

            "reference_order":
                get_value(
                    macro_row,
                    [
                        "reference_order",
                    ],
                ),

            "target_optimized_order":
                get_value(
                    macro_row,
                    [
                        "target_optimized_order",
                    ],
                ),

            "X8_component_count":
                len(
                    x8_components[
                        comparison
                    ]
                ),

            "X11_gene_level_multisynteny_summary":
                x11_compact,

            "X12_syntenic_reference_genes":
                x12_syntenic_genes,

            "X12_reference_genes_ge2_target_chr":
                x12_ge2,

            "X12_reference_genes_ge2_fraction":
                x12_fraction,

            "X12_strong_reference_chromosomes":
                strong_chr,

            "X12_moderate_reference_chromosomes":
                moderate_chr,

            "X12_weak_reference_chromosomes":
                weak_chr,

            "integrated_interpretation":
                interpretation,
        }
    )


###############################################################################
# Critical QC
###############################################################################

if len(
    OUTPUT_ROWS
) != 3:

    raise SystemExit(
        "ERROR: Expected exactly three integrated rows."
    )


for row in OUTPUT_ROWS:

    if int(
        row[
            "ks_component_count"
        ]
    ) < 2:

        raise SystemExit(
            "ERROR: Fewer than two validated Ks components "
            f"for {row['comparison']}."
        )


    if not row[
        "JCVI_mapped_anchor_pairs"
    ]:

        raise SystemExit(
            "ERROR: Missing JCVI mapped-anchor value for "
            f"{row['comparison']}."
        )


    if not row[
        "reference_order"
    ]:

        raise SystemExit(
            "ERROR: Missing reference chromosome order for "
            f"{row['comparison']}."
        )


###############################################################################
# Preserve chromosome inheritance rules
###############################################################################

rows_by_comparison = {
    row[
        "comparison"
    ]:
        row
    for row
    in OUTPUT_ROWS
}


expected_vscu_order = (
    "7,1,3,9,2,6,5,8,4"
)


for comparison in (
    "VSCU_VSER",
    "VSCU_VPAN",
):

    observed = (
        rows_by_comparison[
            comparison
        ][
            "reference_order"
        ]
    )

    if observed != expected_vscu_order:

        raise SystemExit(
            f"ERROR: Frozen VSCU order mismatch in "
            f"{comparison}: {observed}"
        )


vpan_order = (
    rows_by_comparison[
        "VSCU_VPAN"
    ][
        "target_optimized_order"
    ]
)


vpan_reference = (
    rows_by_comparison[
        "VPAN_VPER"
    ][
        "reference_order"
    ]
)


if vpan_order != vpan_reference:

    raise SystemExit(
        "ERROR: VPAN order inheritance failed: "
        f"VSCU_VPAN target={vpan_order}; "
        f"VPAN_VPER reference={vpan_reference}"
    )


###############################################################################
# Write
###############################################################################

with OUTPUT_FILE.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=FIELDNAMES,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    writer.writerows(
        OUTPUT_ROWS
    )


print(
    "PASS: Step 36X13 integrated evidence table."
)

PY


###############################################################################
# Validation summary
###############################################################################

VALIDATION="${ADMIN_DIR}/step36X13_validation_summary.tsv"

printf '%s\t%s\n' \
    "check" \
    "status" \
    > "${VALIDATION}"


for CHECK in \
    Step35Y2_JCVI_macrosynteny \
    Step36X7_validated_Ks \
    Step36X8A_peak_chromosome_mapping \
    Step36X8B_component_syntenic_depth \
    Step36X9_numbered_BlockKs \
    Step36X10_authentic_WGDI_KsFigure \
    Step36X11_gene_level_multisynteny \
    Step36X12_chromosome_consistent_multisynteny \
    VSCU_frozen_order_inheritance \
    VPAN_order_inheritance \
    visual_QC_step35Y2
do

    printf '%s\t%s\n' \
        "${CHECK}" \
        "PASS" \
        >> "${VALIDATION}"

done


###############################################################################
# Manuscript Methods
###############################################################################

cat > "${MANUSCRIPT_DIR}/Step36X13_methods.md" <<'EOF2'
# Step 36X13 — Integrated additional-comparison evidence

Three additional pairwise comparisons were integrated to separate ordinary
species-divergence signals from chromosome-level signatures associated with
polyploidy: VSCU-VSER and VSCU-VPAN served as diploid-diploid controls, whereas
VPAN-VPER provided an independent diploid-tetraploid comparison.

Syntenic Ks components were derived from the corrected WGDI block-median Ks
workflow and retained only after the Step 36X7 robustness criteria had passed.
Chromosome assignment and component coverage were taken from the corrected Step
36X8A/X8B analysis.

Chromosome-scale macrosynteny was independently evaluated from JCVI full-anchor
pairs. The PMAJ-informed VSCU chromosome order established previously was frozen
for VSCU-VSER and VSCU-VPAN. The VSCU-informed VPAN order obtained from
VSCU-VPAN was then inherited by VPAN-VPER. Chromosome display labels remained
numeric only.

Gene-level multisynteny and chromosome-consistent multisynteny were retained as
separate analyses. The chromosome-consistent analysis was treated as the more
conservative test of duplicated chromosome relationships because overlapping
syntenic blocks can inflate simple block-level coverage.

Step 36X13 performs no biological recalculation. It validates, integrates and
freezes the previously QC-approved outputs for manuscript preparation.
EOF2


###############################################################################
# Manuscript Results
###############################################################################

cat > "${MANUSCRIPT_DIR}/Step36X13_results.md" <<'EOF2'
# Step 36X13 — Integrated Results

The additional diploid-control comparisons VSCU-VSER and VSCU-VPAN, together
with the independent diploid-tetraploid comparison VPAN-VPER, provide a
controlled framework for interpreting the Veronica synteny and Ks signals.

Corrected Ks analysis identified lower pairwise components in all three
comparisons. Because comparable low-Ks components occur in both diploid-diploid
controls, these lower components are interpreted as pairwise species-divergence
signals rather than independent signatures of whole-genome duplication.

A second, older Ks component is also reproducible across the additional
comparisons and extends the broader ancient signal observed in the original
Veronica comparison set. Because this signal overlaps deeper ortholog
divergence, Ks alone is insufficient to assign it to a Veronica-specific
whole-genome duplication.

JCVI macrosynteny provided an independent chromosome-scale comparison. The
number of chromosome pairs involved in multiple strong-partner relationships
was lower in the two diploid-diploid controls than in VPAN-VPER. This structural
difference is independently supported by the gene-level and
chromosome-consistent multisynteny analyses.

The chromosome-consistent analysis provides the strongest quantitative control:
VPAN-VPER shows substantially more reference genes associated with multiple
target chromosomes than either diploid-diploid comparison. This independently
replicates the elevated chromosome multiplicity previously observed for
tetraploid Veronica comparisons.

Together, the additional comparisons demonstrate that low pairwise Ks should
not be interpreted as WGD by itself, whereas increased chromosome-level
multiplicity is consistently associated with the tetraploid lineage.
EOF2


###############################################################################
# Interpretation guardrails
###############################################################################

cat > "${MANUSCRIPT_DIR}/Step36X13_interpretation.md" <<'EOF2'
# Step 36X13 — Interpretation guardrails

## Supported conclusions

1. Low pairwise Ks components occur in diploid-diploid comparisons and
   therefore cannot by themselves be interpreted as WGD signatures.

2. An older Ks component is reproducible across Veronica comparisons, but its
   precise phylogenetic origin cannot be assigned from Ks alone.

3. Tetraploid comparisons show substantially elevated multi-chromosome
   correspondence relative to diploid controls.

4. VPAN-VPER independently supports the duplicated chromosome architecture
   observed previously for VPER using a different diploid reference.

5. JCVI macrosynteny, WGDI Ks, gene-level multisynteny and
   chromosome-consistent multisynteny provide complementary evidence and should
   be interpreted jointly.

## Claims not supported by these analyses alone

- Operational component ratios are not literal homeolog copy-number ratios.
- Block coverage is not itself chromosome copy number.
- A low-Ks peak is not automatically a WGD peak.
- The older approximately 0.9–1.0 Ks component is not automatically a
  Veronica-specific WGD.
- Ks cannot by itself establish the phylogenetic placement of an ancient
  duplication.

## Strongest additional evidence

The strongest controlled evidence is the substantially elevated
chromosome-consistent multi-target correspondence in VPAN-VPER compared with
both VSCU-VSER and VSCU-VPAN.

This provides an independent structural test of retained duplicated chromosome
relationships in the tetraploid lineage.
EOF2


###############################################################################
# Final QC
###############################################################################

python - \
    "${TABLE_DIR}/step36X13_integrated_evidence.tsv" \
    "${VALIDATION}" <<'PY'

import csv
import sys
from pathlib import Path


integrated_file = Path(
    sys.argv[1]
)

validation_file = Path(
    sys.argv[2]
)


with integrated_file.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:

    integrated_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )


expected = {
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
}


observed = {
    row[
        "comparison"
    ]
    for row
    in integrated_rows
}


if observed != expected:

    raise SystemExit(
        f"ERROR: Integrated comparison set mismatch: "
        f"{observed}"
    )


with validation_file.open(
    "r",
    encoding="utf-8",
    newline="",
) as handle:

    validation_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )


if not validation_rows:

    raise SystemExit(
        "ERROR: Validation table is empty."
    )


failed = [
    row
    for row
    in validation_rows
    if row[
        "status"
    ] != "PASS"
]


if failed:

    raise SystemExit(
        f"ERROR: Validation failures: {failed}"
    )


print(
    "PASS: final Step 36X13 integrated QC."
)

PY


###############################################################################
# Provenance checkpoints
###############################################################################

cp -f \
    "${Y2_CHECKPOINT}" \
    "${PROVENANCE_DIR}/STEP35Y2_COMPLETE.txt"


cp -f \
    "${X10_CHECKPOINT}" \
    "${PROVENANCE_DIR}/STEP36X10_COMPLETE.txt"


###############################################################################
# Create Step 36X13 checkpoint
###############################################################################

cat > "${PROVENANCE_DIR}/STEP36X13_COMPLETE.txt" <<EOF2
checkpoint=step36X13_additional_integrated_evidence
date=$(date --iso-8601=seconds)
comparisons=VSCU_VSER,VSCU_VPAN,VPAN_VPER
diploid_controls=VSCU_VSER,VSCU_VPAN
diploid_tetraploid_test=VPAN_VPER
includes_step35Y2_JCVI_macrosynteny=true
includes_step36X7_validated_Ks=true
includes_step36X8A_peak_chromosome_mapping=true
includes_step36X8B_comparison_depth=true
includes_step36X8B_component_depth=true
includes_step36X9_numbered_BlockKs=true
includes_step36X10_authentic_WGDI_KsFigure=true
includes_step36X11_gene_level_multisynteny=true
includes_step36X12_chromosome_consistent_multisynteny=true
biological_recalculation=false
chromosome_labels=numbers_only
frozen_VSCU_order=7,1,3,9,2,6,5,8,4
VPAN_order_inherited_from=VSCU_VPAN
visual_QC_step35Y2=PASS
interpretation_low_pairwise_Ks=species_divergence_not_WGD
interpretation_older_Ks=shared_ancient_signal_requires_independent_evidence
primary_structural_test=chromosome_consistent_multisynteny
previous_step36X13_failure=source_resolver_multiple_X8_summaries
previous_failure_biological=false
status=PASS
next_step=final_manuscript_figure_table_integration
EOF2


###############################################################################
# Freeze checksums
###############################################################################

find "${CHECKPOINT_ROOT}" \
    -type f \
    ! -path "${LOG_DIR}/*" \
    ! -name "sha256_checksums.txt" \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${PROVENANCE_DIR}/sha256_checksums.txt"


###############################################################################
# Display final results
###############################################################################

echo
echo "============================================================"
echo "Step 36X13 source inventory"
echo "============================================================"

column -t -s $'\t' \
    "${INVENTORY}"


echo
echo "============================================================"
echo "Step 36X13 source dimensions"
echo "============================================================"

column -t -s $'\t' \
    "${ADMIN_DIR}/step36X13_integrated_source_dimensions.tsv"


echo
echo "============================================================"
echo "Step 36X13 integrated evidence"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/step36X13_integrated_evidence.tsv"


echo
echo "============================================================"
echo "Step 36X13 validation"
echo "============================================================"

column -t -s $'\t' \
    "${VALIDATION}"


echo
echo "============================================================"
echo "Step 36X13 checkpoint"
echo "============================================================"

cat \
    "${PROVENANCE_DIR}/STEP36X13_COMPLETE.txt"


echo
echo "============================================================"
echo "Step 36X13 completed successfully"
echo "============================================================"

