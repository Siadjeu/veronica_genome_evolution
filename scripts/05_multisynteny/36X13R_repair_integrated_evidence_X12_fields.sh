#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=36X13R_fix
#SBATCH --output=11_wgdi/checkpoints/step36X13_integrated_evidence/logs/step36X13R_%j.out
#SBATCH --error=11_wgdi/checkpoints/step36X13_integrated_evidence/logs/step36X13R_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

ROOT="11_wgdi/checkpoints/step36X13_integrated_evidence"

TABLE="${ROOT}/01_tables/step36X13_integrated_evidence.tsv"

X12="11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_comparison_summary.tsv"

PROV="${ROOT}/04_provenance"
LOG_DIR="${ROOT}/logs"

###############################################################################
# Validate
###############################################################################

for FILE in \
    "${TABLE}" \
    "${X12}" \
    "${PROV}/STEP36X13_COMPLETE.txt"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required file missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# Backup current integrated table before repair
###############################################################################

BACKUP="${ROOT}/01_tables/step36X13_integrated_evidence.pre36X13R.tsv"

cp -f \
    "${TABLE}" \
    "${BACKUP}"

###############################################################################
# Rebuild only the X12-derived columns
###############################################################################

python - \
    "${TABLE}" \
    "${X12}" <<'PY'

from __future__ import annotations

import csv
import sys
from pathlib import Path


integrated_file = Path(sys.argv[1])
x12_file = Path(sys.argv[2])


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


integrated_rows = read_tsv(
    integrated_file
)

x12_rows = read_tsv(
    x12_file
)


x12_by_comparison = {
    row["comparison"]: row
    for row in x12_rows
}


expected = {
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
}


observed = {
    row["comparison"]
    for row in integrated_rows
}


if observed != expected:
    raise SystemExit(
        f"ERROR: Unexpected integrated comparison set: {observed}"
    )


if set(x12_by_comparison).intersection(expected) != expected:
    raise SystemExit(
        "ERROR: X12 table does not contain all three comparisons."
    )


###############################################################################
# Populate exact X12 fields
###############################################################################

for row in integrated_rows:

    comparison = row["comparison"]

    source = x12_by_comparison[
        comparison
    ]

    row[
        "X12_syntenic_reference_genes"
    ] = source[
        "reference_syntenic_genes"
    ]

    row[
        "X12_reference_genes_ge2_target_chr"
    ] = source[
        "genes_ge2_target_chromosomes"
    ]

    row[
        "X12_reference_genes_ge2_fraction"
    ] = source[
        "fraction_ge2_target_chromosomes"
    ]

    row[
        "X12_strong_reference_chromosomes"
    ] = source[
        "strong_1to2_candidate_chromosomes"
    ]

    row[
        "X12_moderate_reference_chromosomes"
    ] = source[
        "moderate_balanced_chromosomes"
    ]

    row[
        "X12_weak_reference_chromosomes"
    ] = source[
        "weak_second_partner_chromosomes"
    ]


###############################################################################
# Validate exact expected values
###############################################################################

expected_values = {
    "VSCU_VSER": {
        "X12_syntenic_reference_genes": "1152",
        "X12_reference_genes_ge2_target_chr": "3",
        "X12_reference_genes_ge2_fraction": "0.00260417",
        "X12_strong_reference_chromosomes": "1",
        "X12_moderate_reference_chromosomes": "4",
        "X12_weak_reference_chromosomes": "4",
    },
    "VSCU_VPAN": {
        "X12_syntenic_reference_genes": "1140",
        "X12_reference_genes_ge2_target_chr": "3",
        "X12_reference_genes_ge2_fraction": "0.00263158",
        "X12_strong_reference_chromosomes": "0",
        "X12_moderate_reference_chromosomes": "3",
        "X12_weak_reference_chromosomes": "6",
    },
    "VPAN_VPER": {
        "X12_syntenic_reference_genes": "1466",
        "X12_reference_genes_ge2_target_chr": "240",
        "X12_reference_genes_ge2_fraction": "0.16371078",
        "X12_strong_reference_chromosomes": "0",
        "X12_moderate_reference_chromosomes": "9",
        "X12_weak_reference_chromosomes": "0",
    },
}


for row in integrated_rows:

    comparison = row["comparison"]

    for field, expected_value in (
        expected_values[
            comparison
        ].items()
    ):

        observed_value = row[
            field
        ]

        if observed_value != expected_value:

            raise SystemExit(
                f"ERROR: {comparison} {field}: "
                f"expected {expected_value}, "
                f"observed {observed_value}"
            )


###############################################################################
# Rewrite integrated table
###############################################################################

fieldnames = list(
    integrated_rows[0].keys()
)


with integrated_file.open(
    "w",
    encoding="utf-8",
    newline="",
) as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()

    writer.writerows(
        integrated_rows
    )


print(
    "PASS: Step 36X13 X12-derived fields repaired."
)

PY

###############################################################################
# Add repair provenance
###############################################################################

cat > "${PROV}/STEP36X13R_COMPLETE.txt" <<EOF2
checkpoint=step36X13R_repair_integrated_X12_fields
date=$(date --iso-8601=seconds)
biological_recalculation=false
source_X12=${X12}
repaired_fields=X12_syntenic_reference_genes,X12_reference_genes_ge2_target_chr,X12_reference_genes_ge2_fraction,X12_strong_reference_chromosomes,X12_moderate_reference_chromosomes,X12_weak_reference_chromosomes
backup=${BACKUP}
status=PASS
next_step=final_manuscript_figure_table_integration
EOF2

###############################################################################
# Refresh checksums
###############################################################################

rm -f \
    "${PROV}/sha256_checksums.txt"

find "${ROOT}" \
    -type f \
    ! -path "${LOG_DIR}/*" \
    ! -name "sha256_checksums.txt" \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${PROV}/sha256_checksums.txt"

###############################################################################
# Display corrected integrated table
###############################################################################

echo
echo "============================================================"
echo "Step 36X13R corrected integrated evidence"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE}"

echo
echo "============================================================"
echo "Step 36X13R checkpoint"
echo "============================================================"

cat \
    "${PROV}/STEP36X13R_COMPLETE.txt"

echo
echo "============================================================"
echo "Step 36X13R completed successfully"
echo "============================================================"
