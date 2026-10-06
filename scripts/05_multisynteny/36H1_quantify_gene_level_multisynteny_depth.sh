#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36H1
#SBATCH --output=11_wgdi/logs/step36H1_%j.out
#SBATCH --error=11_wgdi/logs/step36H1_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

BLOCKINFO_DIR="${WGDI_ROOT}/07_blockinfo/02_results"
LENS_DIR="${WGDI_ROOT}/01_inputs/lens"

OUT_ROOT="${WGDI_ROOT}/12_multisynteny_depth"
PAIR_DIR="${OUT_ROOT}/01_unique_gene_pairs"
GENE_DIR="${OUT_ROOT}/02_gene_partner_depth"
CHR_DIR="${OUT_ROOT}/03_chromosome_partner_depth"
TABLE_DIR="${OUT_ROOT}/04_summary_tables"

QC_DIR="${WGDI_ROOT}/02_qc/multisynteny_depth"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36H1"

PAIR_SUMMARY="${QC_DIR}/step36H1_unique_pair_summary.tsv"
DEPTH_SUMMARY="${TABLE_DIR}/step36H1_gene_partner_depth_summary.tsv"
CHR_SUMMARY="${TABLE_DIR}/step36H1_chromosome_partner_summary.tsv"
OVERALL="${TABLE_DIR}/step36H1_overall_summary.tsv"

mkdir -p \
    "${PAIR_DIR}" \
    "${GENE_DIR}" \
    "${CHR_DIR}" \
    "${TABLE_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

###############################################################################
# TARGETED CLEANUP
###############################################################################

rm -f \
    "${PAIR_SUMMARY}" \
    "${DEPTH_SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/STEP36H1_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

find "${PAIR_DIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.tsv' \
    -delete

find "${GENE_DIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.tsv' \
    -delete

find "${CHR_DIR}" \
    -maxdepth 1 \
    -type f \
    -name '*.tsv' \
    -delete

###############################################################################
# ACTIVATE WGDI ENVIRONMENT
###############################################################################

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

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

python - <<'PY'
import pandas
print("pandas:", pandas.__version__)
PY

###############################################################################
# REQUIRE PREVIOUS WGDI CHECKPOINT
###############################################################################

if [[ ! -s \
    "${WGDI_ROOT}/checkpoints/step36E/STEP36E_COMPLETE.txt" ]]
then
    echo "ERROR: Missing Step 36E checkpoint." >&2
    exit 1
fi

###############################################################################
# COMPARISONS
#
# Direction is biologically important:
#
# PMAJ -> VSCU = baseline/outgroup to diploid Veronica
# VSCU -> VANA = diploid reference to tetraploid
# VSCU -> VPER = diploid reference to tetraploid
###############################################################################

python - \
    "${PROJECT_ROOT}" \
    "${BLOCKINFO_DIR}" \
    "${LENS_DIR}" \
    "${PAIR_DIR}" \
    "${GENE_DIR}" \
    "${CHR_DIR}" \
    "${PAIR_SUMMARY}" \
    "${DEPTH_SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${OVERALL}" <<'PY'
from __future__ import annotations

import ast
import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


project_root = Path(sys.argv[1]).resolve()
blockinfo_dir = Path(sys.argv[2]).resolve()
lens_dir = Path(sys.argv[3]).resolve()
pair_dir = Path(sys.argv[4]).resolve()
gene_dir = Path(sys.argv[5]).resolve()
chr_dir = Path(sys.argv[6]).resolve()
pair_summary_path = Path(sys.argv[7])
depth_summary_path = Path(sys.argv[8])
chr_summary_path = Path(sys.argv[9])
overall_path = Path(sys.argv[10])


###############################################################################
# PRIMARY COMPARISONS
###############################################################################

comparisons = [
    {
        "comparison": "PMAJ_VSCU",
        "reference": "PMAJ",
        "target": "VSCU",
        "role": "baseline_outgroup_to_diploid",
    },
    {
        "comparison": "VSCU_VANA",
        "reference": "VSCU",
        "target": "VANA",
        "role": "diploid_to_tetraploid",
    },
    {
        "comparison": "VSCU_VPER",
        "reference": "VSCU",
        "target": "VPER",
        "role": "diploid_to_tetraploid",
    },
]


###############################################################################
# PARSING HELPERS
###############################################################################

def parse_gene_list(value) -> list[str]:
    """
    Robust parser for WGDI block1/block2 cells.

    Expected formats may include:
      ['gene1','gene2']
      gene1,gene2
      gene1 gene2
      gene1;gene2

    Empty/NA values return [].
    """

    if value is None:
        return []

    if isinstance(value, float) and math.isnan(value):
        return []

    text = str(value).strip()

    if not text or text.lower() in {
        "nan",
        "none",
        "na",
    }:
        return []

    # First try a Python-style list/tuple.
    if (
        (text.startswith("[") and text.endswith("]"))
        or
        (text.startswith("(") and text.endswith(")"))
    ):
        try:
            parsed = ast.literal_eval(text)

            if isinstance(parsed, (list, tuple)):
                return [
                    str(item).strip()
                    for item in parsed
                    if str(item).strip()
                ]
        except Exception:
            pass

    # Remove outer brackets if still present.
    text = text.strip("[]()")

    # WGDI outputs can use commas, whitespace or semicolons.
    items = re.split(
        r"[\s,;|]+",
        text,
    )

    cleaned = []

    for item in items:
        item = item.strip().strip("'\"")

        if item:
            cleaned.append(item)

    return cleaned


def load_lens(species: str):
    path = lens_dir / f"{species}.wgdi.lens"

    if not path.is_file():
        raise SystemExit(
            f"ERROR: Missing lens: {path}"
        )

    rows = []

    with path.open(
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            if not line.strip():
                continue

            fields = line.split()

            if len(fields) < 3:
                raise SystemExit(
                    f"ERROR: Malformed lens row "
                    f"{path}:{line_number}"
                )

            rows.append(
                {
                    "chromosome": fields[0],
                    "bp_length": int(float(fields[1])),
                    "ordered_gene_count": int(float(fields[2])),
                }
            )

    return rows


###############################################################################
# LOAD LENS COUNTS
###############################################################################

lens_cache = {}

for species in {
    row["reference"]
    for row in comparisons
}.union(
    {
        row["target"]
        for row in comparisons
    }
):
    lens_cache[species] = load_lens(species)


###############################################################################
# OUTPUT ACCUMULATORS
###############################################################################

pair_summary_rows = []
depth_summary_rows = []
chr_summary_rows = []
overall_rows = []


###############################################################################
# PROCESS COMPARISONS
###############################################################################

for item in comparisons:

    comparison = item["comparison"]
    reference = item["reference"]
    target = item["target"]
    role = item["role"]

    blockinfo_path = (
        blockinfo_dir
        / f"{comparison}.blockinfo.csv"
    )

    if not blockinfo_path.is_file():
        raise SystemExit(
            f"ERROR: Missing blockinfo: {blockinfo_path}"
        )

    print()
    print("=" * 72)
    print(comparison)
    print("=" * 72)

    df = pd.read_csv(blockinfo_path)

    required = {
        "chr1",
        "chr2",
        "block1",
        "block2",
    }

    missing = required - set(df.columns)

    if missing:
        raise SystemExit(
            f"ERROR: {comparison} missing columns: "
            f"{sorted(missing)}"
        )

    ###########################################################################
    # EXPAND BLOCK GENE LISTS
    ###########################################################################

    raw_pair_records = []
    mismatched_blocks = []

    for row_index, row in df.iterrows():

        genes1 = parse_gene_list(
            row["block1"]
        )

        genes2 = parse_gene_list(
            row["block2"]
        )

        if len(genes1) != len(genes2):

            mismatched_blocks.append(
                {
                    "row_index": row_index,
                    "chr1": str(row["chr1"]),
                    "chr2": str(row["chr2"]),
                    "n_block1": len(genes1),
                    "n_block2": len(genes2),
                }
            )

            continue

        if len(genes1) == 0:
            continue

        chr1 = str(row["chr1"])
        chr2 = str(row["chr2"])

        ks_median = (
            row["ks_median"]
            if "ks_median" in df.columns
            else ""
        )

        for gene1, gene2 in zip(
            genes1,
            genes2,
        ):
            raw_pair_records.append(
                {
                    "comparison": comparison,
                    "reference_species": reference,
                    "target_species": target,
                    "reference_chr": chr1,
                    "target_chr": chr2,
                    "reference_gene": gene1,
                    "target_gene": gene2,
                    "block_row": row_index,
                    "block_ks_median": ks_median,
                }
            )

    if mismatched_blocks:
        mismatch_path = (
            pair_dir
            / f"{comparison}.block_length_mismatches.tsv"
        )

        pd.DataFrame(
            mismatched_blocks
        ).to_csv(
            mismatch_path,
            sep="\t",
            index=False,
        )

        raise SystemExit(
            f"ERROR: {comparison}: "
            f"{len(mismatched_blocks)} blocks have "
            f"unequal block1/block2 gene counts. "
            f"See {mismatch_path}"
        )

    if not raw_pair_records:
        raise SystemExit(
            f"ERROR: No gene pairs recovered for "
            f"{comparison}."
        )

    raw_pairs = pd.DataFrame(
        raw_pair_records
    )

    raw_pair_count = len(raw_pairs)

    ###########################################################################
    # DEDUPLICATE EXACT GENE PAIRS
    #
    # Same gene pair may appear in overlapping WGDI blocks.
    # Partner multiplicity must count it only once.
    ###########################################################################

    unique_pairs = (
        raw_pairs
        .sort_values(
            [
                "reference_chr",
                "reference_gene",
                "target_chr",
                "target_gene",
                "block_row",
            ]
        )
        .drop_duplicates(
            subset=[
                "reference_species",
                "target_species",
                "reference_gene",
                "target_gene",
            ],
            keep="first",
        )
        .reset_index(drop=True)
    )

    unique_pair_count = len(unique_pairs)

    duplicate_pair_records = (
        raw_pair_count
        - unique_pair_count
    )

    pair_output = (
        pair_dir
        / f"{comparison}.unique_gene_pairs.tsv"
    )

    unique_pairs.to_csv(
        pair_output,
        sep="\t",
        index=False,
    )

    ###########################################################################
    # PARTNER MULTIPLICITY PER REFERENCE GENE
    ###########################################################################

    grouped = (
        unique_pairs
        .groupby(
            [
                "reference_chr",
                "reference_gene",
            ],
            sort=False,
        )
    )

    depth_rows = []

    for (
        ref_chr,
        ref_gene,
    ), group in grouped:

        partner_genes = sorted(
            set(
                group["target_gene"]
                .astype(str)
            )
        )

        partner_chromosomes = sorted(
            set(
                group["target_chr"]
                .astype(str)
            )
        )

        n_partners = len(
            partner_genes
        )

        n_partner_chromosomes = len(
            partner_chromosomes
        )

        if n_partners == 1:
            depth_class = "1"
        elif n_partners == 2:
            depth_class = "2"
        elif n_partners > 2:
            depth_class = ">2"
        else:
            depth_class = "0"

        if n_partner_chromosomes == 1:
            chromosome_depth_class = "1_chr"
        elif n_partner_chromosomes == 2:
            chromosome_depth_class = "2_chr"
        elif n_partner_chromosomes > 2:
            chromosome_depth_class = ">2_chr"
        else:
            chromosome_depth_class = "0_chr"

        two_distinct_chr = (
            n_partners >= 2
            and n_partner_chromosomes >= 2
        )

        depth_rows.append(
            {
                "comparison": comparison,
                "reference_species": reference,
                "target_species": target,
                "reference_chr": ref_chr,
                "reference_gene": ref_gene,
                "partner_gene_count": n_partners,
                "partner_chromosome_count": n_partner_chromosomes,
                "partner_depth_class": depth_class,
                "chromosome_depth_class": chromosome_depth_class,
                "two_or_more_partners_on_distinct_chromosomes": (
                    "YES"
                    if two_distinct_chr
                    else "NO"
                ),
                "target_chromosomes": ",".join(
                    partner_chromosomes
                ),
                "target_genes": ",".join(
                    partner_genes
                ),
            }
        )

    depth_df = pd.DataFrame(
        depth_rows
    )

    depth_output = (
        gene_dir
        / f"{comparison}.gene_partner_depth.tsv"
    )

    depth_df.to_csv(
        depth_output,
        sep="\t",
        index=False,
    )

    ###########################################################################
    # IMPORTANT DENOMINATOR
    #
    # We report:
    #
    # 1. among reference genes represented in syntenic blocks
    # 2. relative to total ordered genes in reference lens
    #
    # The second quantity implicitly includes genes with zero detected
    # partners.
    ###########################################################################

    total_reference_ordered_genes = sum(
        row["ordered_gene_count"]
        for row in lens_cache[
            reference
        ]
    )

    represented_reference_genes = (
        depth_df[
            "reference_gene"
        ].nunique()
    )

    counts = (
        depth_df[
            "partner_depth_class"
        ]
        .value_counts()
        .to_dict()
    )

    n1 = int(
        counts.get("1", 0)
    )

    n2 = int(
        counts.get("2", 0)
    )

    ngt2 = int(
        counts.get(">2", 0)
    )

    n0_inferred = max(
        total_reference_ordered_genes
        - represented_reference_genes,
        0,
    )

    two_chr_count = int(
        (
            depth_df[
                "two_or_more_partners_on_distinct_chromosomes"
            ]
            == "YES"
        ).sum()
    )

    ###########################################################################
    # SUMMARY TABLE
    ###########################################################################

    depth_summary_rows.append(
        {
            "comparison": comparison,
            "role": role,
            "reference_species": reference,
            "target_species": target,
            "reference_ordered_genes": total_reference_ordered_genes,
            "reference_genes_with_any_syntenic_partner": represented_reference_genes,
            "reference_genes_with_0_detected_partners_inferred": n0_inferred,
            "reference_genes_with_1_partner": n1,
            "reference_genes_with_2_partners": n2,
            "reference_genes_with_gt2_partners": ngt2,
            "reference_genes_with_ge2_partners_on_distinct_target_chromosomes": two_chr_count,
            "fraction_any_partner_of_reference": (
                represented_reference_genes
                / total_reference_ordered_genes
                if total_reference_ordered_genes
                else 0
            ),
            "fraction_1_partner_among_syntenic_reference_genes": (
                n1
                / represented_reference_genes
                if represented_reference_genes
                else 0
            ),
            "fraction_2_partners_among_syntenic_reference_genes": (
                n2
                / represented_reference_genes
                if represented_reference_genes
                else 0
            ),
            "fraction_gt2_partners_among_syntenic_reference_genes": (
                ngt2
                / represented_reference_genes
                if represented_reference_genes
                else 0
            ),
            "fraction_ge2_distinct_chr_among_syntenic_reference_genes": (
                two_chr_count
                / represented_reference_genes
                if represented_reference_genes
                else 0
            ),
            "status": "PASS",
        }
    )

    ###########################################################################
    # CHROMOSOME-PAIR DEPTH
    #
    # Count distinct target genes and reference genes by chr1 -> chr2.
    ###########################################################################

    chr_grouped = (
        unique_pairs
        .groupby(
            [
                "reference_chr",
                "target_chr",
            ]
        )
        .agg(
            unique_reference_genes=(
                "reference_gene",
                "nunique",
            ),
            unique_target_genes=(
                "target_gene",
                "nunique",
            ),
            unique_gene_pairs=(
                "target_gene",
                "size",
            ),
        )
        .reset_index()
    )

    chr_grouped.insert(
        0,
        "comparison",
        comparison,
    )

    chr_grouped.insert(
        1,
        "reference_species",
        reference,
    )

    chr_grouped.insert(
        2,
        "target_species",
        target,
    )

    chr_output = (
        chr_dir
        / f"{comparison}.chromosome_partner_depth.tsv"
    )

    chr_grouped.to_csv(
        chr_output,
        sep="\t",
        index=False,
    )

    for _, row in chr_grouped.iterrows():
        chr_summary_rows.append(
            {
                "comparison": comparison,
                "reference_species": reference,
                "target_species": target,
                "reference_chr": row[
                    "reference_chr"
                ],
                "target_chr": row[
                    "target_chr"
                ],
                "unique_reference_genes": int(
                    row[
                        "unique_reference_genes"
                    ]
                ),
                "unique_target_genes": int(
                    row[
                        "unique_target_genes"
                    ]
                ),
                "unique_gene_pairs": int(
                    row[
                        "unique_gene_pairs"
                    ]
                ),
                "status": "PASS",
            }
        )

    pair_summary_rows.append(
        {
            "comparison": comparison,
            "reference_species": reference,
            "target_species": target,
            "blockinfo_rows": len(df),
            "raw_gene_pair_records": raw_pair_count,
            "unique_gene_pairs": unique_pair_count,
            "duplicate_pair_records_removed": duplicate_pair_records,
            "duplicate_fraction": (
                duplicate_pair_records
                / raw_pair_count
                if raw_pair_count
                else 0
            ),
            "reference_genes_represented": represented_reference_genes,
            "status": "PASS",
        }
    )

    overall_rows.append(
        {
            "comparison": comparison,
            "role": role,
            "reference_species": reference,
            "target_species": target,
            "unique_gene_pairs": unique_pair_count,
            "reference_genes_with_any_partner": represented_reference_genes,
            "one_partner": n1,
            "two_partners": n2,
            "gt2_partners": ngt2,
            "ge2_partners_distinct_target_chr": two_chr_count,
            "status": "PASS",
        }
    )

    print(
        "raw pairs:",
        raw_pair_count,
    )

    print(
        "unique pairs:",
        unique_pair_count,
    )

    print(
        "reference genes:",
        represented_reference_genes,
    )

    print(
        "1 partner:",
        n1,
    )

    print(
        "2 partners:",
        n2,
    )

    print(
        ">2 partners:",
        ngt2,
    )

    print(
        ">=2 partners on distinct target chromosomes:",
        two_chr_count,
    )


###############################################################################
# WRITE COMBINED TABLES
###############################################################################

pd.DataFrame(
    pair_summary_rows
).to_csv(
    pair_summary_path,
    sep="\t",
    index=False,
)

pd.DataFrame(
    depth_summary_rows
).to_csv(
    depth_summary_path,
    sep="\t",
    index=False,
)

pd.DataFrame(
    chr_summary_rows
).to_csv(
    chr_summary_path,
    sep="\t",
    index=False,
)

pd.DataFrame(
    overall_rows
).to_csv(
    overall_path,
    sep="\t",
    index=False,
)


###############################################################################
# FINAL INTERNAL QC
###############################################################################

if len(pair_summary_rows) != 3:
    raise SystemExit(
        "ERROR: Expected 3 comparisons."
    )

if not all(
    row["status"] == "PASS"
    for row in pair_summary_rows
):
    raise SystemExit(
        "ERROR: Pair QC failed."
    )

print()
print("=" * 72)
print("Step 36H1 Python analysis PASS")
print("=" * 72)
PY

###############################################################################
# SHELL QC
###############################################################################

PAIR_PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            n++
        }
        END {
            print n + 0
        }
    ' "${PAIR_SUMMARY}"
)"

DEPTH_PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            n++
        }
        END {
            print n + 0
        }
    ' "${DEPTH_SUMMARY}"
)"

if [[ "${PAIR_PASS_COUNT}" -ne 3 ]]; then
    echo "ERROR: ${PAIR_PASS_COUNT}/3 pair summaries passed." >&2
    exit 1
fi

if [[ "${DEPTH_PASS_COUNT}" -ne 3 ]]; then
    echo "ERROR: ${DEPTH_PASS_COUNT}/3 depth summaries passed." >&2
    exit 1
fi

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36H1_COMPLETE.txt" <<EOF2
checkpoint=step36H1_gene_level_multisynteny_depth
date=$(date --iso-8601=seconds)
analysis=reference_gene_syntenic_partner_multiplicity
source_blockinfo=11_wgdi/07_blockinfo/02_results
comparisons=3
comparison_1=PMAJ_VSCU
comparison_2=VSCU_VANA
comparison_3=VSCU_VPER
reference_direction_1=PMAJ_to_VSCU
reference_direction_2=VSCU_to_VANA
reference_direction_3=VSCU_to_VPER
gene_pair_deduplication=true
overlapping_block_duplicate_pairs_counted_once=true
partner_depth_classes=0,1,2,>2
chromosome_consistency_test=true
zero_partner_count=inferred_from_reference_lens_ordered_gene_count
primary_question=does_diploid_reference_show_enrichment_for_two_chromosome_consistent_partners_in_tetraploids
pair_summary=11_wgdi/02_qc/multisynteny_depth/step36H1_unique_pair_summary.tsv
depth_summary=11_wgdi/12_multisynteny_depth/04_summary_tables/step36H1_gene_partner_depth_summary.tsv
chromosome_summary=11_wgdi/12_multisynteny_depth/04_summary_tables/step36H1_chromosome_partner_summary.tsv
overall_summary=11_wgdi/12_multisynteny_depth/04_summary_tables/step36H1_overall_summary.tsv
status=PASS
next_step=step36H2_chromosome_consistent_1to2_multisynteny
EOF2

cp -f \
    "${PAIR_SUMMARY}" \
    "${DEPTH_SUMMARY}" \
    "${OVERALL}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# DISPLAY RESULTS
###############################################################################

echo
echo "============================================================"
echo "UNIQUE GENE PAIR QC"
echo "============================================================"

column -t -s $'\t' \
    "${PAIR_SUMMARY}"

echo
echo "============================================================"
echo "GENE-LEVEL PARTNER DEPTH"
echo "============================================================"

column -t -s $'\t' \
    "${DEPTH_SUMMARY}"

echo
echo "============================================================"
echo "COMPACT OVERALL SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${OVERALL}"

echo
echo "============================================================"
echo "CHECKPOINT"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36H1_COMPLETE.txt"

echo
echo "Step 36H1 completed successfully."
