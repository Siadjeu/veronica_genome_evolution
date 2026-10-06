#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36H2
#SBATCH --output=11_wgdi/logs/step36H2_%j.out
#SBATCH --error=11_wgdi/logs/step36H2_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

H1_PAIR_DIR="${WGDI_ROOT}/12_multisynteny_depth/01_unique_gene_pairs"

MAP_FILE="${WGDI_ROOT}/10_blockks_dotplots_numbered/01_chromosome_maps/ALL_SPECIES.chromosome_number_map.tsv"

OUT_ROOT="${WGDI_ROOT}/12_multisynteny_depth/05_chromosome_consistent_1to2"
MATRIX_DIR="${OUT_ROOT}/01_matrices"
RANK_DIR="${OUT_ROOT}/02_ranked_partners"
GENE_DIR="${OUT_ROOT}/03_top2_gene_support"
TABLE_DIR="${OUT_ROOT}/04_summary_tables"
FIGURE_DIR="${OUT_ROOT}/05_figures"

QC_DIR="${WGDI_ROOT}/02_qc/multisynteny_depth"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36H2"

RANK_SUMMARY="${TABLE_DIR}/step36H2_ranked_chromosome_partners.tsv"
TOP2_SUMMARY="${TABLE_DIR}/step36H2_top2_multisynteny_summary.tsv"
COMPARISON_SUMMARY="${TABLE_DIR}/step36H2_comparison_summary.tsv"
GENE_SUPPORT_SUMMARY="${TABLE_DIR}/step36H2_gene_support_summary.tsv"
QC_SUMMARY="${QC_DIR}/step36H2_qc_summary.tsv"

mkdir -p \
    "${MATRIX_DIR}" \
    "${RANK_DIR}" \
    "${GENE_DIR}" \
    "${TABLE_DIR}" \
    "${FIGURE_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

###############################################################################
# REQUIRE STEP 36H1
###############################################################################

H1_CHECKPOINT="${WGDI_ROOT}/checkpoints/step36H1/STEP36H1_COMPLETE.txt"

if [[ ! -s "${H1_CHECKPOINT}" ]]; then
    echo "ERROR: Missing Step 36H1 checkpoint." >&2
    exit 1
fi

if ! grep -q '^status=PASS$' "${H1_CHECKPOINT}"; then
    echo "ERROR: Step 36H1 checkpoint is not PASS." >&2
    exit 1
fi

if [[ ! -s "${MAP_FILE}" ]]; then
    echo "ERROR: Missing chromosome map: ${MAP_FILE}" >&2
    exit 1
fi

###############################################################################
# REQUIRED H1 INPUTS
###############################################################################

COMPARISONS=(
    PMAJ_VSCU
    VSCU_VANA
    VSCU_VPER
)

for COMPARISON in "${COMPARISONS[@]}"
do
    INPUT="${H1_PAIR_DIR}/${COMPARISON}.unique_gene_pairs.tsv"

    if [[ ! -s "${INPUT}" ]]; then
        echo "ERROR: Missing H1 unique-pair table: ${INPUT}" >&2
        exit 1
    fi
done

###############################################################################
# TARGETED CLEANUP OF STEP 36H2 ONLY
###############################################################################

rm -f \
    "${RANK_SUMMARY}" \
    "${TOP2_SUMMARY}" \
    "${COMPARISON_SUMMARY}" \
    "${GENE_SUPPORT_SUMMARY}" \
    "${QC_SUMMARY}" \
    "${CHECKPOINT_DIR}/STEP36H2_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

find "${MATRIX_DIR}" \
    -maxdepth 1 \
    -type f \
    -delete

find "${RANK_DIR}" \
    -maxdepth 1 \
    -type f \
    -delete

find "${GENE_DIR}" \
    -maxdepth 1 \
    -type f \
    -delete

find "${FIGURE_DIR}" \
    -maxdepth 1 \
    -type f \
    \( -name '*.pdf' -o -name '*.svg' -o -name '*.png' \) \
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

export MPLBACKEND=Agg
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-2}"

python - <<'PY'
import pandas
import numpy
import matplotlib

print("pandas:", pandas.__version__)
print("numpy:", numpy.__version__)
print("matplotlib:", matplotlib.__version__)
PY

###############################################################################
# ANALYSIS
###############################################################################

python - \
    "${H1_PAIR_DIR}" \
    "${MAP_FILE}" \
    "${MATRIX_DIR}" \
    "${RANK_DIR}" \
    "${GENE_DIR}" \
    "${FIGURE_DIR}" \
    "${RANK_SUMMARY}" \
    "${TOP2_SUMMARY}" \
    "${COMPARISON_SUMMARY}" \
    "${GENE_SUPPORT_SUMMARY}" \
    "${QC_SUMMARY}" <<'PY'
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


pair_dir = Path(sys.argv[1])
map_file = Path(sys.argv[2])
matrix_dir = Path(sys.argv[3])
rank_dir = Path(sys.argv[4])
gene_dir = Path(sys.argv[5])
figure_dir = Path(sys.argv[6])

rank_summary_path = Path(sys.argv[7])
top2_summary_path = Path(sys.argv[8])
comparison_summary_path = Path(sys.argv[9])
gene_support_summary_path = Path(sys.argv[10])
qc_summary_path = Path(sys.argv[11])


###############################################################################
# ANALYSIS DESIGN
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
# LOAD CHROMOSOME ACCESSION -> ChrN MAP
###############################################################################

map_df = pd.read_csv(
    map_file,
    sep="\t",
    dtype=str,
)

required_map_columns = {
    "species",
    "display_order",
    "original_chromosome_id",
    "display_chromosome_id",
}

missing = required_map_columns - set(map_df.columns)

if missing:
    raise SystemExit(
        f"ERROR: chromosome map missing columns: {sorted(missing)}"
    )

map_df["display_order"] = map_df["display_order"].astype(int)

display_map = {}
chromosome_order = {}

for species, group in map_df.groupby("species"):

    group = group.sort_values("display_order")

    display_map[species] = dict(
        zip(
            group["original_chromosome_id"],
            group["display_chromosome_id"],
        )
    )

    chromosome_order[species] = group[
        "original_chromosome_id"
    ].tolist()


###############################################################################
# GLOBAL OUTPUT ACCUMULATORS
###############################################################################

rank_rows = []
top2_rows = []
comparison_rows = []
gene_support_rows = []
qc_rows = []


###############################################################################
# PROCESS EACH COMPARISON
###############################################################################

for spec in comparisons:

    comparison = spec["comparison"]
    reference = spec["reference"]
    target = spec["target"]
    role = spec["role"]

    print()
    print("=" * 80)
    print(comparison)
    print("=" * 80)

    pair_path = (
        pair_dir
        / f"{comparison}.unique_gene_pairs.tsv"
    )

    pairs = pd.read_csv(
        pair_path,
        sep="\t",
        dtype=str,
    )

    required_columns = {
        "reference_chr",
        "target_chr",
        "reference_gene",
        "target_gene",
    }

    missing = required_columns - set(pairs.columns)

    if missing:
        raise SystemExit(
            f"ERROR: {comparison} missing columns: {sorted(missing)}"
        )

    ###########################################################################
    # VERIFY H1 PAIRS ARE UNIQUE
    ###########################################################################

    duplicate_gene_pairs = int(
        pairs.duplicated(
            subset=[
                "reference_gene",
                "target_gene",
            ]
        ).sum()
    )

    if duplicate_gene_pairs != 0:
        raise SystemExit(
            f"ERROR: {comparison} contains "
            f"{duplicate_gene_pairs} duplicate gene pairs."
        )

    ###########################################################################
    # VERIFY CHROMOSOME IDENTIFIERS
    ###########################################################################

    ref_unknown = sorted(
        set(pairs["reference_chr"])
        - set(display_map[reference])
    )

    target_unknown = sorted(
        set(pairs["target_chr"])
        - set(display_map[target])
    )

    if ref_unknown:
        raise SystemExit(
            f"ERROR: {comparison}: unmapped reference chromosomes: "
            f"{ref_unknown}"
        )

    if target_unknown:
        raise SystemExit(
            f"ERROR: {comparison}: unmapped target chromosomes: "
            f"{target_unknown}"
        )

    ###########################################################################
    # CHROMOSOME-PAIR SUPPORT
    #
    # Primary support = number of UNIQUE REFERENCE GENES connecting
    # one reference chromosome to one target chromosome.
    ###########################################################################

    support = (
        pairs.groupby(
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

    ref_order = [
        chrom
        for chrom in chromosome_order[reference]
        if chrom in set(pairs["reference_chr"])
    ]

    target_order = [
        chrom
        for chrom in chromosome_order[target]
        if chrom in set(pairs["target_chr"])
    ]

    ###########################################################################
    # SUPPORT MATRIX
    ###########################################################################

    matrix = (
        support
        .pivot(
            index="reference_chr",
            columns="target_chr",
            values="unique_reference_genes",
        )
        .reindex(
            index=ref_order,
            columns=target_order,
        )
        .fillna(0)
        .astype(int)
    )

    numbered_matrix = matrix.copy()

    numbered_matrix.index = [
        display_map[reference][chrom]
        for chrom in matrix.index
    ]

    numbered_matrix.columns = [
        display_map[target][chrom]
        for chrom in matrix.columns
    ]

    numbered_matrix.to_csv(
        matrix_dir
        / f"{comparison}.chromosome_support_matrix.tsv",
        sep="\t",
    )

    comparison_top2 = []

    ###########################################################################
    # RANK TARGET CHROMOSOMES FOR EVERY REFERENCE CHROMOSOME
    ###########################################################################

    for ref_chr in ref_order:

        ref_pairs = pairs[
            pairs["reference_chr"] == ref_chr
        ].copy()

        total_ref_syntenic_genes = int(
            ref_pairs["reference_gene"].nunique()
        )

        target_stats = []

        for target_chr, group in ref_pairs.groupby(
            "target_chr"
        ):

            target_stats.append(
                {
                    "target_chr": target_chr,
                    "unique_reference_genes": int(
                        group["reference_gene"].nunique()
                    ),
                    "unique_target_genes": int(
                        group["target_gene"].nunique()
                    ),
                    "unique_gene_pairs": int(
                        len(group)
                    ),
                }
            )

        target_stats.sort(
            key=lambda row: (
                -row["unique_reference_genes"],
                -row["unique_gene_pairs"],
                chromosome_order[target].index(
                    row["target_chr"]
                ),
            )
        )

        if not target_stats:
            continue

        #######################################################################
        # WRITE ALL RANKED PARTNERS
        #######################################################################

        for rank, row in enumerate(
            target_stats,
            start=1,
        ):

            rank_rows.append(
                {
                    "comparison": comparison,
                    "role": role,
                    "reference_species": reference,
                    "target_species": target,
                    "reference_chr_original": ref_chr,
                    "reference_chr": display_map[
                        reference
                    ][ref_chr],
                    "target_rank": rank,
                    "target_chr_original": row["target_chr"],
                    "target_chr": display_map[
                        target
                    ][row["target_chr"]],
                    "unique_reference_genes": row[
                        "unique_reference_genes"
                    ],
                    "unique_target_genes": row[
                        "unique_target_genes"
                    ],
                    "unique_gene_pairs": row[
                        "unique_gene_pairs"
                    ],
                    "reference_chr_syntenic_genes_total": (
                        total_ref_syntenic_genes
                    ),
                    "support_fraction": (
                        row["unique_reference_genes"]
                        / total_ref_syntenic_genes
                        if total_ref_syntenic_genes
                        else 0
                    ),
                    "status": "PASS",
                }
            )

        #######################################################################
        # TOP 1 AND TOP 2 CHROMOSOME PARTNERS
        #######################################################################

        top1 = target_stats[0]

        top2 = (
            target_stats[1]
            if len(target_stats) >= 2
            else None
        )

        top1_chr = top1["target_chr"]

        top2_chr = (
            top2["target_chr"]
            if top2 is not None
            else None
        )

        top1_genes = set(
            ref_pairs.loc[
                ref_pairs["target_chr"] == top1_chr,
                "reference_gene",
            ]
        )

        top2_genes = (
            set(
                ref_pairs.loc[
                    ref_pairs["target_chr"] == top2_chr,
                    "reference_gene",
                ]
            )
            if top2_chr is not None
            else set()
        )

        union_genes = (
            top1_genes
            | top2_genes
        )

        dual_genes = (
            top1_genes
            & top2_genes
        )

        top1_fraction = (
            len(top1_genes)
            / total_ref_syntenic_genes
            if total_ref_syntenic_genes
            else 0
        )

        top2_fraction = (
            len(top2_genes)
            / total_ref_syntenic_genes
            if total_ref_syntenic_genes
            else 0
        )

        second_to_first_ratio = (
            len(top2_genes)
            / len(top1_genes)
            if top1_genes
            else 0
        )

        top12_union_fraction = (
            len(union_genes)
            / total_ref_syntenic_genes
            if total_ref_syntenic_genes
            else 0
        )

        dual_fraction = (
            len(dual_genes)
            / total_ref_syntenic_genes
            if total_ref_syntenic_genes
            else 0
        )

        #######################################################################
        # DESCRIPTIVE CLASSIFICATION
        #
        # Strong 1:2 candidate requires:
        #   second partner >=50% of first-partner support
        #   top two together cover >=70% of syntenic reference genes
        #######################################################################

        if top2 is None:
            pattern = "single_dominant_partner"

        elif (
            second_to_first_ratio >= 0.50
            and top12_union_fraction >= 0.70
        ):
            pattern = "strong_1to2_candidate"

        elif second_to_first_ratio >= 0.25:
            pattern = "moderate_secondary_partner"

        else:
            pattern = "weak_secondary_partner"

        row = {
            "comparison": comparison,
            "role": role,
            "reference_species": reference,
            "target_species": target,
            "reference_chr_original": ref_chr,
            "reference_chr": display_map[
                reference
            ][ref_chr],
            "reference_chr_syntenic_genes": total_ref_syntenic_genes,
            "top1_target_chr_original": top1_chr,
            "top1_target_chr": display_map[
                target
            ][top1_chr],
            "top1_supported_reference_genes": len(top1_genes),
            "top1_fraction": top1_fraction,
            "top2_target_chr_original": (
                top2_chr
                if top2_chr is not None
                else ""
            ),
            "top2_target_chr": (
                display_map[target][top2_chr]
                if top2_chr is not None
                else ""
            ),
            "top2_supported_reference_genes": len(top2_genes),
            "top2_fraction": top2_fraction,
            "second_to_first_support_ratio": second_to_first_ratio,
            "top1_top2_union_reference_genes": len(union_genes),
            "top1_top2_union_fraction": top12_union_fraction,
            "reference_genes_on_both_top_target_chr": len(dual_genes),
            "dual_top2_fraction": dual_fraction,
            "pattern_class": pattern,
            "status": "PASS",
        }

        top2_rows.append(row)
        comparison_top2.append(row)

        #######################################################################
        # GENE-LEVEL TOP-2 SUPPORT TABLE
        #######################################################################

        for gene in sorted(
            set(ref_pairs["reference_gene"])
        ):

            gene_targets = set(
                ref_pairs.loc[
                    ref_pairs["reference_gene"] == gene,
                    "target_chr",
                ]
            )

            on_top1 = (
                top1_chr in gene_targets
            )

            on_top2 = (
                top2_chr in gene_targets
                if top2_chr is not None
                else False
            )

            if on_top1 and on_top2:
                support_class = "both_top1_top2"
            elif on_top1:
                support_class = "top1_only"
            elif on_top2:
                support_class = "top2_only"
            else:
                support_class = "other_target_chr"

            gene_support_rows.append(
                {
                    "comparison": comparison,
                    "reference_species": reference,
                    "target_species": target,
                    "reference_chr_original": ref_chr,
                    "reference_chr": display_map[
                        reference
                    ][ref_chr],
                    "reference_gene": gene,
                    "top1_target_chr": display_map[
                        target
                    ][top1_chr],
                    "top2_target_chr": (
                        display_map[target][top2_chr]
                        if top2_chr is not None
                        else ""
                    ),
                    "gene_target_chromosome_count": len(
                        gene_targets
                    ),
                    "support_class": support_class,
                    "status": "PASS",
                }
            )

    ###########################################################################
    # WRITE COMPARISON-SPECIFIC TOP-2 TABLE
    ###########################################################################

    comparison_top2_df = pd.DataFrame(
        comparison_top2
    )

    if comparison_top2_df.empty:
        raise SystemExit(
            f"ERROR: No chromosome-level results for {comparison}."
        )

    comparison_top2_df.to_csv(
        rank_dir
        / f"{comparison}.top2_chromosome_partners.tsv",
        sep="\t",
        index=False,
    )

    ###########################################################################
    # COMPARISON-LEVEL SUMMARY
    ###########################################################################

    gene_target_chr_counts = (
        pairs.groupby(
            "reference_gene"
        )["target_chr"]
        .nunique()
    )

    n_reference_genes = int(
        len(gene_target_chr_counts)
    )

    genes_exactly_2_chr = int(
        (
            gene_target_chr_counts == 2
        ).sum()
    )

    genes_ge2_chr = int(
        (
            gene_target_chr_counts >= 2
        ).sum()
    )

    class_counts = (
        comparison_top2_df[
            "pattern_class"
        ]
        .value_counts()
        .to_dict()
    )

    n_ref_chr = int(
        len(comparison_top2_df)
    )

    n_strong = int(
        class_counts.get(
            "strong_1to2_candidate",
            0,
        )
    )

    comparison_rows.append(
        {
            "comparison": comparison,
            "role": role,
            "reference_species": reference,
            "target_species": target,
            "reference_chromosomes_tested": n_ref_chr,
            "reference_syntenic_genes": n_reference_genes,
            "reference_genes_on_exactly_2_target_chromosomes": (
                genes_exactly_2_chr
            ),
            "reference_genes_on_ge2_target_chromosomes": (
                genes_ge2_chr
            ),
            "fraction_reference_genes_on_ge2_target_chromosomes": (
                genes_ge2_chr / n_reference_genes
                if n_reference_genes
                else 0
            ),
            "strong_1to2_candidate_reference_chromosomes": (
                n_strong
            ),
            "moderate_secondary_partner_reference_chromosomes": int(
                class_counts.get(
                    "moderate_secondary_partner",
                    0,
                )
            ),
            "weak_secondary_partner_reference_chromosomes": int(
                class_counts.get(
                    "weak_secondary_partner",
                    0,
                )
            ),
            "single_dominant_partner_reference_chromosomes": int(
                class_counts.get(
                    "single_dominant_partner",
                    0,
                )
            ),
            "fraction_reference_chromosomes_strong_1to2": (
                n_strong / n_ref_chr
                if n_ref_chr
                else 0
            ),
            "median_second_to_first_support_ratio": float(
                comparison_top2_df[
                    "second_to_first_support_ratio"
                ].median()
            ),
            "median_top1_top2_union_fraction": float(
                comparison_top2_df[
                    "top1_top2_union_fraction"
                ].median()
            ),
            "median_dual_top2_fraction": float(
                comparison_top2_df[
                    "dual_top2_fraction"
                ].median()
            ),
            "status": "PASS",
        }
    )

    ###########################################################################
    # HEATMAP
    ###########################################################################

    fig_width = max(
        8,
        0.55 * numbered_matrix.shape[1] + 3,
    )

    fig_height = max(
        6,
        0.55 * numbered_matrix.shape[0] + 2.5,
    )

    fig, ax = plt.subplots(
        figsize=(
            fig_width,
            fig_height,
        )
    )

    image = ax.imshow(
        numbered_matrix.values,
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_xticks(
        np.arange(
            numbered_matrix.shape[1]
        )
    )

    ax.set_xticklabels(
        numbered_matrix.columns,
        rotation=90,
    )

    ax.set_yticks(
        np.arange(
            numbered_matrix.shape[0]
        )
    )

    ax.set_yticklabels(
        numbered_matrix.index
    )

    ax.set_xlabel(target)
    ax.set_ylabel(reference)

    ax.set_title(
        f"{reference} → {target}: chromosome-level syntenic support"
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
    )

    colorbar.set_label(
        "Unique syntenic reference genes"
    )

    max_value = (
        int(numbered_matrix.values.max())
        if numbered_matrix.size
        else 0
    )

    for i in range(
        numbered_matrix.shape[0]
    ):
        for j in range(
            numbered_matrix.shape[1]
        ):

            value = int(
                numbered_matrix.iloc[i, j]
            )

            if value == 0:
                continue

            ax.text(
                j,
                i,
                str(value),
                ha="center",
                va="center",
                fontsize=7,
            )

    ax.set_xticks(
        np.arange(
            -0.5,
            numbered_matrix.shape[1],
            1,
        ),
        minor=True,
    )

    ax.set_yticks(
        np.arange(
            -0.5,
            numbered_matrix.shape[0],
            1,
        ),
        minor=True,
    )

    ax.grid(
        which="minor",
        linewidth=0.5,
    )

    ax.tick_params(
        which="minor",
        bottom=False,
        left=False,
    )

    fig.tight_layout()

    for extension in [
        "pdf",
        "svg",
        "png",
    ]:

        output = (
            figure_dir
            / (
                f"{comparison}."
                f"chromosome_1to2_support_heatmap."
                f"{extension}"
            )
        )

        fig.savefig(
            output,
            dpi=300,
            bbox_inches="tight",
        )

    plt.close(fig)

    ###########################################################################
    # QC RECORD
    ###########################################################################

    qc_rows.append(
        {
            "comparison": comparison,
            "unique_gene_pairs": int(len(pairs)),
            "reference_genes": int(
                pairs[
                    "reference_gene"
                ].nunique()
            ),
            "reference_chromosomes": int(
                pairs[
                    "reference_chr"
                ].nunique()
            ),
            "target_chromosomes": int(
                pairs[
                    "target_chr"
                ].nunique()
            ),
            "chromosome_matrix_rows": int(
                numbered_matrix.shape[0]
            ),
            "chromosome_matrix_columns": int(
                numbered_matrix.shape[1]
            ),
            "top2_rows": int(
                len(comparison_top2_df)
            ),
            "duplicate_gene_pairs": duplicate_gene_pairs,
            "status": "PASS",
        }
    )


###############################################################################
# WRITE COMBINED OUTPUTS
###############################################################################

rank_df = pd.DataFrame(rank_rows)
top2_df = pd.DataFrame(top2_rows)
comparison_df = pd.DataFrame(comparison_rows)
gene_support_df = pd.DataFrame(gene_support_rows)
qc_df = pd.DataFrame(qc_rows)

rank_df.to_csv(
    rank_summary_path,
    sep="\t",
    index=False,
)

top2_df.to_csv(
    top2_summary_path,
    sep="\t",
    index=False,
)

comparison_df.to_csv(
    comparison_summary_path,
    sep="\t",
    index=False,
)

gene_support_df.to_csv(
    gene_support_summary_path,
    sep="\t",
    index=False,
)

qc_df.to_csv(
    qc_summary_path,
    sep="\t",
    index=False,
)


###############################################################################
# FINAL INTERNAL VALIDATION
###############################################################################

if len(comparison_df) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 comparisons; found {len(comparison_df)}."
    )

if len(qc_df) != 3:
    raise SystemExit(
        f"ERROR: Expected 3 QC rows; found {len(qc_df)}."
    )

if not (comparison_df["status"] == "PASS").all():
    raise SystemExit(
        "ERROR: Comparison summary contains non-PASS rows."
    )

if not (qc_df["status"] == "PASS").all():
    raise SystemExit(
        "ERROR: QC summary contains non-PASS rows."
    )

print()
print("=" * 80)
print("Step 36H2 analysis PASS")
print("=" * 80)
print()
print(
    comparison_df.to_string(
        index=False
    )
)
PY

###############################################################################
# VALIDATE OUTPUT TABLES
###############################################################################

for FILE in \
    "${RANK_SUMMARY}" \
    "${TOP2_SUMMARY}" \
    "${COMPARISON_SUMMARY}" \
    "${GENE_SUPPORT_SUMMARY}" \
    "${QC_SUMMARY}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing output table: ${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# VALIDATE FIGURES
###############################################################################

for COMPARISON in "${COMPARISONS[@]}"
do
    for EXT in pdf svg png
    do
        FIGURE="${FIGURE_DIR}/${COMPARISON}.chromosome_1to2_support_heatmap.${EXT}"

        if [[ ! -s "${FIGURE}" ]]; then
            echo "ERROR: Missing figure: ${FIGURE}" >&2
            exit 1
        fi
    done
done

PDF_COUNT="$(
    find "${FIGURE_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.pdf' \
        | wc -l
)"

SVG_COUNT="$(
    find "${FIGURE_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.svg' \
        | wc -l
)"

PNG_COUNT="$(
    find "${FIGURE_DIR}" \
        -maxdepth 1 \
        -type f \
        -name '*.png' \
        | wc -l
)"

if [[ "${PDF_COUNT}" -ne 3 ]]; then
    echo "ERROR: Expected 3 PDFs; found ${PDF_COUNT}." >&2
    exit 1
fi

if [[ "${SVG_COUNT}" -ne 3 ]]; then
    echo "ERROR: Expected 3 SVGs; found ${SVG_COUNT}." >&2
    exit 1
fi

if [[ "${PNG_COUNT}" -ne 3 ]]; then
    echo "ERROR: Expected 3 PNGs; found ${PNG_COUNT}." >&2
    exit 1
fi

###############################################################################
# DISPLAY RESULTS
###############################################################################

echo
echo "============================================================"
echo "STEP 36H2 COMPARISON SUMMARY"
echo "============================================================"

column -t -s $'\t' \
    "${COMPARISON_SUMMARY}"

echo
echo "============================================================"
echo "TOP-2 CHROMOSOME PARTNERS"
echo "============================================================"

column -t -s $'\t' \
    "${TOP2_SUMMARY}"

echo
echo "============================================================"
echo "QC"
echo "============================================================"

column -t -s $'\t' \
    "${QC_SUMMARY}"

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36H2_COMPLETE.txt" <<EOF2
checkpoint=step36H2_chromosome_consistent_1to2_multisynteny
date=$(date --iso-8601=seconds)
analysis=chromosome_consistent_multisyntenic_partner_depth
source=step36H1_unique_gene_pairs
comparisons=3
comparison_1=PMAJ_to_VSCU
comparison_2=VSCU_to_VANA
comparison_3=VSCU_to_VPER
unique_gene_pair_input=true
overlapping_block_duplicate_pairs_removed=true
primary_chromosome_support_metric=unique_reference_genes
target_chromosomes_ranked_per_reference_chromosome=true
top1_top2_analysis=true
second_to_first_support_ratio=true
dual_top2_gene_support=true
strong_1to2_candidate_definition=second_to_first_ratio_ge_0.50_and_top1_top2_union_fraction_ge_0.70
classification_is_descriptive_not_final_WGD_call=true
chromosome_display_labels=Chr1,Chr2,...
original_accession_ids_retained_in_tables=true
comparison_summary=11_wgdi/12_multisynteny_depth/05_chromosome_consistent_1to2/04_summary_tables/step36H2_comparison_summary.tsv
top2_summary=11_wgdi/12_multisynteny_depth/05_chromosome_consistent_1to2/04_summary_tables/step36H2_top2_multisynteny_summary.tsv
rank_summary=11_wgdi/12_multisynteny_depth/05_chromosome_consistent_1to2/04_summary_tables/step36H2_ranked_chromosome_partners.tsv
gene_support_summary=11_wgdi/12_multisynteny_depth/05_chromosome_consistent_1to2/04_summary_tables/step36H2_gene_support_summary.tsv
figures=11_wgdi/12_multisynteny_depth/05_chromosome_consistent_1to2/05_figures
status=PASS
next_step=step36H3_statistical_and_manuscript_integration
EOF2

cp -f \
    "${COMPARISON_SUMMARY}" \
    "${TOP2_SUMMARY}" \
    "${QC_SUMMARY}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "============================================================"
echo "CHECKPOINT"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36H2_COMPLETE.txt"

echo
echo "Step 36H2 completed successfully."
