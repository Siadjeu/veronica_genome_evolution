#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36X12
#SBATCH --output=11_wgdi/logs/step36X12_%j.out
#SBATCH --error=11_wgdi/logs/step36X12_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

X11_ROOT="${WGDI_ROOT}/12_multisynteny_depth_additional"
PAIR_DIR="${X11_ROOT}/01_unique_gene_pairs"

OUT_ROOT="${WGDI_ROOT}/13_chromosome_consistent_multisynteny_additional"
TABLE_DIR="${OUT_ROOT}/01_tables"
MATRIX_DIR="${OUT_ROOT}/02_matrices"
FIGURE_DIR="${OUT_ROOT}/03_figures"

SUMMARY="${TABLE_DIR}/step36X12_comparison_summary.tsv"
CHR_SUMMARY="${TABLE_DIR}/step36X12_reference_chromosome_summary.tsv"
MATRIX_SUMMARY="${TABLE_DIR}/step36X12_chromosome_pair_support.tsv"

QC_DIR="${WGDI_ROOT}/02_qc/additional_chromosome_consistent_multisynteny"
QC_TABLE="${QC_DIR}/step36X12_qc.tsv"

CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X12"
LENS_DIR="${WGDI_ROOT}/01_inputs/lens"
LOG_DIR="${WGDI_ROOT}/logs"

mkdir -p \
    "${TABLE_DIR}" \
    "${MATRIX_DIR}" \
    "${FIGURE_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

rm -f \
    "${SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${MATRIX_SUMMARY}" \
    "${QC_TABLE}" \
    "${CHECKPOINT_DIR}/STEP36X12_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

for comparison in VSCU_VSER VSCU_VPAN VPAN_VPER
do
    rm -f \
        "${MATRIX_DIR}/${comparison}.chromosome_support_matrix.tsv" \
        "${FIGURE_DIR}/${comparison}.chromosome_support_heatmap.pdf" \
        "${FIGURE_DIR}/${comparison}.chromosome_support_heatmap.svg"
done

python - \
    "${PROJECT_ROOT}" \
    "${PAIR_DIR}" \
    "${LENS_DIR}" \
    "${TABLE_DIR}" \
    "${MATRIX_DIR}" \
    "${FIGURE_DIR}" \
    "${SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${MATRIX_SUMMARY}" \
    "${QC_TABLE}" <<'PY'
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

project_root = Path(sys.argv[1])
pair_dir = Path(sys.argv[2])
lens_dir = Path(sys.argv[3])
table_dir = Path(sys.argv[4])
matrix_dir = Path(sys.argv[5])
figure_dir = Path(sys.argv[6])
summary_path = Path(sys.argv[7])
chr_summary_path = Path(sys.argv[8])
matrix_summary_path = Path(sys.argv[9])
qc_path = Path(sys.argv[10])

COMPARISONS = [
    {
        "comparison": "VSCU_VSER",
        "reference": "VSCU",
        "target": "VSER",
        "comparison_class": "diploid_diploid",
    },
    {
        "comparison": "VSCU_VPAN",
        "reference": "VSCU",
        "target": "VPAN",
        "comparison_class": "diploid_diploid",
    },
    {
        "comparison": "VPAN_VPER",
        "reference": "VPAN",
        "target": "VPER",
        "comparison_class": "diploid_tetraploid",
    },
]

STRONG_SECOND_FIRST_RATIO = 0.50
STRONG_TOP12_UNION_FRACTION = 0.70


def read_lens(path):
    order = []
    counts = {}

    with path.open(
        encoding="utf-8-sig"
    ) as handle:

        for line_number, line in enumerate(
            handle,
            start=1,
        ):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                raise SystemExit(
                    f"ERROR: malformed lens row: "
                    f"{path}:{line_number}"
                )

            chromosome = fields[0]
            count = int(
                round(float(fields[2]))
            )

            order.append(chromosome)
            counts[chromosome] = count

    if not order:
        raise SystemExit(
            f"ERROR: empty lens: {path}"
        )

    return order, counts


def write_tsv(path, rows):
    if not rows:
        raise SystemExit(
            f"ERROR: no rows generated for {path}"
        )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


def display_map(chromosomes):
    return {
        chromosome: f"Chr{index}"
        for index, chromosome in enumerate(
            chromosomes,
            start=1,
        )
    }


all_chr_rows = []
all_matrix_rows = []
comparison_rows = []
qc_rows = []

for cfg in COMPARISONS:

    comparison = cfg["comparison"]
    reference = cfg["reference"]
    target = cfg["target"]

    pair_path = (
        pair_dir
        / f"{comparison}.unique_gene_pairs.tsv"
    )

    ref_lens_path = (
        lens_dir
        / f"{reference}.wgdi.lens"
    )

    target_lens_path = (
        lens_dir
        / f"{target}.wgdi.lens"
    )

    for required in [
        pair_path,
        ref_lens_path,
        target_lens_path,
    ]:
        if (
            not required.is_file()
            or required.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: missing X12 input: {required}"
            )

    ref_order, ref_counts = read_lens(
        ref_lens_path
    )

    target_order, target_counts = read_lens(
        target_lens_path
    )

    ref_display = display_map(
        ref_order
    )

    target_display = display_map(
        target_order
    )

    with pair_path.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:

        pairs = list(
            csv.DictReader(
                handle,
                delimiter="\t",
            )
        )

    if not pairs:
        raise SystemExit(
            f"ERROR: empty X11 pair table "
            f"for {comparison}"
        )

    ref_gene_to_targets = defaultdict(set)
    ref_gene_to_target_chrs = defaultdict(set)
    ref_chr_to_genes = defaultdict(set)
    ref_chr_target_chr_genes = defaultdict(set)

    for row in pairs:

        ref_gene = row["reference_gene"]
        ref_chr = row["reference_chromosome"]
        target_gene = row["target_gene"]
        target_chr = row["target_chromosome"]

        if ref_chr not in ref_counts:
            raise SystemExit(
                f"ERROR: {ref_chr} not in "
                f"{reference} lens"
            )

        if target_chr not in target_counts:
            raise SystemExit(
                f"ERROR: {target_chr} not in "
                f"{target} lens"
            )

        ref_gene_to_targets[
            ref_gene
        ].add(target_gene)

        ref_gene_to_target_chrs[
            ref_gene
        ].add(target_chr)

        ref_chr_to_genes[
            ref_chr
        ].add(ref_gene)

        ref_chr_target_chr_genes[
            (ref_chr, target_chr)
        ].add(ref_gene)

    # Matrix support rows
    comparison_matrix_rows = []

    for ref_chr in ref_order:

        for target_chr in target_order:

            genes = ref_chr_target_chr_genes.get(
                (ref_chr, target_chr),
                set(),
            )

            support = len(genes)

            comparison_matrix_rows.append(
                {
                    "comparison": comparison,
                    "reference_species": reference,
                    "target_species": target,
                    "reference_chromosome": ref_chr,
                    "reference_display": (
                        ref_display[ref_chr]
                    ),
                    "target_chromosome": target_chr,
                    "target_display": (
                        target_display[target_chr]
                    ),
                    "supporting_reference_genes": str(
                        support
                    ),
                }
            )

            all_matrix_rows.append(
                comparison_matrix_rows[-1]
            )

    matrix_path = (
        matrix_dir
        / f"{comparison}.chromosome_support_matrix.tsv"
    )

    write_tsv(
        matrix_path,
        comparison_matrix_rows,
    )

    strong_count = 0
    moderate_count = 0
    weak_count = 0
    single_count = 0

    chromosome_rows = []

    for ref_chr in ref_order:

        syntenic_genes = ref_chr_to_genes.get(
            ref_chr,
            set(),
        )

        if not syntenic_genes:
            continue

        supports = []

        for target_chr in target_order:

            support = len(
                ref_chr_target_chr_genes.get(
                    (ref_chr, target_chr),
                    set(),
                )
            )

            if support > 0:
                supports.append(
                    (
                        target_chr,
                        support,
                    )
                )

        supports.sort(
            key=lambda item: (
                -item[1],
                target_order.index(
                    item[0]
                ),
            )
        )

        top1_chr = supports[0][0]
        top1_support = supports[0][1]

        if len(supports) >= 2:
            top2_chr = supports[1][0]
            top2_support = supports[1][1]
        else:
            top2_chr = "NA"
            top2_support = 0

        second_first_ratio = (
            top2_support / top1_support
            if top1_support > 0
            else 0.0
        )

        top1_genes = (
            ref_chr_target_chr_genes[
                (ref_chr, top1_chr)
            ]
        )

        if top2_chr != "NA":
            top2_genes = (
                ref_chr_target_chr_genes[
                    (ref_chr, top2_chr)
                ]
            )
        else:
            top2_genes = set()

        union_genes = (
            top1_genes | top2_genes
        )

        dual_top2_genes = (
            top1_genes & top2_genes
        )

        union_fraction = (
            len(union_genes)
            / len(syntenic_genes)
        )

        dual_top2_fraction = (
            len(dual_top2_genes)
            / len(syntenic_genes)
        )

        genes_exactly_2_target_chr = sum(
            len(
                ref_gene_to_target_chrs[
                    gene
                ]
            ) == 2
            for gene in syntenic_genes
        )

        genes_ge2_target_chr = sum(
            len(
                ref_gene_to_target_chrs[
                    gene
                ]
            ) >= 2
            for gene in syntenic_genes
        )

        if top2_support == 0:
            classification = (
                "single_target"
            )
            single_count += 1

        elif (
            second_first_ratio
            >= STRONG_SECOND_FIRST_RATIO
            and union_fraction
            >= STRONG_TOP12_UNION_FRACTION
        ):
            classification = (
                "strong_1to2_candidate"
            )
            strong_count += 1

        elif (
            second_first_ratio
            >= STRONG_SECOND_FIRST_RATIO
        ):
            classification = (
                "moderate_balanced_two_partner"
            )
            moderate_count += 1

        else:
            classification = (
                "weak_second_partner"
            )
            weak_count += 1

        chromosome_rows.append(
            {
                "comparison": comparison,
                "comparison_class": (
                    cfg["comparison_class"]
                ),
                "reference_species": reference,
                "target_species": target,
                "reference_chromosome": ref_chr,
                "reference_display": (
                    ref_display[ref_chr]
                ),
                "reference_ordered_genes": str(
                    ref_counts[ref_chr]
                ),
                "reference_syntenic_genes": str(
                    len(syntenic_genes)
                ),
                "top1_target_chromosome": (
                    top1_chr
                ),
                "top1_target_display": (
                    target_display[top1_chr]
                ),
                "top1_support_genes": str(
                    top1_support
                ),
                "top2_target_chromosome": (
                    top2_chr
                ),
                "top2_target_display": (
                    target_display[top2_chr]
                    if top2_chr != "NA"
                    else "NA"
                ),
                "top2_support_genes": str(
                    top2_support
                ),
                "second_to_first_support_ratio": (
                    f"{second_first_ratio:.8f}"
                ),
                "top1_top2_union_genes": str(
                    len(union_genes)
                ),
                "top1_top2_union_fraction": (
                    f"{union_fraction:.8f}"
                ),
                "dual_top1_top2_genes": str(
                    len(dual_top2_genes)
                ),
                "dual_top1_top2_fraction": (
                    f"{dual_top2_fraction:.8f}"
                ),
                "genes_exactly_2_target_chromosomes": str(
                    genes_exactly_2_target_chr
                ),
                "genes_ge2_target_chromosomes": str(
                    genes_ge2_target_chr
                ),
                "classification": classification,
                "status": "PASS",
            }
        )

        all_chr_rows.append(
            chromosome_rows[-1]
        )

    if not chromosome_rows:
        raise SystemExit(
            f"ERROR: no chromosome rows "
            f"for {comparison}"
        )

    all_syntenic_genes = set(
        ref_gene_to_targets
    )

    total_syntenic = len(
        all_syntenic_genes
    )

    exactly_two_chr = sum(
        len(
            ref_gene_to_target_chrs[
                gene
            ]
        ) == 2
        for gene in all_syntenic_genes
    )

    ge2_chr = sum(
        len(
            ref_gene_to_target_chrs[
                gene
            ]
        ) >= 2
        for gene in all_syntenic_genes
    )

    ratios = [
        float(
            row[
                "second_to_first_support_ratio"
            ]
        )
        for row in chromosome_rows
    ]

    unions = [
        float(
            row[
                "top1_top2_union_fraction"
            ]
        )
        for row in chromosome_rows
    ]

    duals = [
        float(
            row[
                "dual_top1_top2_fraction"
            ]
        )
        for row in chromosome_rows
    ]

    def median(values):
        ordered = sorted(values)
        n = len(ordered)

        if n % 2:
            return ordered[n // 2]

        return (
            ordered[n // 2 - 1]
            + ordered[n // 2]
        ) / 2

    comparison_rows.append(
        {
            "comparison": comparison,
            "comparison_class": (
                cfg["comparison_class"]
            ),
            "reference_species": reference,
            "target_species": target,
            "reference_chromosomes_tested": str(
                len(chromosome_rows)
            ),
            "reference_syntenic_genes": str(
                total_syntenic
            ),
            "genes_exactly_2_target_chromosomes": str(
                exactly_two_chr
            ),
            "genes_ge2_target_chromosomes": str(
                ge2_chr
            ),
            "fraction_ge2_target_chromosomes": (
                f"{ge2_chr / total_syntenic:.8f}"
            ),
            "strong_1to2_candidate_chromosomes": str(
                strong_count
            ),
            "moderate_balanced_chromosomes": str(
                moderate_count
            ),
            "weak_second_partner_chromosomes": str(
                weak_count
            ),
            "single_target_chromosomes": str(
                single_count
            ),
            "strong_candidate_fraction": (
                f"{strong_count / len(chromosome_rows):.8f}"
            ),
            "median_second_to_first_ratio": (
                f"{median(ratios):.8f}"
            ),
            "median_top1_top2_union_fraction": (
                f"{median(unions):.8f}"
            ),
            "median_dual_top1_top2_fraction": (
                f"{median(duals):.8f}"
            ),
            "strong_definition": (
                "second:first>=0.50_AND_union>=0.70"
            ),
            "status": "PASS",
        }
    )

    qc_rows.append(
        {
            "comparison": comparison,
            "unique_pair_rows": str(
                len(pairs)
            ),
            "unique_reference_genes": str(
                total_syntenic
            ),
            "reference_chromosomes_tested": str(
                len(chromosome_rows)
            ),
            "reference_chromosome_count_lens": str(
                len(ref_order)
            ),
            "target_chromosome_count_lens": str(
                len(target_order)
            ),
            "status": "PASS",
        }
    )

    print(
        f"{comparison}: "
        f"{total_syntenic} syntenic reference genes; "
        f"{ge2_chr} with >=2 target chromosomes "
        f"({ge2_chr / total_syntenic:.4%}); "
        f"strong={strong_count}, "
        f"moderate={moderate_count}, "
        f"weak={weak_count}.",
        flush=True,
    )

# Main TSV outputs
write_tsv(
    summary_path,
    comparison_rows,
)

write_tsv(
    chr_summary_path,
    all_chr_rows,
)

write_tsv(
    matrix_summary_path,
    all_matrix_rows,
)

write_tsv(
    qc_path,
    qc_rows,
)

if len(comparison_rows) != 3:
    raise SystemExit(
        f"ERROR: expected 3 X12 comparisons; "
        f"found {len(comparison_rows)}"
    )

###########################################################################
# Heatmaps
###########################################################################

try:
    import matplotlib
    matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    import numpy as np

except Exception as error:
    raise SystemExit(
        "ERROR: matplotlib/numpy unavailable "
        f"for X12 heatmaps: {error}"
    )

for cfg in COMPARISONS:

    comparison = cfg["comparison"]
    reference = cfg["reference"]
    target = cfg["target"]

    ref_order, _ = read_lens(
        lens_dir / f"{reference}.wgdi.lens"
    )

    target_order, _ = read_lens(
        lens_dir / f"{target}.wgdi.lens"
    )

    ref_display = display_map(
        ref_order
    )

    target_display = display_map(
        target_order
    )

    rows = [
        row
        for row in all_matrix_rows
        if row["comparison"] == comparison
    ]

    matrix = np.zeros(
        (
            len(ref_order),
            len(target_order),
        ),
        dtype=float,
    )

    ref_index = {
        chromosome: i
        for i, chromosome in enumerate(
            ref_order
        )
    }

    target_index = {
        chromosome: i
        for i, chromosome in enumerate(
            target_order
        )
    }

    for row in rows:
        matrix[
            ref_index[
                row["reference_chromosome"]
            ],
            target_index[
                row["target_chromosome"]
            ],
        ] = float(
            row["supporting_reference_genes"]
        )

    fig_width = max(
        8.0,
        len(target_order) * 0.70,
    )

    fig_height = max(
        6.0,
        len(ref_order) * 0.65,
    )

    fig, ax = plt.subplots(
        figsize=(
            fig_width,
            fig_height,
        )
    )

    image = ax.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
    )

    ax.set_xticks(
        range(len(target_order))
    )

    ax.set_xticklabels(
        [
            target_display[c]
            for c in target_order
        ],
        rotation=45,
        ha="right",
    )

    ax.set_yticks(
        range(len(ref_order))
    )

    ax.set_yticklabels(
        [
            ref_display[c]
            for c in ref_order
        ]
    )

    ax.set_xlabel(
        f"{target} target chromosomes"
    )

    ax.set_ylabel(
        f"{reference} reference chromosomes"
    )

    ax.set_title(
        f"{comparison}: syntenic reference-gene support"
    )

    for i in range(
        len(ref_order)
    ):
        for j in range(
            len(target_order)
        ):
            value = int(
                matrix[i, j]
            )

            if value > 0:
                ax.text(
                    j,
                    i,
                    str(value),
                    ha="center",
                    va="center",
                    fontsize=7,
                )

    ax.set_xticks(
        [
            x - 0.5
            for x in range(
                1,
                len(target_order),
            )
        ],
        minor=True,
    )

    ax.set_yticks(
        [
            y - 0.5
            for y in range(
                1,
                len(ref_order),
            )
        ],
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

    cbar = fig.colorbar(
        image,
        ax=ax,
    )

    cbar.set_label(
        "Reference genes with a syntenic partner"
    )

    fig.tight_layout()

    fig.savefig(
        figure_dir
        / f"{comparison}.chromosome_support_heatmap.pdf",
        bbox_inches="tight",
    )

    fig.savefig(
        figure_dir
        / f"{comparison}.chromosome_support_heatmap.svg",
        bbox_inches="tight",
    )

    plt.close(fig)

print(
    "Step 36X12 chromosome-consistent multisynteny: PASS"
)
PY

PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            n++
        }
        END {
            print n + 0
        }
    ' "${SUMMARY}"
)"

if [[ "${PASS_COUNT}" -ne 3 ]]; then
    echo "ERROR: ${PASS_COUNT}/3 X12 comparisons passed." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36X12_COMPLETE.txt" <<EOF2
checkpoint=step36X12_additional_chromosome_consistent_multisynteny
date=$(date --iso-8601=seconds)
comparisons_expected=3
comparisons_complete=${PASS_COUNT}
comparisons=VSCU_VSER,VSCU_VPAN,VPAN_VPER
source=step36X11_exact_unique_gene_pairs
reference_orientation=VSCU_to_VSER,VSCU_to_VPAN,VPAN_to_VPER
top_partner_measure=unique_reference_genes_supporting_reference_chr_target_chr_pair
strong_candidate_second_to_first_threshold=0.50
strong_candidate_top1_top2_union_threshold=0.70
strong_classification=descriptive_not_final_WGD_call
display_chromosome_labels=Chr1_ChrN_in_WGDI_lens_order
summary=11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_comparison_summary.tsv
chromosome_summary=11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_reference_chromosome_summary.tsv
matrix_summary=11_wgdi/13_chromosome_consistent_multisynteny_additional/01_tables/step36X12_chromosome_pair_support.tsv
status=PASS
next_step=step36H3_Ks_component_specific_multisynteny
EOF2

cp -f \
    "${SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${MATRIX_SUMMARY}" \
    "${QC_TABLE}" \
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
echo "Step 36X12 comparison summary"
echo "============================================================"

column -t -s $'\t' "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36X12 checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36X12_COMPLETE.txt"
