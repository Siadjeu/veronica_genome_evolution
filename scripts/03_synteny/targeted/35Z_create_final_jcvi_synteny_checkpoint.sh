#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=jcvi_check
#SBATCH --output=10_synteny/logs/jcvi_checkpoint_%j.out
#SBATCH --error=10_synteny/logs/jcvi_checkpoint_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

CHECKPOINT_ROOT="${SYNTENY_DIR}/final_checkpoint_jcvi_synteny"

METHODS_DIR="${CHECKPOINT_ROOT}/01_methods"
SCRIPTS_DIR="${CHECKPOINT_ROOT}/02_scripts"
MANIFEST_DIR="${CHECKPOINT_ROOT}/03_manifests"
QC_DIR="${CHECKPOINT_ROOT}/04_qc"
TABLES_DIR="${CHECKPOINT_ROOT}/05_tables"
FIGURES_DIR="${CHECKPOINT_ROOT}/06_figures"
ANCHOR_SUMMARY_DIR="${CHECKPOINT_ROOT}/07_anchor_summaries"
ENVIRONMENT_DIR="${CHECKPOINT_ROOT}/08_environment"
INVENTORY_DIR="${CHECKPOINT_ROOT}/09_inventory"

LOG_DIR="${SYNTENY_DIR}/logs"

mkdir -p \
    "${LOG_DIR}" \
    "${CHECKPOINT_ROOT}"

cd "${PROJECT_DIR}"

# ============================================================
# Activate JCVI environment
# ============================================================

module purge
module load hpc-env/13.1

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
elif [[ -x "${HOME}/anaconda3/bin/conda" ]]; then
    CONDA_EXE_PATH="${HOME}/anaconda3/bin/conda"
else
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate jcvi_env

echo "Python:"
python --version

echo "JCVI:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY

# ============================================================
# Validate major checkpoints before freezing results
# ============================================================

REQUIRED_CHECKPOINTS=(
    "10_synteny/checkpoint_jcvi_inputs/JCVI_INPUTS_COMPLETE.txt"
    "10_synteny/checkpoint_jcvi_step34/JCVI_STEP34_COMPLETE.txt"
    "10_synteny/checkpoint_macrosynteny_all_species/ALL_SPECIES_MACROSYNTENY_COMPLETE.txt"
    "10_synteny/checkpoint_refined_macrosynteny_story/REFINED_MACROSYNTENY_STORY_COMPLETE.txt"
    "10_synteny/checkpoint_missing_tree_order_comparisons/MISSING_TREE_ORDER_COMPARISONS_COMPLETE.txt"
    "10_synteny/checkpoint_refined_macrosynteny_plots/REFINED_MACROSYNTENY_PLOTS_COMPLETE.txt"
    "10_synteny/checkpoint_publication_macrosynteny_plots/PUBLICATION_MACROSYNTENY_PLOTS_COMPLETE.txt"
)

for FILE in "${REQUIRED_CHECKPOINTS[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required synteny checkpoint missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi

    if ! grep -q '^status=PASS$' "${FILE}"; then
        echo "ERROR: Checkpoint does not contain status=PASS:" >&2
        echo "${FILE}" >&2
        cat "${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Rebuild checkpoint cleanly
# ============================================================

rm -rf "${CHECKPOINT_ROOT}"

mkdir -p \
    "${METHODS_DIR}" \
    "${SCRIPTS_DIR}" \
    "${MANIFEST_DIR}" \
    "${QC_DIR}" \
    "${TABLES_DIR}" \
    "${FIGURES_DIR}/overview" \
    "${FIGURES_DIR}/publication_style" \
    "${FIGURES_DIR}/pairwise_dotplots" \
    "${ANCHOR_SUMMARY_DIR}" \
    "${ENVIRONMENT_DIR}" \
    "${INVENTORY_DIR}"

# ============================================================
# Helper functions
# ============================================================

copy_required_file() {
    local SOURCE="$1"
    local DESTINATION_DIR="$2"

    if [[ ! -s "${SOURCE}" ]]; then
        echo "ERROR: Required file missing or empty:" >&2
        echo "${SOURCE}" >&2
        exit 1
    fi

    cp -p "${SOURCE}" "${DESTINATION_DIR}/"
}

copy_optional_file() {
    local SOURCE="$1"
    local DESTINATION_DIR="$2"

    if [[ -s "${SOURCE}" ]]; then
        cp -p "${SOURCE}" "${DESTINATION_DIR}/"
    else
        printf '%s\t%s\n' \
            "${SOURCE}" \
            "MISSING_OPTIONAL" \
            >> "${INVENTORY_DIR}/missing_optional_files.tsv"
    fi
}

copy_optional_pattern() {
    local PATTERN="$1"
    local DESTINATION_DIR="$2"

    shopt -s nullglob

    local FILES=( ${PATTERN} )

    shopt -u nullglob

    if [[ "${#FILES[@]}" -eq 0 ]]; then
        printf '%s\t%s\n' \
            "${PATTERN}" \
            "NO_MATCH_OPTIONAL" \
            >> "${INVENTORY_DIR}/missing_optional_files.tsv"

        return
    fi

    for SOURCE in "${FILES[@]}"
    do
        if [[ -f "${SOURCE}" && -s "${SOURCE}" ]]; then
            cp -p "${SOURCE}" "${DESTINATION_DIR}/"
        fi
    done
}

printf 'path\tstatus\n' \
    > "${INVENTORY_DIR}/missing_optional_files.tsv"

# ============================================================
# Copy workflow scripts
# ============================================================

SCRIPT_CANDIDATES=(
    "05_scripts/32A_discover_synteny_inputs.sh"
    "05_scripts/32B_prepare_chromosome_filtered_synteny_inputs.sh"
    "05_scripts/32B2_map_gtf_seqids_to_genome_by_md5.sh"
    "05_scripts/32B3_rebuild_verified_chromosome_filtered_inputs.sh"
    "05_scripts/33_prepare_validate_jcvi_inputs.sh"
    "05_scripts/34A_prepare_jcvi_comparison_manifest.sh"
    "05_scripts/34B_run_jcvi_comparisons.sh"
    "05_scripts/34C_summarize_jcvi_comparisons.sh"
    "05_scripts/34d_generate_jcvi_synteny_dotplots.sh"
    "05_scripts/34e_all_species_macrosynteny_numeric_chromosomes.sh"
    "05_scripts/34e2_resume_all_species_macrosynteny_plot.sh"
    "05_scripts/34f_prepare_refined_macrosynteny_story.sh"
    "05_scripts/34g_run_missing_tree_order_comparisons.sh"
    "05_scripts/34g2_summarize_missing_tree_order_comparisons.sh"
    "05_scripts/34h_generate_refined_macrosynteny_plots.sh"
    "05_scripts/34i_generate_publication_macrosynteny_plots.sh"
    "05_scripts/35_quantify_syntenic_depth_and_coverage.sh"
)

for SCRIPT in "${SCRIPT_CANDIDATES[@]}"
do
    copy_optional_file \
        "${SCRIPT}" \
        "${SCRIPTS_DIR}"
done

# Copy preserved script versions, where present.
copy_optional_pattern \
    "05_scripts/34f_prepare_refined_macrosynteny_story.sh.*" \
    "${SCRIPTS_DIR}"

copy_optional_pattern \
    "05_scripts/34h_generate_refined_macrosynteny_plots.sh.*" \
    "${SCRIPTS_DIR}"

copy_optional_pattern \
    "05_scripts/34i_generate_publication_macrosynteny_plots.sh.*" \
    "${SCRIPTS_DIR}"

# ============================================================
# Copy manifests and orders
# ============================================================

copy_required_file \
    "10_synteny/manifests/jcvi_input_manifest.tsv" \
    "${MANIFEST_DIR}"

copy_optional_file \
    "10_synteny/manifests/synteny_input_manifest.verified.tsv" \
    "${MANIFEST_DIR}"

copy_optional_pattern \
    "10_synteny/manifests/*.tsv" \
    "${MANIFEST_DIR}"

copy_required_file \
    "10_synteny/refined_macrosynteny/orders/full_9_species_plot_order.tsv" \
    "${MANIFEST_DIR}"

copy_required_file \
    "10_synteny/refined_macrosynteny/orders/veronica_only_plot_order.tsv" \
    "${MANIFEST_DIR}"

copy_required_file \
    "10_synteny/refined_macrosynteny/missing_comparisons/manifests/step34g_comparisons.tsv" \
    "${MANIFEST_DIR}"

copy_optional_pattern \
    "10_synteny/jcvi_step34/manifests/*.tsv" \
    "${MANIFEST_DIR}"

# ============================================================
# Copy checkpoint and QC files
# ============================================================

for CHECKPOINT_FILE in "${REQUIRED_CHECKPOINTS[@]}"
do
    cp -p \
        "${CHECKPOINT_FILE}" \
        "${QC_DIR}/"
done

copy_optional_pattern \
    "10_synteny/checkpoint_jcvi_inputs/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/checkpoint_jcvi_step34/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/checkpoint_macrosynteny_all_species/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/checkpoint_refined_macrosynteny_story/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/checkpoint_missing_tree_order_comparisons/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/checkpoint_refined_macrosynteny_plots/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/checkpoint_publication_macrosynteny_plots/*" \
    "${QC_DIR}"

copy_optional_pattern \
    "10_synteny/refined_macrosynteny/tables/*mapping_qc*.tsv" \
    "${QC_DIR}"

copy_optional_file \
    "10_synteny/refined_macrosynteny/tables/refined_numeric_chromosome_mapping.tsv" \
    "${QC_DIR}"

copy_optional_file \
    "10_synteny/jcvi_step34/tables/jcvi_step34_results.tsv" \
    "${QC_DIR}"

copy_optional_file \
    "10_synteny/jcvi_step34/tables/jcvi_step34_summary.tsv" \
    "${QC_DIR}"

copy_optional_file \
    "10_synteny/refined_macrosynteny/missing_comparisons/tables/step34g_missing_comparison_results.tsv" \
    "${QC_DIR}"

copy_optional_file \
    "10_synteny/refined_macrosynteny/missing_comparisons/tables/step34g_summary.tsv" \
    "${QC_DIR}"

# ============================================================
# Copy major result tables
# ============================================================

RESULT_TABLES=(
    "10_synteny/jcvi_step34/tables/jcvi_step34_results.tsv"
    "10_synteny/jcvi_step34/tables/jcvi_step34_summary.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_storyline_summary.tsv"
    "10_synteny/refined_macrosynteny/tables/full_9_species_required_adjacent_comparisons.tsv"
    "10_synteny/refined_macrosynteny/tables/veronica_only_required_adjacent_comparisons.tsv"
    "10_synteny/refined_macrosynteny/tables/missing_tree_order_adjacent_comparisons.tsv"
    "10_synteny/refined_macrosynteny/tables/chromosome_pair_anchor_counts.tsv"
    "10_synteny/refined_macrosynteny/tables/synteny_ratio_summary.tsv"
    "10_synteny/refined_macrosynteny/tables/chromosome_rearrangement_candidates.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_full_plot_edge_mapping_qc.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_veronica_only_edge_mapping_qc.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_full_plot_synteny_ratios.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_veronica_only_plot_synteny_ratios.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_full_plot_rearrangement_candidates.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_veronica_only_plot_rearrangement_candidates.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_macrosynteny_plot_summary.tsv"
    "10_synteny/refined_macrosynteny/tables/publication_full_plot_synteny_ratios.tsv"
    "10_synteny/refined_macrosynteny/tables/publication_veronica_plot_synteny_ratios.tsv"
    "10_synteny/refined_macrosynteny/tables/publication_full_plot_rearrangement_pairs.tsv"
    "10_synteny/refined_macrosynteny/tables/publication_veronica_plot_rearrangement_pairs.tsv"
    "10_synteny/refined_macrosynteny/tables/publication_macrosynteny_plot_summary.tsv"
    "10_synteny/refined_macrosynteny/tables/refined_numeric_chromosome_mapping.tsv"
)

for TABLE in "${RESULT_TABLES[@]}"
do
    copy_optional_file \
        "${TABLE}" \
        "${TABLES_DIR}"
done

# Copy Step 35 outputs if they already exist.
copy_optional_pattern \
    "10_synteny/syntenic_depth/tables/*.tsv" \
    "${TABLES_DIR}"

copy_optional_pattern \
    "10_synteny/jcvi_step35/tables/*.tsv" \
    "${TABLES_DIR}"

copy_optional_pattern \
    "10_synteny/*depth*/tables/*.tsv" \
    "${TABLES_DIR}"

# ============================================================
# Copy final figures
# ============================================================

copy_optional_pattern \
    "10_synteny/jcvi_step34/macrosynteny_all_species/plots/*" \
    "${FIGURES_DIR}/overview"

copy_optional_pattern \
    "10_synteny/refined_macrosynteny/plots/full_9_species.macrosynteny.refined.*" \
    "${FIGURES_DIR}/overview"

copy_optional_pattern \
    "10_synteny/refined_macrosynteny/plots/veronica_only_8_species.macrosynteny.refined.*" \
    "${FIGURES_DIR}/overview"

copy_required_file \
    "10_synteny/refined_macrosynteny/publication_plots/full_9_species.macrosynteny.publication_style.pdf" \
    "${FIGURES_DIR}/publication_style"

copy_required_file \
    "10_synteny/refined_macrosynteny/publication_plots/full_9_species.macrosynteny.publication_style.svg" \
    "${FIGURES_DIR}/publication_style"

copy_required_file \
    "10_synteny/refined_macrosynteny/publication_plots/full_9_species.macrosynteny.publication_style.png" \
    "${FIGURES_DIR}/publication_style"

copy_required_file \
    "10_synteny/refined_macrosynteny/publication_plots/veronica_only_8_species.macrosynteny.publication_style.pdf" \
    "${FIGURES_DIR}/publication_style"

copy_required_file \
    "10_synteny/refined_macrosynteny/publication_plots/veronica_only_8_species.macrosynteny.publication_style.svg" \
    "${FIGURES_DIR}/publication_style"

copy_required_file \
    "10_synteny/refined_macrosynteny/publication_plots/veronica_only_8_species.macrosynteny.publication_style.png" \
    "${FIGURES_DIR}/publication_style"

copy_optional_pattern \
    "10_synteny/jcvi_step34/dotplots/**/*" \
    "${FIGURES_DIR}/pairwise_dotplots"

copy_optional_pattern \
    "10_synteny/jcvi_step34/plots/*" \
    "${FIGURES_DIR}/pairwise_dotplots"

# ============================================================
# Build compact anchor comparison summaries
# ============================================================

python - \
    "${PROJECT_DIR}" \
    "${ANCHOR_SUMMARY_DIR}" \
    "${TABLES_DIR}" <<'PY'
from __future__ import annotations

import csv
import sys
from pathlib import Path

project_dir = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
tables_dir = Path(sys.argv[3])

step34_file = (
    project_dir
    / "10_synteny"
    / "jcvi_step34"
    / "tables"
    / "jcvi_step34_results.tsv"
)

step34g_file = (
    project_dir
    / "10_synteny"
    / "refined_macrosynteny"
    / "missing_comparisons"
    / "tables"
    / "step34g_missing_comparison_results.tsv"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True,
)


def read_tsv(path):
    if not path.is_file():
        return []

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


def write_tsv(
    path,
    rows,
    fieldnames,
):
    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


rows = []

for row in read_tsv(step34_file):
    if row.get("status") != "PASS":
        continue

    rows.append(
        {
            "source_step": "34",
            "comparison_id": row.get(
                "comparison_id",
                "",
            ),
            "comparison_type": row.get(
                "comparison_type",
                "",
            ),
            "query_species": row.get(
                "query_species",
                "",
            ),
            "subject_species": row.get(
                "subject_species",
                "",
            ),
            "anchor_file": row.get(
                "anchor_file",
                "",
            ),
            "anchor_pair_count": row.get(
                "anchor_pair_count",
                row.get(
                    "anchor_pairs",
                    "",
                ),
            ),
            "anchor_block_count": row.get(
                "anchor_block_count",
                row.get(
                    "anchor_blocks",
                    "",
                ),
            ),
            "simple_anchor_file": row.get(
                "simple_anchor_file",
                "",
            ),
            "simple_anchor_pair_count": row.get(
                "simple_anchor_pair_count",
                row.get(
                    "simple_anchor_pairs",
                    "",
                ),
            ),
            "status": "PASS",
        }
    )

for row in read_tsv(step34g_file):
    if row.get("status") != "PASS":
        continue

    rows.append(
        {
            "source_step": "34G",
            "comparison_id": row.get(
                "comparison_id",
                "",
            ),
            "comparison_type": (
                "tree_order_additional"
            ),
            "query_species": row.get(
                "query_species",
                "",
            ),
            "subject_species": row.get(
                "subject_species",
                "",
            ),
            "anchor_file": row.get(
                "anchor_file",
                "",
            ),
            "anchor_pair_count": row.get(
                "anchor_pair_count",
                "",
            ),
            "anchor_block_count": row.get(
                "anchor_block_count",
                "",
            ),
            "simple_anchor_file": row.get(
                "simple_anchor_file",
                "",
            ),
            "simple_anchor_pair_count": row.get(
                "simple_anchor_pair_count",
                "",
            ),
            "status": "PASS",
        }
    )

rows.sort(
    key=lambda row: (
        row["source_step"],
        row["comparison_id"],
    )
)

fields = [
    "source_step",
    "comparison_id",
    "comparison_type",
    "query_species",
    "subject_species",
    "anchor_file",
    "anchor_pair_count",
    "anchor_block_count",
    "simple_anchor_file",
    "simple_anchor_pair_count",
    "status",
]

write_tsv(
    output_dir
    / "all_successful_jcvi_comparisons.tsv",
    rows,
    fields,
)

summary_rows = [
    {
        "metric": "successful_comparisons",
        "value": len(rows),
    },
    {
        "metric": "step34_comparisons",
        "value": sum(
            row["source_step"] == "34"
            for row in rows
        ),
    },
    {
        "metric": "step34g_additional_comparisons",
        "value": sum(
            row["source_step"] == "34G"
            for row in rows
        ),
    },
    {
        "metric": "status",
        "value": "PASS",
    },
]

write_tsv(
    output_dir
    / "anchor_comparison_summary.tsv",
    summary_rows,
    [
        "metric",
        "value",
    ],
)

print(
    "Successful anchor comparisons:",
    len(rows),
)
PY

# ============================================================
# Create methods and results documents
# ============================================================

cat > "${METHODS_DIR}/JCVI_SYNTENY_METHODS.md" <<'EOF2'
# JCVI synteny and macrosynteny methods

## Study design

Chromosome-scale synteny was analysed across nine species:

- PMAJ: Plantago major, 2x, six chromosomes
- VPAN: Veronica panormitana, 2x, nine chromosomes
- VSCU: Veronica scutellata, 2x, nine chromosomes
- VANA: Veronica anagallis-aquatica, 4x, eighteen chromosomes
- VARV: Veronica arvensis, 2x, eight chromosomes
- VVER: Veronica verna, 2x, eight chromosomes
- VPER: Veronica persica, 4x, fourteen chromosomes
- VSER: Veronica serpyllifolia, 2x, seven chromosomes
- VTRI: Veronica triphyllos, 2x, seven chromosomes

Plantago major was used as the diploid outgroup.

## Input preparation

One representative protein-coding transcript per locus was retained for each species. Gene coordinates and protein identifiers were validated for exact agreement before synteny analysis. Chromosome identifiers from genome FASTA files and BRAKER-derived annotations were reconciled using exact sequence MD5 matching. Only verified chromosome-associated genes were retained.

For each species, JCVI-compatible BED and protein FASTA files were generated. The BED files were checked for:

- duplicate gene identifiers;
- invalid coordinates;
- chromosome-order consistency;
- within-chromosome gene-order consistency;
- BED/protein identifier agreement;
- expected chromosome count.

## Homology search and collinearity

Pairwise and self-collinearity analyses were performed with JCVI version 1.6.5. Protein similarities were calculated with DIAMOND version 2.2.4 using the JCVI `diamond_blastp` backend.

JCVI ortholog cataloguing was run with:

- protein sequence input;
- `diamond_blastp` alignment backend;
- C-score threshold of 0.70 for the additional tree-order comparisons;
- no stripping of gene names;
- sixteen CPUs for computational comparisons.

Full JCVI anchor files were retained for gene-level and chromosome-level analyses. Simplified anchor files were generated using a minimum block span of five genes for plotting and block summaries.

## Comparison design

The primary JCVI workflow included:

- nine within-species self-comparisons;
- sixteen ordinary pairwise comparisons among Veronica species;
- three comparisons involving the outgroup;
- three additional comparisons required to connect adjacent taxa in tree-based plotting order.

The three additional comparisons were:

- VARV versus VPAN;
- VPAN versus VTRI;
- VTRI versus VSCU.

## Phylogenetic display order

The rooted OrthoFinder species tree was used to derive a reproducible left-to-right terminal-taxon order. Plantago major was placed first as the outgroup in the full figure.

The full plotting order was:

PMAJ, VANA, VSER, VPER, VVER, VARV, VPAN, VTRI, VSCU.

The Veronica-only plotting order was:

VANA, VSER, VPER, VVER, VARV, VPAN, VTRI, VSCU.

This is a display order derived from the rooted tree. Branch rotations could alter the left-to-right order without changing the underlying topology.

## Chromosome-number standardization

Analytical chromosome identifiers were preserved in the original validated inputs. A separate plotting-only mapping replaced chromosome identifiers with consecutive numeric labels. Chromosome lengths and gene positions were retained.

## Full-anchor mapping

Full JCVI anchor pairs were mapped back to gene coordinates in the validated BED files. For each adjacent species comparison, the proportion of full anchor pairs successfully assigned to both genomes was calculated. Plotting required a minimum mapping fraction of 0.80. All final adjacent comparisons achieved a mapping fraction of 1.00.

## Multiple strong-partner relationships

For each chromosome in a species comparison, anchor counts were calculated against every chromosome in the second species.

A chromosome partner was classified as strong when:

1. the chromosome pair contained at least 20 full-anchor gene pairs; and
2. the pair contained at least 25 percent of the anchor count of the strongest partner for that focal chromosome.

A chromosome pair was classified as participating in a multiple strong-partner relationship when either chromosome had more than one strong partner.

These relationships may result from:

- chromosome fusion;
- chromosome fission;
- translocation;
- retained duplicated or homeologous chromosomes;
- ancient duplicated regions;
- post-polyploid chromosome restructuring.

They were therefore not automatically classified as specific fusion, fission or inversion events.

## Operational chromosome-component ratios

Strong chromosome-partner relationships were represented as a bipartite network. Connected components were identified, and the number of chromosomes from the upper and lower species in each component was counted.

The component receiving the greatest support from anchor-pair counts was reported as the operational chromosome-component ratio.

These ratios are exploratory summaries of chromosome connectivity. They are not treated as direct ploidy ratios and require integration with:

- self-synteny;
- syntenic depth;
- Ks distributions;
- gene-tree duplication mapping;
- ancestral chromosome reconstruction.

## Macrosynteny visualization

Publication-style macrosynteny figures were generated from mapped full-anchor pairs.

- chromosome bars were assigned species-specific colors;
- ploidy was displayed explicitly as 2x or 4x before each species abbreviation;
- blue links represented mapped full-anchor synteny pairs;
- red links represented chromosome pairs participating in multiple strong-partner relationships;
- operational chromosome-component ratios were displayed for adjacent comparisons;
- chromosome identifiers were displayed numerically;
- separate full nine-species and Veronica-only figures were generated.

Figures were exported as PDF, SVG and 400-dpi PNG files.

## Reproducibility

All scripts, manifests, QC files, result tables, figures, environment records, inventories and SHA-256 checksums were preserved in the final JCVI synteny checkpoint.
EOF2

cat > "${METHODS_DIR}/JCVI_SYNTENY_RESULTS_SUMMARY.md" <<'EOF2'
# JCVI synteny results summary

## Input validation

All nine species passed the final JCVI input validation. Across the validated JCVI input set, 413,482 BED genes had matching protein identifiers. No duplicate identifiers, invalid coordinates or gene-order failures were detected.

## JCVI comparison completion

The primary JCVI comparison set completed successfully:

- 28 comparisons passed;
- zero comparisons failed;
- zero comparisons were missing;
- nine self-comparisons passed;
- sixteen ordinary pairwise comparisons passed;
- three outgroup comparisons passed.

The primary comparison set contained:

- 940,442 full anchor pairs;
- 39,982 anchor blocks;
- 39,377 simplified anchor pairs.

Three additional tree-order comparisons also passed:

- VARV-VPAN: 21,961 full anchor pairs and 488 simplified anchor pairs;
- VPAN-VTRI: 21,521 full anchor pairs and 513 simplified anchor pairs;
- VTRI-VSCU: 21,721 full anchor pairs and 725 simplified anchor pairs.

Together, these three additional comparisons contained 65,203 full anchor pairs and 1,726 simplified anchor pairs.

## Anchor mapping

All full anchor pairs used for the final tree-ordered nine-species plot mapped successfully to both validated BED files.

The final mapping fractions were 1.00 for all eight adjacent comparisons:

- PMAJ-VANA;
- VANA-VSER;
- VSER-VPER;
- VPER-VVER;
- VVER-VARV;
- VARV-VPAN;
- VPAN-VTRI;
- VTRI-VSCU.

This confirmed that the initial empty plot was a plotting-input issue rather than an identifier or synteny failure.

## Macrosynteny patterns

The final macrosynteny analysis identified extensive conserved gene-level synteny among all adjacent comparisons. Multiple strong chromosome-partner relationships were also identified, indicating complex chromosome connectivity.

These relationships are consistent with one or more of:

- retained homeologous relationships;
- chromosome fusion or fission;
- translocation;
- differential diploidization;
- post-polyploid chromosome restructuring.

Specific chromosome rearrangement classes were not assigned from the overview alone.

## Operational ratios

The operational chromosome-component ratios summarize the dominant connected strong-synteny component between adjacent taxa.

Examples observed in the full plot included:

- PMAJ-VANA: 6:18;
- VANA-VSER: 14:5;
- VSER-VPER: 6:12;
- VPER-VVER: 12:7;
- VVER-VARV: 1:1;
- VARV-VPAN: 4:5;
- VPAN-VTRI: 7:5;
- VTRI-VSCU: 5:7.

The PMAJ-VANA and VSER-VPER relationships are compatible with elevated chromosome multiplicity in the tetraploid species, but these component ratios are not direct estimates of syntenic depth or ploidy.

## Relevance to the research question

The JCVI analysis demonstrates that Veronica contains extensive conserved synteny together with substantial chromosome restructuring. VANA and VPER display complex chromosome relationships consistent with elevated ploidy, while nominally diploid Veronica species also exhibit variable chromosome connectivity.

The results support testing the following model:

An ancestral duplication may have occurred within Veronica, followed by extensive diploidization and chromosome restructuring, with additional lineage-specific polyploidization in VANA and VPER.

The JCVI results alone do not prove a genus-wide ancestral WGD. The hypothesis requires integration with:

- WGDI self-collinearity;
- self-syntenic Ks distributions;
- outgroup-relative syntenic depth;
- block-level Ks;
- ancestral karyotype reconstruction;
- phylogenetic duplication mapping.
EOF2

# ============================================================
# Create environment information
# ============================================================

{
    echo -e "field\tvalue"
    echo -e "checkpoint_date\t$(date --iso-8601=seconds)"
    echo -e "hostname\t$(hostname)"
    echo -e "project_directory\t${PROJECT_DIR}"
    echo -e "conda_environment\t${CONDA_DEFAULT_ENV:-unknown}"
    echo -e "python_executable\t$(command -v python)"
    echo -e "python_version\t$(python --version 2>&1)"
    echo -e "jcvi_version\t$(python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY
)"
    echo -e "diamond_executable\t$(command -v diamond || echo unavailable)"
    echo -e "diamond_version\t$(diamond version 2>/dev/null || echo unavailable)"
} > "${ENVIRONMENT_DIR}/software_environment.tsv"

conda list \
    > "${ENVIRONMENT_DIR}/conda_list_jcvi_env.txt"

python -m pip freeze \
    > "${ENVIRONMENT_DIR}/pip_freeze_jcvi_env.txt" \
    2>/dev/null || true

# ============================================================
# Create file inventory
# ============================================================

python - \
    "${CHECKPOINT_ROOT}" \
    "${INVENTORY_DIR}/checkpoint_file_inventory.tsv" <<'PY'
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

checkpoint_root = Path(sys.argv[1])
output_file = Path(sys.argv[2])


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


rows = []

for path in sorted(
    checkpoint_root.rglob("*")
):
    if not path.is_file():
        continue

    if path == output_file:
        continue

    relative = path.relative_to(
        checkpoint_root
    )

    rows.append(
        {
            "relative_path": str(
                relative
            ),
            "category": (
                relative.parts[0]
                if relative.parts
                else ""
            ),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    )

with output_file.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "relative_path",
            "category",
            "bytes",
            "sha256",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)

print(
    "Checkpoint files inventoried:",
    len(rows),
)
PY

# ============================================================
# Create high-level checkpoint summary
# ============================================================

SCRIPT_COUNT=$(
    find "${SCRIPTS_DIR}" \
        -maxdepth 1 \
        -type f \
        | wc -l
)

TABLE_COUNT=$(
    find "${TABLES_DIR}" \
        -maxdepth 1 \
        -type f \
        | wc -l
)

FIGURE_COUNT=$(
    find "${FIGURES_DIR}" \
        -type f \
        | wc -l
)

QC_COUNT=$(
    find "${QC_DIR}" \
        -maxdepth 1 \
        -type f \
        | wc -l
)

TOTAL_FILES=$(
    find "${CHECKPOINT_ROOT}" \
        -type f \
        | wc -l
)

cat > "${CHECKPOINT_ROOT}/README.md" <<EOF2
# Final JCVI synteny checkpoint

This checkpoint freezes the complete JCVI synteny and macrosynteny workflow before the WGDI analysis begins.

## Contents

- 01_methods: manuscript-oriented methods and results summaries
- 02_scripts: scripts used during synteny preprocessing, JCVI analysis, QC and plotting
- 03_manifests: validated species, input, comparison and plotting-order manifests
- 04_qc: workflow checkpoints, input QC and anchor-mapping QC
- 05_tables: synteny, ratio, rearrangement and plotting-result tables
- 06_figures: overview, publication-style and optional pairwise figures
- 07_anchor_summaries: compact inventory of successful JCVI comparisons
- 08_environment: software and Conda environment information
- 09_inventory: complete checkpoint file inventory and missing optional files

## Frozen workflow status

- species analysed: 9
- primary JCVI comparisons passed: 28
- additional tree-order comparisons passed: 3
- final full-anchor mapping fraction: 1.00 for every adjacent plotting comparison
- full nine-species publication figure: PASS
- Veronica-only publication figure: PASS
- checkpoint status: PASS

This checkpoint should not be modified by later WGDI analyses.
EOF2

cat > "${CHECKPOINT_ROOT}/JCVI_SYNTENY_FINAL_CHECKPOINT.txt" <<EOF2
checkpoint=final_jcvi_synteny_methods_scripts_results
date=$(date --iso-8601=seconds)
species=9
primary_jcvi_comparisons_pass=28
additional_tree_order_comparisons_pass=3
scripts_copied=${SCRIPT_COUNT}
tables_copied=${TABLE_COUNT}
figures_copied=${FIGURE_COUNT}
qc_files_copied=${QC_COUNT}
total_checkpoint_files=${TOTAL_FILES}
full_anchor_mapping_qc=PASS
publication_macrosynteny_plots=PASS
status=PASS
next_step=step36A_prepare_validate_wgdi_inputs
EOF2

# Regenerate inventory after README and checkpoint were added.
python - \
    "${CHECKPOINT_ROOT}" \
    "${INVENTORY_DIR}/checkpoint_file_inventory.tsv" <<'PY'
from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

checkpoint_root = Path(sys.argv[1])
output_file = Path(sys.argv[2])


def sha256(path):
    digest = hashlib.sha256()

    with path.open("rb") as handle:
        while True:
            block = handle.read(
                1024 * 1024
            )

            if not block:
                break

            digest.update(block)

    return digest.hexdigest()


rows = []

for path in sorted(
    checkpoint_root.rglob("*")
):
    if not path.is_file():
        continue

    if path == output_file:
        continue

    relative = path.relative_to(
        checkpoint_root
    )

    rows.append(
        {
            "relative_path": str(
                relative
            ),
            "category": (
                relative.parts[0]
                if relative.parts
                else ""
            ),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    )

with output_file.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "relative_path",
            "category",
            "bytes",
            "sha256",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(rows)
PY

# Final whole-checkpoint checksum file.
find "${CHECKPOINT_ROOT}" \
    -type f \
    ! -name "FINAL_SHA256SUMS.txt" \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_ROOT}/FINAL_SHA256SUMS.txt"

# ============================================================
# Final validation
# ============================================================

REQUIRED_FINAL_FILES=(
    "${CHECKPOINT_ROOT}/README.md"
    "${CHECKPOINT_ROOT}/JCVI_SYNTENY_FINAL_CHECKPOINT.txt"
    "${METHODS_DIR}/JCVI_SYNTENY_METHODS.md"
    "${METHODS_DIR}/JCVI_SYNTENY_RESULTS_SUMMARY.md"
    "${MANIFEST_DIR}/jcvi_input_manifest.tsv"
    "${TABLES_DIR}/jcvi_step34_results.tsv"
    "${TABLES_DIR}/publication_full_plot_synteny_ratios.tsv"
    "${FIGURES_DIR}/publication_style/full_9_species.macrosynteny.publication_style.pdf"
    "${FIGURES_DIR}/publication_style/veronica_only_8_species.macrosynteny.publication_style.pdf"
    "${ANCHOR_SUMMARY_DIR}/all_successful_jcvi_comparisons.tsv"
    "${ENVIRONMENT_DIR}/software_environment.tsv"
    "${INVENTORY_DIR}/checkpoint_file_inventory.tsv"
    "${CHECKPOINT_ROOT}/FINAL_SHA256SUMS.txt"
)

for FILE in "${REQUIRED_FINAL_FILES[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required final checkpoint file missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

if ! grep -q '^status=PASS$' \
    "${CHECKPOINT_ROOT}/JCVI_SYNTENY_FINAL_CHECKPOINT.txt"
then
    echo "ERROR: Final checkpoint status is not PASS." >&2
    exit 1
fi

echo
echo "============================================================"
echo "Final JCVI synteny checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_ROOT}/JCVI_SYNTENY_FINAL_CHECKPOINT.txt"

echo
echo "Checkpoint size:"

du -sh \
    "${CHECKPOINT_ROOT}"

echo
echo "Top-level contents:"

find "${CHECKPOINT_ROOT}" \
    -mindepth 1 \
    -maxdepth 1 \
    -printf '%f\n' \
    | sort

echo
echo "Step 35Z completed successfully."
