#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=JCVI35Y2
#SBATCH --output=10_synteny/synteny_optimized_pairwise_additional/logs/step35Y2_%j.out
#SBATCH --error=10_synteny/synteny_optimized_pairwise_additional/logs/step35Y2_%j.err

set -euo pipefail

###############################################################################
# STEP 35Y2
# Additional synteny-optimized pairwise macrosynteny plots
#
# Comparisons:
#   VSCU -> VSER   diploid -> diploid
#   VSCU -> VPAN   diploid -> diploid
#   VPAN -> VPER   diploid -> tetraploid
#
# Reproduces Step 35Y:
#   - numeric chromosome labels only
#   - JCVI full-anchor pairs
#   - target chromosome order optimized by full-anchor similarity
#   - target orientation inferred from anchor-position correlation
#   - blue = normal mapped full-anchor links
#   - red = chromosome pairs with multiple strong partners
#
# Frozen chromosome-order dependency:
#
# PMAJ -> VSCU previously established:
#   VSCU = 7,1,3,9,2,6,5,8,4
#
# That exact VSCU order is reused for:
#   VSCU -> VSER
#   VSCU -> VPAN
#
# VSCU -> VPAN establishes the optimized VPAN order.
# That exact VPAN order is then reused for:
#   VPAN -> VPER
###############################################################################

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

SYNTENY_DIR="10_synteny"

OUTPUT_ROOT="${SYNTENY_DIR}/synteny_optimized_pairwise_additional"
TABLE_DIR="${OUTPUT_ROOT}/tables"
ORDER_DIR="${OUTPUT_ROOT}/orders"
PLOT_DIR="${OUTPUT_ROOT}/plots"
LOG_DIR="${OUTPUT_ROOT}/logs"

CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_step35Y2"

JCVI_BED_DIR="${SYNTENY_DIR}/jcvi_inputs/bed"

CHR_MAPPING="${SYNTENY_DIR}/refined_macrosynteny/tables/refined_numeric_chromosome_mapping.tsv"

STEP34_RESULTS="${SYNTENY_DIR}/jcvi_step34/tables/jcvi_step34_results.tsv"

STEP35Y1_RESULTS="${SYNTENY_DIR}/synteny_optimized_pairwise/additional_prerequisite_comparisons/tables/step35Y1_additional_prerequisite_results.tsv"

mkdir -p \
    "${TABLE_DIR}" \
    "${ORDER_DIR}" \
    "${PLOT_DIR}" \
    "${LOG_DIR}" \
    "${CHECKPOINT_DIR}"

###############################################################################
# Environment
###############################################################################

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

export MPLBACKEND=Agg

echo "============================================================"
echo "Step 35Y2: additional synteny-optimized pairwise plots"
echo "============================================================"

date --iso-8601=seconds

###############################################################################
# Validate required inputs
###############################################################################

REQUIRED_FILES=(
    "${CHR_MAPPING}"
    "${STEP34_RESULTS}"
    "${STEP35Y1_RESULTS}"
    "${JCVI_BED_DIR}/VSCU.bed"
    "${JCVI_BED_DIR}/VSER.bed"
    "${JCVI_BED_DIR}/VPAN.bed"
    "${JCVI_BED_DIR}/VPER.bed"
    "${SYNTENY_DIR}/checkpoint_step35Y/STEP35Y_COMPLETE.txt"
    "${SYNTENY_DIR}/checkpoint_step35Y1/STEP35Y1_COMPLETE.txt"
)

for FILE in "${REQUIRED_FILES[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

if ! grep -q '^status=PASS$' \
    "${SYNTENY_DIR}/checkpoint_step35Y/STEP35Y_COMPLETE.txt"
then
    echo "ERROR: Original Step 35Y checkpoint is not PASS." >&2
    exit 1
fi

if ! grep -q '^status=PASS$' \
    "${SYNTENY_DIR}/checkpoint_step35Y1/STEP35Y1_COMPLETE.txt"
then
    echo "ERROR: Step 35Y1 VSCU-VSER checkpoint is not PASS." >&2
    exit 1
fi

###############################################################################
# Targeted cleanup only
###############################################################################

rm -f \
    "${TABLE_DIR}/VSCU_VSER.chromosome_anchor_matrix.tsv" \
    "${TABLE_DIR}/VSCU_VSER.chromosome_anchor_matrix.long.tsv" \
    "${TABLE_DIR}/VSCU_VPAN.chromosome_anchor_matrix.tsv" \
    "${TABLE_DIR}/VSCU_VPAN.chromosome_anchor_matrix.long.tsv" \
    "${TABLE_DIR}/VPAN_VPER.chromosome_anchor_matrix.tsv" \
    "${TABLE_DIR}/VPAN_VPER.chromosome_anchor_matrix.long.tsv" \
    "${TABLE_DIR}/step35Y2_mapping_qc.tsv" \
    "${TABLE_DIR}/step35Y2_pairwise_summary.tsv" \
    "${TABLE_DIR}/step35Y2_strong_partner_summary.tsv" \
    "${TABLE_DIR}/step35Y2_operational_components.tsv" \
    "${TABLE_DIR}/step35Y2_overall_summary.tsv"

rm -f \
    "${ORDER_DIR}/VSCU_VSER.optimized_chromosome_order.tsv" \
    "${ORDER_DIR}/VSCU_VPAN.optimized_chromosome_order.tsv" \
    "${ORDER_DIR}/VPAN_VPER.optimized_chromosome_order.tsv" \
    "${ORDER_DIR}/step35Y2_all_optimized_orders.tsv"

rm -f \
    "${PLOT_DIR}/VSCU_VSER.synteny_optimized.pdf" \
    "${PLOT_DIR}/VSCU_VSER.synteny_optimized.svg" \
    "${PLOT_DIR}/VSCU_VSER.synteny_optimized.png" \
    "${PLOT_DIR}/VSCU_VPAN.synteny_optimized.pdf" \
    "${PLOT_DIR}/VSCU_VPAN.synteny_optimized.svg" \
    "${PLOT_DIR}/VSCU_VPAN.synteny_optimized.png" \
    "${PLOT_DIR}/VPAN_VPER.synteny_optimized.pdf" \
    "${PLOT_DIR}/VPAN_VPER.synteny_optimized.svg" \
    "${PLOT_DIR}/VPAN_VPER.synteny_optimized.png"

rm -f \
    "${CHECKPOINT_DIR}/STEP35Y2_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/step35Y2_mapping_qc.tsv" \
    "${CHECKPOINT_DIR}/step35Y2_pairwise_summary.tsv" \
    "${CHECKPOINT_DIR}/step35Y2_strong_partner_summary.tsv" \
    "${CHECKPOINT_DIR}/step35Y2_operational_components.tsv" \
    "${CHECKPOINT_DIR}/step35Y2_all_optimized_orders.tsv" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# Main Python analysis
###############################################################################

python - \
    "${PROJECT_ROOT}" \
    "${OUTPUT_ROOT}" \
    "${CHR_MAPPING}" \
    "${STEP34_RESULTS}" \
    "${STEP35Y1_RESULTS}" \
    "${JCVI_BED_DIR}" <<'PY'

from __future__ import annotations

import csv
import math
import sys

from collections import Counter, defaultdict, deque
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyBboxPatch, PathPatch


PROJECT_ROOT = Path(
    sys.argv[1]
).resolve()

OUTPUT_ROOT = Path(
    sys.argv[2]
)

CHR_MAPPING_FILE = Path(
    sys.argv[3]
)

STEP34_RESULTS = Path(
    sys.argv[4]
)

STEP35Y1_RESULTS = Path(
    sys.argv[5]
)

JCVI_BED_DIR = Path(
    sys.argv[6]
)


TABLE_DIR = (
    OUTPUT_ROOT
    / "tables"
)

ORDER_DIR = (
    OUTPUT_ROOT
    / "orders"
)

PLOT_DIR = (
    OUTPUT_ROOT
    / "plots"
)


for directory in [
    TABLE_DIR,
    ORDER_DIR,
    PLOT_DIR,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )


###############################################################################
# Biological and plotting configuration
###############################################################################

SPECIES_NAMES = {
    "VSCU": "Veronica scutellata",
    "VSER": "Veronica serpyllifolia",
    "VPAN": "Veronica panormitana",
    "VPER": "Veronica persica",
}


PLOIDY = {
    "VSCU": "2x",
    "VSER": "2x",
    "VPAN": "2x",
    "VPER": "4x",
}


SPECIES_COLOURS = {
    "VSCU": "#59A14F",
    "VSER": "#4C78A8",
    "VPAN": "#B279A2",
    "VPER": "#E69F00",
}


BLUE_LINK = "#377EB8"
RED_LINK = "#D73027"


# Exact Step-35Y thresholds
MIN_STRONG_ANCHORS = 20
RELATIVE_STRONG_THRESHOLD = 0.25

MAX_LINKS_PER_PAIRWISE_PLOT = 8000

CHROMOSOME_GAP = 0.025
CHROMOSOME_HEIGHT = 0.12


# ------------------------------------------------------------
# CRITICAL:
# Frozen PMAJ-informed VSCU order established in original Step 35Y.
# ------------------------------------------------------------

FROZEN_VSCU_ORDER = [
    "7",
    "1",
    "3",
    "9",
    "2",
    "6",
    "5",
    "8",
    "4",
]


PAIRWISE_COMPARISONS = [
    {
        "comparison_id": "VSCU_VSER",
        "reference": "VSCU",
        "target": "VSER",
        "reference_order_source":
            "PMAJ_informed_VSCU_frozen_from_step35Y",
        "purpose": "diploid_diploid_VSER",
    },
    {
        "comparison_id": "VSCU_VPAN",
        "reference": "VSCU",
        "target": "VPAN",
        "reference_order_source":
            "PMAJ_informed_VSCU_frozen_from_step35Y",
        "purpose": "diploid_diploid_VPAN",
    },
    {
        "comparison_id": "VPAN_VPER",
        "reference": "VPAN",
        "target": "VPER",
        "reference_order_source":
            "VSCU_informed_VPAN_from_VSCU_VPAN",
        "purpose": "diploid_tetraploid_VPER",
    },
]


###############################################################################
# General functions
###############################################################################

def read_tsv(path: Path):

    if (
        not path.is_file()
        or path.stat().st_size == 0
    ):
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
    path: Path,
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


def numeric_sort_key(value):

    try:
        return (
            0,
            int(value),
        )

    except (
        ValueError,
        TypeError,
    ):
        return (
            1,
            str(value),
        )


def safe_int(
    value,
    default=0,
):

    try:
        return int(value)

    except (
        ValueError,
        TypeError,
    ):
        return default


###############################################################################
# Read verified numeric chromosome mapping
###############################################################################

mapping_rows = read_tsv(
    CHR_MAPPING_FILE
)

if not mapping_rows:
    raise SystemExit(
        "ERROR: Chromosome mapping table is empty."
    )


original_to_numeric = defaultdict(
    dict
)

chromosome_lengths = defaultdict(
    dict
)


for row in mapping_rows:

    species = row[
        "species_code"
    ]

    original = row[
        "original_sequence_id"
    ]

    numeric = str(
        row[
            "numeric_chromosome"
        ]
    )

    length_bp = int(
        row[
            "length_bp"
        ]
    )

    if original in original_to_numeric[
        species
    ]:

        raise SystemExit(
            f"ERROR: Duplicate chromosome mapping "
            f"for {species}: {original}"
        )

    original_to_numeric[
        species
    ][original] = numeric

    chromosome_lengths[
        species
    ][numeric] = length_bp


def numeric_chromosomes(
    species,
):

    chromosomes = sorted(
        chromosome_lengths[
            species
        ],
        key=numeric_sort_key,
    )

    if not chromosomes:
        raise SystemExit(
            f"ERROR: No chromosomes found for "
            f"{species}."
        )

    return chromosomes


EXPECTED_CHROMOSOME_COUNTS = {
    "VSCU": 9,
    "VSER": 7,
    "VPAN": 9,
    "VPER": 14,
}


for species, expected in (
    EXPECTED_CHROMOSOME_COUNTS.items()
):

    observed = len(
        numeric_chromosomes(
            species
        )
    )

    if observed != expected:

        raise SystemExit(
            f"ERROR: {species}: "
            f"expected {expected} chromosomes, "
            f"found {observed}."
        )


if (
    set(
        FROZEN_VSCU_ORDER
    )
    !=
    set(
        numeric_chromosomes(
            "VSCU"
        )
    )
):

    raise SystemExit(
        "ERROR: Frozen PMAJ-informed VSCU "
        "order does not exactly match "
        "the verified VSCU chromosome set."
    )


###############################################################################
# Read validated JCVI BED files
###############################################################################

def read_species_bed(
    species,
):

    bed_file = (
        JCVI_BED_DIR
        / f"{species}.bed"
    )

    if not bed_file.is_file():

        raise SystemExit(
            f"ERROR: BED file missing for "
            f"{species}: {bed_file}"
        )


    genes = {}

    chromosome_genes = defaultdict(
        list
    )


    with bed_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line_number, line in enumerate(
            handle,
            start=1,
        ):

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue


            fields = line.rstrip(
                "\n"
            ).split(
                "\t"
            )


            if len(fields) < 4:

                raise SystemExit(
                    f"ERROR: Invalid BED line "
                    f"{bed_file}:{line_number}"
                )


            original_chromosome = (
                fields[0]
            )


            if (
                original_chromosome
                not in original_to_numeric[
                    species
                ]
            ):
                continue


            try:
                start = int(
                    fields[1]
                )

                end = int(
                    fields[2]
                )

            except ValueError:

                raise SystemExit(
                    f"ERROR: Invalid BED coordinates "
                    f"{bed_file}:{line_number}"
                )


            gene_id = fields[
                3
            ]


            if gene_id in genes:

                raise SystemExit(
                    f"ERROR: Duplicate BED gene ID "
                    f"for {species}: {gene_id}"
                )


            chromosome = (
                original_to_numeric[
                    species
                ][
                    original_chromosome
                ]
            )


            midpoint = (
                start + end
            ) / 2.0


            genes[
                gene_id
            ] = {
                "chromosome":
                    chromosome,
                "start":
                    start,
                "end":
                    end,
                "midpoint":
                    midpoint,
            }


            chromosome_genes[
                chromosome
            ].append(
                gene_id
            )


    if not genes:

        raise SystemExit(
            f"ERROR: No BED genes read "
            f"for {species}."
        )


    for chromosome in chromosome_genes:

        chromosome_genes[
            chromosome
        ].sort(
            key=lambda gene_id: (
                genes[
                    gene_id
                ][
                    "start"
                ],
                genes[
                    gene_id
                ][
                    "end"
                ],
                gene_id,
            )
        )


    return (
        genes,
        chromosome_genes,
    )


species_gene_data = {}


for species in [
    "VSCU",
    "VSER",
    "VPAN",
    "VPER",
]:

    species_gene_data[
        species
    ] = read_species_bed(
        species
    )


###############################################################################
# Register successful JCVI comparisons
###############################################################################

comparison_rows = []


for path in [
    STEP34_RESULTS,
    STEP35Y1_RESULTS,
]:

    comparison_rows.extend(
        read_tsv(
            path
        )
    )


comparison_lookup = {}


for row in comparison_rows:

    if (
        row.get(
            "status",
            "",
        ).strip()
        != "PASS"
    ):
        continue


    query_species = row.get(
        "query_species",
        "",
    ).strip()


    subject_species = row.get(
        "subject_species",
        "",
    ).strip()


    if (
        not query_species
        or not subject_species
    ):
        continue


    anchor_value = row.get(
        "anchor_file",
        "",
    ).strip()


    if not anchor_value:
        continue


    anchor_path = Path(
        anchor_value
    )


    if not anchor_path.is_absolute():

        anchor_path = (
            PROJECT_ROOT
            / anchor_path
        ).resolve()


    if (
        not anchor_path.is_file()
        or anchor_path.stat().st_size == 0
    ):
        continue


    key = tuple(
        sorted(
            [
                query_species,
                subject_species,
            ]
        )
    )


    comparison_lookup[
        key
    ] = {
        "comparison_id":
            row.get(
                "comparison_id",
                f"{query_species}__{subject_species}",
            ),
        "query_species":
            query_species,
        "subject_species":
            subject_species,
        "anchor_file":
            anchor_path,
    }


def get_comparison(
    species_a,
    species_b,
):

    key = tuple(
        sorted(
            [
                species_a,
                species_b,
            ]
        )
    )


    if key not in comparison_lookup:

        available = "\n".join(
            "- " + "__".join(
                item
            )
            for item in sorted(
                comparison_lookup
            )
        )

        raise SystemExit(
            "ERROR: No successful full-anchor "
            f"comparison found for "
            f"{species_a}-{species_b}.\n"
            f"Available comparisons:\n"
            f"{available}"
        )


    return comparison_lookup[
        key
    ]


for species_a, species_b in [
    ("VSCU", "VSER"),
    ("VSCU", "VPAN"),
    ("VPAN", "VPER"),
]:

    get_comparison(
        species_a,
        species_b,
    )


###############################################################################
# Read and orient anchor pairs
###############################################################################

def read_anchor_pairs(
    path: Path,
):

    pairs = []


    with path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:

        for line in handle:

            stripped = line.strip()


            if (
                not stripped
                or stripped.startswith("#")
            ):
                continue


            fields = stripped.split()


            if len(fields) < 2:
                continue


            pairs.append(
                (
                    fields[0],
                    fields[1],
                )
            )


    return pairs


def map_anchor_pairs(
    reference,
    target,
    comparison,
):

    reference_genes = (
        species_gene_data[
            reference
        ][0]
    )

    target_genes = (
        species_gene_data[
            target
        ][0]
    )


    stored_query = (
        comparison[
            "query_species"
        ]
    )

    stored_subject = (
        comparison[
            "subject_species"
        ]
    )


    raw_pairs = read_anchor_pairs(
        comparison[
            "anchor_file"
        ]
    )


    mapped_pairs = []

    unmapped_reference = 0
    unmapped_target = 0


    for (
        query_gene,
        subject_gene,
    ) in raw_pairs:


        if (
            stored_query == reference
            and stored_subject == target
        ):

            reference_gene = (
                query_gene
            )

            target_gene = (
                subject_gene
            )


        elif (
            stored_query == target
            and stored_subject == reference
        ):

            reference_gene = (
                subject_gene
            )

            target_gene = (
                query_gene
            )


        else:

            raise SystemExit(
                "ERROR: Comparison orientation "
                f"cannot be resolved for "
                f"{reference}-{target}."
            )


        reference_present = (
            reference_gene
            in reference_genes
        )

        target_present = (
            target_gene
            in target_genes
        )


        if not reference_present:
            unmapped_reference += 1


        if not target_present:
            unmapped_target += 1


        if (
            not reference_present
            or not target_present
        ):
            continue


        mapped_pairs.append(
            {
                "reference_gene":
                    reference_gene,
                "target_gene":
                    target_gene,
                "reference_chromosome":
                    reference_genes[
                        reference_gene
                    ][
                        "chromosome"
                    ],
                "target_chromosome":
                    target_genes[
                        target_gene
                    ][
                        "chromosome"
                    ],
                "reference_midpoint":
                    reference_genes[
                        reference_gene
                    ][
                        "midpoint"
                    ],
                "target_midpoint":
                    target_genes[
                        target_gene
                    ][
                        "midpoint"
                    ],
            }
        )


    if not raw_pairs:

        raise SystemExit(
            f"ERROR: No anchor pairs for "
            f"{reference}-{target}."
        )


    if not mapped_pairs:

        raise SystemExit(
            f"ERROR: No mapped anchor pairs for "
            f"{reference}-{target}."
        )


    mapping_fraction = (
        len(
            mapped_pairs
        )
        /
        len(
            raw_pairs
        )
    )


    if mapping_fraction < 0.80:

        raise SystemExit(
            f"ERROR: Mapping fraction below 0.80 "
            f"for {reference}-{target}: "
            f"{mapping_fraction:.6f}"
        )


    return {
        "raw_anchor_count":
            len(
                raw_pairs
            ),
        "mapped_pairs":
            mapped_pairs,
        "mapped_anchor_count":
            len(
                mapped_pairs
            ),
        "unmapped_reference":
            unmapped_reference,
        "unmapped_target":
            unmapped_target,
        "mapping_fraction":
            mapping_fraction,
    }


###############################################################################
# Chromosome anchor matrix
###############################################################################

def calculate_anchor_matrix(
    mapped_pairs,
):

    counts = Counter()


    for pair in mapped_pairs:

        counts[
            (
                pair[
                    "reference_chromosome"
                ],
                pair[
                    "target_chromosome"
                ],
            )
        ] += 1


    return counts


###############################################################################
# Strong chromosome relationships
###############################################################################

def calculate_strong_relationships(
    matrix,
):

    reference_partner_counts = defaultdict(
        Counter
    )

    target_partner_counts = defaultdict(
        Counter
    )


    for (
        reference_chromosome,
        target_chromosome,
    ), count in matrix.items():

        reference_partner_counts[
            reference_chromosome
        ][
            target_chromosome
        ] = count


        target_partner_counts[
            target_chromosome
        ][
            reference_chromosome
        ] = count


    def strong_partners(
        partner_counts,
    ):

        if not partner_counts:
            return set()


        strongest_count = max(
            partner_counts.values()
        )


        return {
            partner
            for (
                partner,
                count,
            ) in partner_counts.items()
            if (
                count >= MIN_STRONG_ANCHORS
                and
                count
                >=
                RELATIVE_STRONG_THRESHOLD
                * strongest_count
            )
        }


    reference_strong = {
        chromosome:
            strong_partners(
                partner_counts
            )
        for (
            chromosome,
            partner_counts,
        ) in reference_partner_counts.items()
    }


    target_strong = {
        chromosome:
            strong_partners(
                partner_counts
            )
        for (
            chromosome,
            partner_counts,
        ) in target_partner_counts.items()
    }


    rearrangement_pairs = set()


    for (
        chromosome,
        partners,
    ) in reference_strong.items():

        if len(
            partners
        ) > 1:

            for partner in partners:

                rearrangement_pairs.add(
                    (
                        chromosome,
                        partner,
                    )
                )


    for (
        chromosome,
        partners,
    ) in target_strong.items():

        if len(
            partners
        ) > 1:

            for partner in partners:

                rearrangement_pairs.add(
                    (
                        partner,
                        chromosome,
                    )
                )


    return {
        "reference_partner_counts":
            reference_partner_counts,
        "target_partner_counts":
            target_partner_counts,
        "reference_strong":
            reference_strong,
        "target_strong":
            target_strong,
        "rearrangement_pairs":
            rearrangement_pairs,
    }


###############################################################################
# Operational chromosome-component ratio
###############################################################################

def operational_component_ratio(
    matrix,
    reference_strong,
):

    graph = defaultdict(
        set
    )


    for (
        reference_chromosome,
        partners,
    ) in reference_strong.items():

        for target_chromosome in partners:

            reference_node = (
                "R",
                reference_chromosome,
            )

            target_node = (
                "T",
                target_chromosome,
            )


            graph[
                reference_node
            ].add(
                target_node
            )


            graph[
                target_node
            ].add(
                reference_node
            )


    visited = set()

    components = []


    for start_node in graph:

        if start_node in visited:
            continue


        queue = deque(
            [
                start_node
            ]
        )

        visited.add(
            start_node
        )


        reference_set = set()
        target_set = set()


        while queue:

            side, chromosome = (
                queue.popleft()
            )


            if side == "R":

                reference_set.add(
                    chromosome
                )

            else:

                target_set.add(
                    chromosome
                )


            for neighbour in graph[
                (
                    side,
                    chromosome,
                )
            ]:

                if neighbour not in visited:

                    visited.add(
                        neighbour
                    )

                    queue.append(
                        neighbour
                    )


        support = sum(
            matrix.get(
                (
                    reference_chromosome,
                    target_chromosome,
                ),
                0,
            )
            for reference_chromosome
            in reference_set
            for target_chromosome
            in target_set
        )


        components.append(
            {
                "reference_count":
                    len(
                        reference_set
                    ),
                "target_count":
                    len(
                        target_set
                    ),
                "support":
                    support,
                "reference_chromosomes":
                    ",".join(
                        sorted(
                            reference_set,
                            key=numeric_sort_key,
                        )
                    ),
                "target_chromosomes":
                    ",".join(
                        sorted(
                            target_set,
                            key=numeric_sort_key,
                        )
                    ),
            }
        )


    if not components:

        return (
            "NA",
            [],
        )


    components.sort(
        key=lambda component: (
            -component[
                "support"
            ],
            -component[
                "reference_count"
            ],
            -component[
                "target_count"
            ],
        )
    )


    dominant = (
        components[
            0
        ]
    )


    ratio = (
        f"{dominant['reference_count']}:"
        f"{dominant['target_count']}"
    )


    return (
        ratio,
        components,
    )


###############################################################################
# Optimize target chromosome order
###############################################################################

def optimize_target_order(
    target,
    reference_order,
    matrix,
):

    reference_rank = {
        chromosome:
            index
        for (
            index,
            chromosome,
        ) in enumerate(
            reference_order,
            start=1,
        )
    }


    target_rows = []


    for target_chromosome in (
        numeric_chromosomes(
            target
        )
    ):

        partner_counts = {
            reference_chromosome:
                matrix.get(
                    (
                        reference_chromosome,
                        target_chromosome,
                    ),
                    0,
                )
            for reference_chromosome
            in reference_order
        }


        total_anchors = sum(
            partner_counts.values()
        )


        strongest_reference = max(
            reference_order,
            key=lambda chromosome: (
                partner_counts[
                    chromosome
                ],
                -reference_rank[
                    chromosome
                ],
            ),
        )


        strongest_count = (
            partner_counts[
                strongest_reference
            ]
        )


        if total_anchors > 0:

            weighted_reference_position = (
                sum(
                    reference_rank[
                        chromosome
                    ]
                    * count
                    for (
                        chromosome,
                        count,
                    ) in partner_counts.items()
                )
                /
                total_anchors
            )

        else:

            weighted_reference_position = (
                len(
                    reference_order
                )
                +
                safe_int(
                    target_chromosome,
                    999,
                )
            )


        target_rows.append(
            {
                "target_chromosome":
                    target_chromosome,
                "strongest_reference_chromosome":
                    strongest_reference,
                "strongest_anchor_count":
                    strongest_count,
                "total_reference_anchor_count":
                    total_anchors,
                "weighted_reference_position":
                    weighted_reference_position,
                "strongest_reference_rank":
                    reference_rank[
                        strongest_reference
                    ],
            }
        )


    target_rows.sort(
        key=lambda row: (
            row[
                "strongest_reference_rank"
            ],
            row[
                "weighted_reference_position"
            ],
            -row[
                "strongest_anchor_count"
            ],
            numeric_sort_key(
                row[
                    "target_chromosome"
                ]
            ),
        )
    )


    target_order = [
        row[
            "target_chromosome"
        ]
        for row in target_rows
    ]


    for (
        plot_rank,
        row,
    ) in enumerate(
        target_rows,
        start=1,
    ):

        row[
            "optimized_plot_rank"
        ] = plot_rank


    return (
        target_order,
        target_rows,
    )


###############################################################################
# Infer target chromosome orientation
###############################################################################

def infer_target_orientation(
    reference,
    target,
    mapped_pairs,
    reference_order,
):

    reference_rank = {
        chromosome:
            index
        for (
            index,
            chromosome,
        ) in enumerate(
            reference_order,
            start=1,
        )
    }


    by_target_chromosome = defaultdict(
        list
    )


    for pair in mapped_pairs:

        by_target_chromosome[
            pair[
                "target_chromosome"
            ]
        ].append(
            pair
        )


    orientations = {}


    for target_chromosome in (
        numeric_chromosomes(
            target
        )
    ):

        pairs = by_target_chromosome.get(
            target_chromosome,
            [],
        )


        if len(
            pairs
        ) < 5:

            orientations[
                target_chromosome
            ] = {
                "orientation":
                    "+",
                "correlation":
                    0.0,
                "anchors_used":
                    len(
                        pairs
                    ),
            }

            continue


        target_length = (
            chromosome_lengths[
                target
            ][
                target_chromosome
            ]
        )


        x_values = []
        y_values = []


        for pair in pairs:

            reference_chromosome = (
                pair[
                    "reference_chromosome"
                ]
            )


            reference_length = (
                chromosome_lengths[
                    reference
                ][
                    reference_chromosome
                ]
            )


            reference_global = (
                reference_rank[
                    reference_chromosome
                ]
                +
                pair[
                    "reference_midpoint"
                ]
                /
                max(
                    1,
                    reference_length,
                )
            )


            target_relative = (
                pair[
                    "target_midpoint"
                ]
                /
                max(
                    1,
                    target_length,
                )
            )


            x_values.append(
                reference_global
            )

            y_values.append(
                target_relative
            )


        x_mean = (
            sum(
                x_values
            )
            /
            len(
                x_values
            )
        )


        y_mean = (
            sum(
                y_values
            )
            /
            len(
                y_values
            )
        )


        numerator = sum(
            (
                x
                -
                x_mean
            )
            *
            (
                y
                -
                y_mean
            )
            for (
                x,
                y,
            ) in zip(
                x_values,
                y_values,
            )
        )


        x_denominator = math.sqrt(
            sum(
                (
                    x
                    -
                    x_mean
                )
                ** 2
                for x in x_values
            )
        )


        y_denominator = math.sqrt(
            sum(
                (
                    y
                    -
                    y_mean
                )
                ** 2
                for y in y_values
            )
        )


        denominator = (
            x_denominator
            *
            y_denominator
        )


        correlation = (
            numerator
            /
            denominator
            if denominator > 0
            else 0.0
        )


        orientations[
            target_chromosome
        ] = {
            "orientation":
                (
                    "-"
                    if correlation < 0
                    else "+"
                ),
            "correlation":
                correlation,
            "anchors_used":
                len(
                    pairs
                ),
        }


    return orientations


###############################################################################
# Matrix tables
###############################################################################

def create_matrix_rows(
    comparison_name,
    reference,
    target,
    reference_order,
    target_order,
    matrix,
):

    rows = []


    for reference_chromosome in (
        reference_order
    ):

        row = {
            "comparison":
                comparison_name,
            "reference_species":
                reference,
            "reference_chromosome":
                reference_chromosome,
        }


        for target_chromosome in (
            target_order
        ):

            row[
                f"{target}_chr{target_chromosome}"
            ] = matrix.get(
                (
                    reference_chromosome,
                    target_chromosome,
                ),
                0,
            )


        rows.append(
            row
        )


    return rows


def create_long_matrix_rows(
    comparison_name,
    reference,
    target,
    reference_order,
    target_order,
    matrix,
    rearrangement_pairs,
):

    rows = []


    for reference_chromosome in (
        reference_order
    ):

        for target_chromosome in (
            target_order
        ):

            anchor_count = matrix.get(
                (
                    reference_chromosome,
                    target_chromosome,
                ),
                0,
            )


            rows.append(
                {
                    "comparison":
                        comparison_name,
                    "reference_species":
                        reference,
                    "reference_chromosome":
                        reference_chromosome,
                    "target_species":
                        target,
                    "target_chromosome":
                        target_chromosome,
                    "anchor_count":
                        anchor_count,
                    "multiple_strong_partner_relationship":
                        (
                            "YES"
                            if (
                                reference_chromosome,
                                target_chromosome,
                            )
                            in rearrangement_pairs
                            else "NO"
                        ),
                }
            )


    return rows


###############################################################################
# Plot layout
###############################################################################

def build_track_layout(
    species,
    chromosome_order,
):

    total_length = sum(
        chromosome_lengths[
            species
        ][
            chromosome
        ]
        for chromosome
        in chromosome_order
    )


    usable_width = (
        1.0
        -
        CHROMOSOME_GAP
        *
        max(
            0,
            len(
                chromosome_order
            )
            -
            1,
        )
    )


    if usable_width <= 0:

        raise SystemExit(
            f"ERROR: Chromosome gap too large "
            f"for {species}."
        )


    cursor = 0.0

    layout = {}


    for chromosome in (
        chromosome_order
    ):

        chromosome_width = (
            usable_width
            *
            chromosome_lengths[
                species
            ][
                chromosome
            ]
            /
            total_length
        )


        layout[
            chromosome
        ] = {
            "start":
                cursor,
            "end":
                cursor
                +
                chromosome_width,
            "midpoint":
                cursor
                +
                chromosome_width
                /
                2.0,
            "length":
                chromosome_lengths[
                    species
                ][
                    chromosome
                ],
        }


        cursor += (
            chromosome_width
            +
            CHROMOSOME_GAP
        )


    return layout


def deterministic_subsample(
    records,
    maximum_records,
):

    if len(
        records
    ) <= maximum_records:

        return records


    step = (
        len(
            records
        )
        /
        maximum_records
    )


    selected = []

    position = 0.0


    while (
        len(
            selected
        )
        <
        maximum_records
        and
        int(
            position
        )
        <
        len(
            records
        )
    ):

        selected.append(
            records[
                int(
                    position
                )
            ]
        )

        position += step


    return selected


def anchor_x(
    chromosome,
    midpoint,
    track_layout,
    orientation="+",
):

    chromosome_layout = (
        track_layout[
            chromosome
        ]
    )


    fraction = (
        midpoint
        /
        max(
            1,
            chromosome_layout[
                "length"
            ],
        )
    )


    fraction = min(
        1.0,
        max(
            0.0,
            fraction,
        ),
    )


    if orientation == "-":

        fraction = (
            1.0
            -
            fraction
        )


    return (
        chromosome_layout[
            "start"
        ]
        +
        fraction
        *
        (
            chromosome_layout[
                "end"
            ]
            -
            chromosome_layout[
                "start"
            ]
        )
    )


###############################################################################
# Pairwise plot
###############################################################################

def draw_pairwise_plot(
    comparison_name,
    purpose,
    reference,
    target,
    reference_order,
    target_order,
    target_orientations,
    mapped_pairs,
    rearrangement_pairs,
    operational_ratio,
):

    reference_layout = build_track_layout(
        reference,
        reference_order,
    )


    target_layout = build_track_layout(
        target,
        target_order,
    )


    figure, axis = plt.subplots(
        figsize=(
            18,
            7.2,
        )
    )


    reference_y = 1.0
    target_y = 0.0


    axis.set_xlim(
        -0.14,
        1.075,
    )


    axis.set_ylim(
        -0.38,
        1.38,
    )


    axis.axis(
        "off"
    )


    plot_pairs = deterministic_subsample(
        mapped_pairs,
        MAX_LINKS_PER_PAIRWISE_PLOT,
    )


    normal_pairs = [
        pair
        for pair
        in plot_pairs
        if (
            pair[
                "reference_chromosome"
            ],
            pair[
                "target_chromosome"
            ],
        )
        not in rearrangement_pairs
    ]


    rearranged_pairs = [
        pair
        for pair
        in plot_pairs
        if (
            pair[
                "reference_chromosome"
            ],
            pair[
                "target_chromosome"
            ],
        )
        in rearrangement_pairs
    ]


    for pair in (
        normal_pairs
        +
        rearranged_pairs
    ):

        reference_chromosome = (
            pair[
                "reference_chromosome"
            ]
        )


        target_chromosome = (
            pair[
                "target_chromosome"
            ]
        )


        reference_x = anchor_x(
            reference_chromosome,
            pair[
                "reference_midpoint"
            ],
            reference_layout,
            orientation="+",
        )


        target_orientation = (
            target_orientations[
                target_chromosome
            ][
                "orientation"
            ]
        )


        target_x = anchor_x(
            target_chromosome,
            pair[
                "target_midpoint"
            ],
            target_layout,
            orientation=target_orientation,
        )


        is_rearrangement = (
            (
                reference_chromosome,
                target_chromosome,
            )
            in rearrangement_pairs
        )


        if is_rearrangement:

            edge_colour = (
                RED_LINK
            )

            alpha = 0.12
            linewidth = 0.42
            zorder = 2

        else:

            edge_colour = (
                BLUE_LINK
            )

            alpha = 0.045
            linewidth = 0.24
            zorder = 1


        vertices = [
            (
                reference_x,
                reference_y
                -
                CHROMOSOME_HEIGHT
                /
                2.0,
            ),
            (
                reference_x,
                0.68,
            ),
            (
                target_x,
                0.32,
            ),
            (
                target_x,
                target_y
                +
                CHROMOSOME_HEIGHT
                /
                2.0,
            ),
        ]


        codes = [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
        ]


        axis.add_patch(
            PathPatch(
                MplPath(
                    vertices,
                    codes,
                ),
                facecolor="none",
                edgecolor=edge_colour,
                linewidth=linewidth,
                alpha=alpha,
                zorder=zorder,
            )
        )


    def draw_species_track(
        species,
        chromosome_order,
        track_layout,
        y,
        orientations=None,
    ):

        colour = (
            SPECIES_COLOURS[
                species
            ]
        )


        axis.text(
            -0.025,
            y + 0.02,
            (
                f"{PLOIDY[species]} "
                f"{species}"
            ),
            fontsize=15,
            fontweight="bold",
            ha="right",
            va="center",
            color=colour,
        )


        axis.text(
            -0.025,
            y - 0.11,
            SPECIES_NAMES[
                species
            ],
            fontsize=9.5,
            fontstyle="italic",
            ha="right",
            va="center",
            color="#444444",
        )


        for chromosome in (
            chromosome_order
        ):

            chromosome_layout = (
                track_layout[
                    chromosome
                ]
            )


            patch = FancyBboxPatch(
                (
                    chromosome_layout[
                        "start"
                    ],
                    y
                    -
                    CHROMOSOME_HEIGHT
                    /
                    2.0,
                ),
                (
                    chromosome_layout[
                        "end"
                    ]
                    -
                    chromosome_layout[
                        "start"
                    ]
                ),
                CHROMOSOME_HEIGHT,
                boxstyle=(
                    "round,pad=0.001,"
                    "rounding_size=0.010"
                ),
                facecolor=colour,
                edgecolor="#222222",
                linewidth=0.8,
                zorder=4,
            )


            axis.add_patch(
                patch
            )


            orientation_text = ""


            if orientations is not None:

                orientation = (
                    orientations[
                        chromosome
                    ][
                        "orientation"
                    ]
                )


                orientation_text = (
                    " <"
                    if orientation == "-"
                    else " >"
                )


            # IMPORTANT:
            # chromosome labels are NUMBERS ONLY.
            # No "Chr" prefix is introduced.
            axis.text(
                chromosome_layout[
                    "midpoint"
                ],
                y + 0.105,
                (
                    f"{chromosome}"
                    f"{orientation_text}"
                ),
                fontsize=10,
                fontweight="bold",
                ha="center",
                va="bottom",
                color="#222222",
                zorder=5,
            )


    draw_species_track(
        reference,
        reference_order,
        reference_layout,
        reference_y,
        orientations=None,
    )


    draw_species_track(
        target,
        target_order,
        target_layout,
        target_y,
        orientations=target_orientations,
    )


    TITLE_PURPOSE = {
        "diploid_diploid_VSER":
            "Diploid-diploid comparison",
        "diploid_diploid_VPAN":
            "Diploid-diploid comparison",
        "diploid_tetraploid_VPER":
            "Diploid-tetraploid comparison",
    }


    axis.set_title(
        (
            f"{reference} versus {target}: "
            "synteny-optimized chromosome order\n"
            f"{TITLE_PURPOSE[purpose]}"
        ),
        fontsize=16,
        fontweight="bold",
        pad=20,
    )


    axis.text(
        1.015,
        0.50,
        (
            "Operational\n"
            "component ratio\n"
            f"{operational_ratio}"
        ),
        fontsize=9,
        ha="left",
        va="center",
        color="#333333",
    )


    figure.text(
        0.015,
        0.025,
        (
            "Blue links: mapped full-anchor synteny pairs. "
            "Red links: chromosome pairs participating in "
            "multiple strong-partner relationships. "
            "Target chromosomes are ordered by syntenic "
            "similarity to the reference. The < and > symbols "
            "indicate plotting orientation inferred from "
            "anchor-position correlation."
        ),
        fontsize=8.4,
        ha="left",
        va="bottom",
        color="#333333",
    )


    figure.tight_layout(
        rect=(
            0.025,
            0.075,
            0.93,
            0.96,
        )
    )


    output_prefix = (
        PLOT_DIR
        /
        (
            f"{comparison_name}."
            "synteny_optimized"
        )
    )


    figure.savefig(
        Path(
            str(
                output_prefix
            )
            +
            ".pdf"
        ),
        bbox_inches="tight",
    )


    figure.savefig(
        Path(
            str(
                output_prefix
            )
            +
            ".svg"
        ),
        bbox_inches="tight",
    )


    figure.savefig(
        Path(
            str(
                output_prefix
            )
            +
            ".png"
        ),
        dpi=400,
        bbox_inches="tight",
    )


    plt.close(
        figure
    )


###############################################################################
# Run comparisons in dependency order
###############################################################################

all_qc_rows = []

all_summary_rows = []

all_order_rows = []

all_strong_partner_rows = []

all_component_rows = []


# This will be generated by VSCU_VPAN
# and then inherited by VPAN_VPER.
vpan_vscu_informed_order = None


for comparison_definition in (
    PAIRWISE_COMPARISONS
):

    comparison_name = (
        comparison_definition[
            "comparison_id"
        ]
    )


    reference = (
        comparison_definition[
            "reference"
        ]
    )


    target = (
        comparison_definition[
            "target"
        ]
    )


    purpose = (
        comparison_definition[
            "purpose"
        ]
    )


    ###########################################################################
    # Reference-order dependency
    ###########################################################################

    if reference == "VSCU":

        reference_order = list(
            FROZEN_VSCU_ORDER
        )


    elif comparison_name == "VPAN_VPER":

        if (
            vpan_vscu_informed_order
            is None
        ):

            raise SystemExit(
                "ERROR: VSCU-informed VPAN order "
                "was not generated before VPAN_VPER."
            )


        reference_order = list(
            vpan_vscu_informed_order
        )


    else:

        raise SystemExit(
            f"ERROR: Unexpected reference-order "
            f"dependency for {comparison_name}."
        )


    comparison = get_comparison(
        reference,
        target,
    )


    mapping_result = map_anchor_pairs(
        reference,
        target,
        comparison,
    )


    mapped_pairs = (
        mapping_result[
            "mapped_pairs"
        ]
    )


    matrix = (
        calculate_anchor_matrix(
            mapped_pairs
        )
    )


    strong_results = (
        calculate_strong_relationships(
            matrix
        )
    )


    target_order, target_order_rows = (
        optimize_target_order(
            target,
            reference_order,
            matrix,
        )
    )


    # ---------------------------------------------------------
    # Freeze VPAN order from VSCU -> VPAN.
    # This exact order becomes the reference for VPAN -> VPER.
    # ---------------------------------------------------------

    if comparison_name == "VSCU_VPAN":

        vpan_vscu_informed_order = list(
            target_order
        )


    target_orientations = (
        infer_target_orientation(
            reference,
            target,
            mapped_pairs,
            reference_order,
        )
    )


    component_ratio, components = (
        operational_component_ratio(
            matrix,
            strong_results[
                "reference_strong"
            ],
        )
    )


    ###########################################################################
    # Wide chromosome-anchor matrix
    ###########################################################################

    matrix_rows = create_matrix_rows(
        comparison_name,
        reference,
        target,
        reference_order,
        target_order,
        matrix,
    )


    matrix_fieldnames = [
        "comparison",
        "reference_species",
        "reference_chromosome",
    ] + [
        (
            f"{target}_"
            f"chr{chromosome}"
        )
        for chromosome
        in target_order
    ]


    write_tsv(
        TABLE_DIR
        /
        (
            f"{comparison_name}."
            "chromosome_anchor_matrix.tsv"
        ),
        matrix_rows,
        matrix_fieldnames,
    )


    ###########################################################################
    # Long chromosome-anchor matrix
    ###########################################################################

    long_matrix_rows = (
        create_long_matrix_rows(
            comparison_name,
            reference,
            target,
            reference_order,
            target_order,
            matrix,
            strong_results[
                "rearrangement_pairs"
            ],
        )
    )


    write_tsv(
        TABLE_DIR
        /
        (
            f"{comparison_name}."
            "chromosome_anchor_matrix.long.tsv"
        ),
        long_matrix_rows,
        [
            "comparison",
            "reference_species",
            "reference_chromosome",
            "target_species",
            "target_chromosome",
            "anchor_count",
            "multiple_strong_partner_relationship",
        ],
    )


    ###########################################################################
    # Optimized chromosome-order table
    ###########################################################################

    reference_order_rows = []


    for (
        rank,
        chromosome,
    ) in enumerate(
        reference_order,
        start=1,
    ):

        reference_order_rows.append(
            {
                "comparison":
                    comparison_name,
                "species_role":
                    "reference",
                "species_code":
                    reference,
                "chromosome":
                    chromosome,
                "plot_rank":
                    rank,
                "orientation":
                    "+",
                "orientation_correlation":
                    "",
                "orientation_anchors":
                    "",
                "strongest_partner":
                    "",
                "strongest_anchor_count":
                    "",
                "weighted_partner_position":
                    "",
                "order_source":
                    comparison_definition[
                        "reference_order_source"
                    ],
            }
        )


    target_rows_lookup = {
        row[
            "target_chromosome"
        ]:
        row
        for row
        in target_order_rows
    }


    target_output_rows = []


    for (
        rank,
        chromosome,
    ) in enumerate(
        target_order,
        start=1,
    ):

        order_row = (
            target_rows_lookup[
                chromosome
            ]
        )


        orientation_row = (
            target_orientations[
                chromosome
            ]
        )


        target_output_rows.append(
            {
                "comparison":
                    comparison_name,
                "species_role":
                    "target",
                "species_code":
                    target,
                "chromosome":
                    chromosome,
                "plot_rank":
                    rank,
                "orientation":
                    orientation_row[
                        "orientation"
                    ],
                "orientation_correlation":
                    (
                        f"{orientation_row['correlation']:.6f}"
                    ),
                "orientation_anchors":
                    orientation_row[
                        "anchors_used"
                    ],
                "strongest_partner":
                    order_row[
                        "strongest_reference_chromosome"
                    ],
                "strongest_anchor_count":
                    order_row[
                        "strongest_anchor_count"
                    ],
                "weighted_partner_position":
                    (
                        f"{order_row['weighted_reference_position']:.6f}"
                    ),
                "order_source":
                    "optimized_by_full_anchor_similarity",
            }
        )


    combined_order_rows = (
        reference_order_rows
        +
        target_output_rows
    )


    write_tsv(
        ORDER_DIR
        /
        (
            f"{comparison_name}."
            "optimized_chromosome_order.tsv"
        ),
        combined_order_rows,
        [
            "comparison",
            "species_role",
            "species_code",
            "chromosome",
            "plot_rank",
            "orientation",
            "orientation_correlation",
            "orientation_anchors",
            "strongest_partner",
            "strongest_anchor_count",
            "weighted_partner_position",
            "order_source",
        ],
    )


    all_order_rows.extend(
        combined_order_rows
    )


    ###########################################################################
    # Strong-partner summary — reference and target species
    ###########################################################################

    for focal_role in [
        "reference",
        "target",
    ]:

        if focal_role == "reference":

            focal_species = reference

            partner_species = target

            focal_order = (
                reference_order
            )

            partner_counts_lookup = (
                strong_results[
                    "reference_partner_counts"
                ]
            )

            strong_lookup = (
                strong_results[
                    "reference_strong"
                ]
            )


        else:

            focal_species = target

            partner_species = reference

            focal_order = (
                target_order
            )

            partner_counts_lookup = (
                strong_results[
                    "target_partner_counts"
                ]
            )

            strong_lookup = (
                strong_results[
                    "target_strong"
                ]
            )


        for focal_chromosome in (
            focal_order
        ):

            partner_counts = (
                partner_counts_lookup.get(
                    focal_chromosome,
                    {},
                )
            )


            strong_partners = sorted(
                strong_lookup.get(
                    focal_chromosome,
                    set(),
                ),
                key=numeric_sort_key,
            )


            strongest_partner = ""


            if partner_counts:

                strongest_partner = max(
                    partner_counts,
                    key=partner_counts.get,
                )


            all_strong_partner_rows.append(
                {
                    "comparison":
                        comparison_name,
                    "focal_species":
                        focal_species,
                    "focal_chromosome":
                        focal_chromosome,
                    "partner_species":
                        partner_species,
                    "strongest_partner":
                        strongest_partner,
                    "strongest_anchor_count":
                        (
                            partner_counts.get(
                                strongest_partner,
                                0,
                            )
                            if strongest_partner
                            else 0
                        ),
                    "strong_partner_count":
                        len(
                            strong_partners
                        ),
                    "strong_partners":
                        ",".join(
                            strong_partners
                        ),
                    "relationship":
                        (
                            "MULTIPLE_STRONG_PARTNERS"
                            if len(
                                strong_partners
                            ) > 1
                            else (
                                "ONE_STRONG_PARTNER"
                                if len(
                                    strong_partners
                                ) == 1
                                else "NO_STRONG_PARTNER"
                            )
                        ),
                }
            )


    ###########################################################################
    # Operational components
    ###########################################################################

    for (
        component_rank,
        component,
    ) in enumerate(
        components,
        start=1,
    ):

        all_component_rows.append(
            {
                "comparison":
                    comparison_name,
                "component_rank":
                    component_rank,
                "reference_species":
                    reference,
                "target_species":
                    target,
                "reference_chromosome_count":
                    component[
                        "reference_count"
                    ],
                "target_chromosome_count":
                    component[
                        "target_count"
                    ],
                "component_ratio":
                    (
                        f"{component['reference_count']}:"
                        f"{component['target_count']}"
                    ),
                "anchor_support":
                    component[
                        "support"
                    ],
                "reference_chromosomes":
                    component[
                        "reference_chromosomes"
                    ],
                "target_chromosomes":
                    component[
                        "target_chromosomes"
                    ],
                "dominant_component":
                    (
                        "YES"
                        if component_rank == 1
                        else "NO"
                    ),
            }
        )


    ###########################################################################
    # Mapping QC
    ###########################################################################

    all_qc_rows.append(
        {
            "comparison":
                comparison_name,
            "reference_species":
                reference,
            "target_species":
                target,
            "anchor_file":
                str(
                    comparison[
                        "anchor_file"
                    ]
                ),
            "raw_anchor_pairs":
                mapping_result[
                    "raw_anchor_count"
                ],
            "mapped_anchor_pairs":
                mapping_result[
                    "mapped_anchor_count"
                ],
            "unmapped_reference_genes":
                mapping_result[
                    "unmapped_reference"
                ],
            "unmapped_target_genes":
                mapping_result[
                    "unmapped_target"
                ],
            "mapping_fraction":
                (
                    f"{mapping_result['mapping_fraction']:.6f}"
                ),
            "status":
                "PASS",
        }
    )


    ###########################################################################
    # Comparison summary
    ###########################################################################

    all_summary_rows.append(
        {
            "comparison":
                comparison_name,
            "purpose":
                purpose,
            "reference_species":
                reference,
            "target_species":
                target,
            "reference_chromosome_count":
                len(
                    reference_order
                ),
            "target_chromosome_count":
                len(
                    target_order
                ),
            "mapped_anchor_pairs":
                len(
                    mapped_pairs
                ),
            "chromosome_pairs_with_anchors":
                len(
                    matrix
                ),
            "multiple_strong_partner_pairs":
                len(
                    strong_results[
                        "rearrangement_pairs"
                    ]
                ),
            "operational_component_ratio":
                component_ratio,
            "reference_order":
                ",".join(
                    reference_order
                ),
            "target_optimized_order":
                ",".join(
                    target_order
                ),
            "status":
                "PASS",
        }
    )


    ###########################################################################
    # Generate figure
    ###########################################################################

    draw_pairwise_plot(
        comparison_name,
        purpose,
        reference,
        target,
        reference_order,
        target_order,
        target_orientations,
        mapped_pairs,
        strong_results[
            "rearrangement_pairs"
        ],
        component_ratio,
    )


    print(
        f"{comparison_name}: "
        f"mapped={len(mapped_pairs)}, "
        f"ratio={component_ratio}, "
        f"reference_order="
        f"{','.join(reference_order)}, "
        f"target_order="
        f"{','.join(target_order)}"
    )


###############################################################################
# Write combined tables
###############################################################################

write_tsv(
    TABLE_DIR
    /
    "step35Y2_mapping_qc.tsv",
    all_qc_rows,
    [
        "comparison",
        "reference_species",
        "target_species",
        "anchor_file",
        "raw_anchor_pairs",
        "mapped_anchor_pairs",
        "unmapped_reference_genes",
        "unmapped_target_genes",
        "mapping_fraction",
        "status",
    ],
)


write_tsv(
    TABLE_DIR
    /
    "step35Y2_pairwise_summary.tsv",
    all_summary_rows,
    [
        "comparison",
        "purpose",
        "reference_species",
        "target_species",
        "reference_chromosome_count",
        "target_chromosome_count",
        "mapped_anchor_pairs",
        "chromosome_pairs_with_anchors",
        "multiple_strong_partner_pairs",
        "operational_component_ratio",
        "reference_order",
        "target_optimized_order",
        "status",
    ],
)


write_tsv(
    ORDER_DIR
    /
    "step35Y2_all_optimized_orders.tsv",
    all_order_rows,
    [
        "comparison",
        "species_role",
        "species_code",
        "chromosome",
        "plot_rank",
        "orientation",
        "orientation_correlation",
        "orientation_anchors",
        "strongest_partner",
        "strongest_anchor_count",
        "weighted_partner_position",
        "order_source",
    ],
)


write_tsv(
    TABLE_DIR
    /
    "step35Y2_strong_partner_summary.tsv",
    all_strong_partner_rows,
    [
        "comparison",
        "focal_species",
        "focal_chromosome",
        "partner_species",
        "strongest_partner",
        "strongest_anchor_count",
        "strong_partner_count",
        "strong_partners",
        "relationship",
    ],
)


write_tsv(
    TABLE_DIR
    /
    "step35Y2_operational_components.tsv",
    all_component_rows,
    [
        "comparison",
        "component_rank",
        "reference_species",
        "target_species",
        "reference_chromosome_count",
        "target_chromosome_count",
        "component_ratio",
        "anchor_support",
        "reference_chromosomes",
        "target_chromosomes",
        "dominant_component",
    ],
)


###############################################################################
# Strong QC of inherited chromosome orders
###############################################################################

overall_status = (
    "PASS"
    if (
        len(
            all_summary_rows
        ) == 3

        and

        all(
            row[
                "status"
            ] == "PASS"
            for row
            in all_summary_rows
        )

        and

        all(
            float(
                row[
                    "mapping_fraction"
                ]
            ) >= 0.80
            for row
            in all_qc_rows
        )

        and

        all_summary_rows[
            0
        ][
            "reference_order"
        ]
        ==
        ",".join(
            FROZEN_VSCU_ORDER
        )

        and

        all_summary_rows[
            1
        ][
            "reference_order"
        ]
        ==
        ",".join(
            FROZEN_VSCU_ORDER
        )

        and

        all_summary_rows[
            2
        ][
            "reference_order"
        ]
        ==
        all_summary_rows[
            1
        ][
            "target_optimized_order"
        ]
    )

    else "FAIL"
)


overall_rows = [
    {
        "metric":
            "comparisons_expected",
        "value":
            3,
    },
    {
        "metric":
            "comparisons_completed",
        "value":
            len(
                all_summary_rows
            ),
    },
    {
        "metric":
            "comparisons_pass",
        "value":
            sum(
                row[
                    "status"
                ] == "PASS"
                for row
                in all_summary_rows
            ),
    },
    {
        "metric":
            "figures_expected",
        "value":
            9,
    },
    {
        "metric":
            "frozen_VSCU_order",
        "value":
            ",".join(
                FROZEN_VSCU_ORDER
            ),
    },
    {
        "metric":
            "VPAN_reference_order_source",
        "value":
            "VSCU_VPAN_target_optimized_order",
    },
    {
        "metric":
            "status",
        "value":
            overall_status,
    },
]


write_tsv(
    TABLE_DIR
    /
    "step35Y2_overall_summary.tsv",
    overall_rows,
    [
        "metric",
        "value",
    ],
)


if overall_status != "PASS":

    raise SystemExit(
        "ERROR: Step 35Y2 Python analysis failed QC."
    )


print(
    "Step 35Y2 Python analysis: PASS"
)

PY

###############################################################################
# Validate generated figures
###############################################################################

EXPECTED_PLOTS=(
    "${PLOT_DIR}/VSCU_VSER.synteny_optimized.pdf"
    "${PLOT_DIR}/VSCU_VSER.synteny_optimized.svg"
    "${PLOT_DIR}/VSCU_VSER.synteny_optimized.png"
    "${PLOT_DIR}/VSCU_VPAN.synteny_optimized.pdf"
    "${PLOT_DIR}/VSCU_VPAN.synteny_optimized.svg"
    "${PLOT_DIR}/VSCU_VPAN.synteny_optimized.png"
    "${PLOT_DIR}/VPAN_VPER.synteny_optimized.pdf"
    "${PLOT_DIR}/VPAN_VPER.synteny_optimized.svg"
    "${PLOT_DIR}/VPAN_VPER.synteny_optimized.png"
)

for FILE in "${EXPECTED_PLOTS[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Expected plot missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

###############################################################################
# Validate file signatures
###############################################################################

python - "${PLOT_DIR}" <<'PY'

from pathlib import Path
import sys


plot_dir = Path(
    sys.argv[1]
)


comparisons = [
    "VSCU_VSER",
    "VSCU_VPAN",
    "VPAN_VPER",
]


for comparison in comparisons:

    pdf = (
        plot_dir
        /
        f"{comparison}.synteny_optimized.pdf"
    )

    svg = (
        plot_dir
        /
        f"{comparison}.synteny_optimized.svg"
    )

    png = (
        plot_dir
        /
        f"{comparison}.synteny_optimized.png"
    )


    if (
        pdf.read_bytes()[
            :5
        ]
        != b"%PDF-"
    ):

        raise SystemExit(
            f"ERROR: Invalid PDF signature: {pdf}"
        )


    svg_head = (
        svg.read_text(
            encoding="utf-8",
            errors="ignore",
        )[
            :2000
        ].lower()
    )


    if "<svg" not in svg_head:

        raise SystemExit(
            f"ERROR: Invalid SVG signature: {svg}"
        )


    if (
        png.read_bytes()[
            :8
        ]
        !=
        b"\x89PNG\r\n\x1a\n"
    ):

        raise SystemExit(
            f"ERROR: Invalid PNG signature: {png}"
        )


print(
    "PASS: all nine figure files validated."
)

PY

###############################################################################
# Validate tables
###############################################################################

REQUIRED_TABLES=(
    "${TABLE_DIR}/VSCU_VSER.chromosome_anchor_matrix.tsv"
    "${TABLE_DIR}/VSCU_VSER.chromosome_anchor_matrix.long.tsv"
    "${TABLE_DIR}/VSCU_VPAN.chromosome_anchor_matrix.tsv"
    "${TABLE_DIR}/VSCU_VPAN.chromosome_anchor_matrix.long.tsv"
    "${TABLE_DIR}/VPAN_VPER.chromosome_anchor_matrix.tsv"
    "${TABLE_DIR}/VPAN_VPER.chromosome_anchor_matrix.long.tsv"
    "${TABLE_DIR}/step35Y2_mapping_qc.tsv"
    "${TABLE_DIR}/step35Y2_pairwise_summary.tsv"
    "${TABLE_DIR}/step35Y2_strong_partner_summary.tsv"
    "${TABLE_DIR}/step35Y2_operational_components.tsv"
    "${TABLE_DIR}/step35Y2_overall_summary.tsv"
    "${ORDER_DIR}/VSCU_VSER.optimized_chromosome_order.tsv"
    "${ORDER_DIR}/VSCU_VPAN.optimized_chromosome_order.tsv"
    "${ORDER_DIR}/VPAN_VPER.optimized_chromosome_order.tsv"
    "${ORDER_DIR}/step35Y2_all_optimized_orders.tsv"
)

for FILE in "${REQUIRED_TABLES[@]}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required Step 35Y2 output missing:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

if ! grep -Pq \
    '^status\tPASS$' \
    "${TABLE_DIR}/step35Y2_overall_summary.tsv"
then

    echo "ERROR: Step 35Y2 overall status is not PASS." >&2

    cat \
        "${TABLE_DIR}/step35Y2_overall_summary.tsv" \
        >&2

    exit 1
fi

###############################################################################
# Create checkpoint
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP35Y2_COMPLETE.txt" <<EOF2
checkpoint=step35Y2_additional_synteny_optimized_pairwise_plots
date=$(date --iso-8601=seconds)
comparisons=VSCU_VSER,VSCU_VPAN,VPAN_VPER
VSCU_VSER_role=diploid_diploid_comparison
VSCU_VPAN_role=diploid_diploid_comparison
VPAN_VPER_role=diploid_tetraploid_comparison
VSCU_reference_order=7,1,3,9,2,6,5,8,4
VSCU_reference_order_source=PMAJ_informed_step35Y_frozen
VPAN_reference_order_source=VSCU_VPAN_target_optimized_order
chromosome_labels=numbers_only
target_order=optimized_by_full_anchor_similarity
target_orientation=anchor_position_correlation
normal_links=blue_mapped_full_anchor_pairs
rearrangement_links=red_multiple_strong_partner_pairs
strong_partner_min_anchors=20
strong_partner_relative_threshold=0.25
max_links_per_plot=8000
formats=PDF,SVG,PNG
cleanup=targeted_step35Y2_outputs_only
existing_step34_outputs_modified=false
existing_step35Y_outputs_modified=false
existing_step35Y1_outputs_modified=false
status=PASS
next_step=inspect_additional_macrosynteny_plots_and_integrate_with_Ks_and_multisynteny
EOF2

###############################################################################
# Freeze important outputs
###############################################################################

cp -f \
    "${TABLE_DIR}/step35Y2_mapping_qc.tsv" \
    "${TABLE_DIR}/step35Y2_pairwise_summary.tsv" \
    "${TABLE_DIR}/step35Y2_strong_partner_summary.tsv" \
    "${TABLE_DIR}/step35Y2_operational_components.tsv" \
    "${ORDER_DIR}/step35Y2_all_optimized_orders.tsv" \
    "${CHECKPOINT_DIR}/"

###############################################################################
# Checksums
###############################################################################

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# Display final results
###############################################################################

echo
echo "============================================================"
echo "Step 35Y2 mapping QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/step35Y2_mapping_qc.tsv"

echo
echo "============================================================"
echo "Step 35Y2 pairwise summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/step35Y2_pairwise_summary.tsv"

echo
echo "============================================================"
echo "Step 35Y2 optimized chromosome orders"
echo "============================================================"

column -t -s $'\t' \
    "${ORDER_DIR}/step35Y2_all_optimized_orders.tsv"

echo
echo "============================================================"
echo "Step 35Y2 checkpoint"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP35Y2_COMPLETE.txt"

echo
echo "============================================================"
echo "Step 35Y2 completed successfully"
echo "============================================================"
