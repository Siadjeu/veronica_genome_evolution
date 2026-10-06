#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=macro34i
#SBATCH --output=10_synteny/logs/macro34i_%j.out
#SBATCH --error=10_synteny/logs/macro34i_%j.err

set -euo pipefail

export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

REFINED_DIR="${SYNTENY_DIR}/refined_macrosynteny"
TABLE_DIR="${REFINED_DIR}/tables"
PLOT_DIR="${REFINED_DIR}/publication_plots"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_publication_macrosynteny_plots"

mkdir -p \
    "${SYNTENY_DIR}/logs" \
    "${TABLE_DIR}" \
    "${PLOT_DIR}" \
    "${CHECKPOINT_DIR}"

cd "${PROJECT_DIR}"

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
else
    echo "ERROR: Conda could not be located." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate jcvi_env

export MPLBACKEND=Agg

echo "Python:"
python --version

rm -f "${CHECKPOINT_DIR}"/*

python - <<'PY'
from __future__ import annotations

import csv
import math
from collections import Counter, defaultdict, deque
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyBboxPatch, PathPatch


PROJECT_DIR = Path(
    os.environ["PROJECT_DIR"]
).resolve()

SYNTENY_DIR = PROJECT_DIR / "10_synteny"

REFINED_DIR = (
    SYNTENY_DIR
    / "refined_macrosynteny"
)

TABLE_DIR = (
    REFINED_DIR
    / "tables"
)

PLOT_DIR = (
    REFINED_DIR
    / "publication_plots"
)

ORDER_DIR = (
    REFINED_DIR
    / "orders"
)

JCVI_BED_DIR = (
    SYNTENY_DIR
    / "jcvi_inputs"
    / "bed"
)

INPUT_MANIFEST = (
    SYNTENY_DIR
    / "manifests"
    / "synteny_input_manifest.verified.tsv"
)

STEP34_RESULTS = (
    SYNTENY_DIR
    / "jcvi_step34"
    / "tables"
    / "jcvi_step34_results.tsv"
)

STEP34G_RESULTS = (
    REFINED_DIR
    / "missing_comparisons"
    / "tables"
    / "step34g_missing_comparison_results.tsv"
)

FULL_ORDER_FILE = (
    ORDER_DIR
    / "full_9_species_plot_order.tsv"
)

VERONICA_ORDER_FILE = (
    ORDER_DIR
    / "veronica_only_plot_order.tsv"
)

CHROMOSOME_MAPPING_FILE = (
    TABLE_DIR
    / "refined_numeric_chromosome_mapping.tsv"
)

required_files = [
    INPUT_MANIFEST,
    STEP34_RESULTS,
    STEP34G_RESULTS,
    FULL_ORDER_FILE,
    VERONICA_ORDER_FILE,
    CHROMOSOME_MAPPING_FILE,
]

for path in required_files:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(
            f"ERROR: Required input missing or empty: {path}"
        )

PLOT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


# ============================================================
# Plot configuration
# ============================================================

SPECIES_COLOURS = {
    "PMAJ": "#4C78A8",
    "VANA": "#8E5BD9",
    "VSER": "#2AA198",
    "VPER": "#E69F00",
    "VVER": "#CC79A7",
    "VARV": "#0072B2",
    "VPAN": "#D55E00",
    "VTRI": "#E64B35",
    "VSCU": "#59A14F",
}

PLOIDY_LABELS = {
    "PMAJ": "2x",
    "VANA": "4x",
    "VSER": "2x",
    "VPER": "4x",
    "VVER": "2x",
    "VARV": "2x",
    "VPAN": "2x",
    "VTRI": "2x",
    "VSCU": "2x",
}

NORMAL_LINK_COLOUR = "#4A90C2"
REARRANGEMENT_LINK_COLOUR = "#D73027"

NORMAL_LINK_ALPHA = 0.055
REARRANGEMENT_LINK_ALPHA = 0.115

NORMAL_LINK_WIDTH = 0.26
REARRANGEMENT_LINK_WIDTH = 0.42

MIN_STRONG_ANCHORS = 20
RELATIVE_STRONG_THRESHOLD = 0.25

MAXIMUM_LINKS_PER_EDGE = 4200

CHROMOSOME_GAP = 0.018
CHROMOSOME_HEIGHT = 0.105

LEFT_LABEL_X = -0.105
RATIO_LABEL_X = 1.018


# ============================================================
# General table functions
# ============================================================

def read_tsv(path):
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
        )

        writer.writeheader()
        writer.writerows(rows)


def read_plot_order(path):
    rows = read_tsv(path)

    rows.sort(
        key=lambda row: int(
            row["plot_rank"]
        )
    )

    return [
        row["species_code"]
        for row in rows
    ]


# ============================================================
# Read species metadata
# ============================================================

manifest_rows = read_tsv(
    INPUT_MANIFEST
)

species_metadata = {
    row["species_code"]: row
    for row in manifest_rows
}

full_order = read_plot_order(
    FULL_ORDER_FILE
)

veronica_order = read_plot_order(
    VERONICA_ORDER_FILE
)

for code in set(
    full_order
    + veronica_order
):
    if code not in species_metadata:
        raise SystemExit(
            f"ERROR: {code} absent from verified manifest."
        )


# ============================================================
# Read chromosome-number mapping generated by Step 34H
# ============================================================

mapping_rows = read_tsv(
    CHROMOSOME_MAPPING_FILE
)

original_to_numeric = defaultdict(
    dict
)

chromosome_lengths = defaultdict(
    dict
)

ordered_chromosomes = defaultdict(
    list
)

for row in mapping_rows:
    code = row["species_code"]
    original_id = row[
        "original_sequence_id"
    ]
    numeric_id = row[
        "numeric_chromosome"
    ]
    length_bp = int(
        row["length_bp"]
    )

    original_to_numeric[
        code
    ][original_id] = numeric_id

    chromosome_lengths[
        code
    ][numeric_id] = length_bp

for code in original_to_numeric:
    ordered_chromosomes[
        code
    ] = sorted(
        chromosome_lengths[
            code
        ],
        key=lambda value: int(value),
    )


# ============================================================
# Read JCVI BED files
# ============================================================

def read_numeric_gene_coordinates(
    species_code,
):
    bed_file = (
        JCVI_BED_DIR
        / f"{species_code}.bed"
    )

    if not bed_file.is_file():
        raise SystemExit(
            f"ERROR: BED missing for {species_code}: {bed_file}"
        )

    genes = {}

    with bed_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) < 4:
                continue

            original_chromosome = fields[0]

            if (
                original_chromosome
                not in original_to_numeric[
                    species_code
                ]
            ):
                continue

            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue

            gene_id = fields[3]

            genes[gene_id] = {
                "chromosome": (
                    original_to_numeric[
                        species_code
                    ][original_chromosome]
                ),
                "midpoint": (
                    start + end
                ) / 2.0,
            }

    if not genes:
        raise SystemExit(
            f"ERROR: No genes loaded for {species_code}."
        )

    return genes


species_genes = {
    code: read_numeric_gene_coordinates(
        code
    )
    for code in full_order
}


# ============================================================
# Read JCVI comparison files
# ============================================================

step34_rows = read_tsv(
    STEP34_RESULTS
)

step34g_rows = read_tsv(
    STEP34G_RESULTS
)

comparison_lookup = {}

for row in step34_rows + step34g_rows:
    if row.get("status", "") != "PASS":
        continue

    query = row["query_species"]
    subject = row["subject_species"]

    anchor_file = row.get(
        "anchor_file",
        ""
    )

    if not anchor_file:
        continue

    anchor_path = Path(
        anchor_file
    )

    if (
        not anchor_path.is_file()
        or anchor_path.stat().st_size == 0
    ):
        continue

    key = tuple(
        sorted(
            [
                query,
                subject,
            ]
        )
    )

    comparison_lookup[key] = {
        "comparison_id": row[
            "comparison_id"
        ],
        "stored_query": query,
        "stored_subject": subject,
        "anchor_file": anchor_path,
    }


def build_edges(
    plot_name,
    order,
):
    edges = []

    for index in range(
        len(order) - 1
    ):
        upper = order[index]
        lower = order[index + 1]

        key = tuple(
            sorted(
                [
                    upper,
                    lower,
                ]
            )
        )

        if key not in comparison_lookup:
            raise SystemExit(
                "ERROR: Missing anchor comparison for "
                f"{upper}-{lower}"
            )

        comparison = comparison_lookup[
            key
        ]

        edges.append(
            {
                "plot_name": plot_name,
                "edge_rank": index + 1,
                "upper_species": upper,
                "lower_species": lower,
                **comparison,
            }
        )

    return edges


full_edges = build_edges(
    "FULL_9_SPECIES",
    full_order,
)

veronica_edges = build_edges(
    "VERONICA_ONLY",
    veronica_order,
)


# ============================================================
# Read full JCVI anchor pairs
# ============================================================

def read_anchor_pairs(path):
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


def strong_partners(
    partner_counts,
):
    if not partner_counts:
        return set()

    strongest = max(
        partner_counts.values()
    )

    return {
        partner
        for partner, count
        in partner_counts.items()
        if (
            count
            >= MIN_STRONG_ANCHORS
            and count
            >= RELATIVE_STRONG_THRESHOLD
            * strongest
        )
    }


def analyse_edge(edge):
    upper = edge[
        "upper_species"
    ]
    lower = edge[
        "lower_species"
    ]

    stored_query = edge[
        "stored_query"
    ]

    full_pairs = read_anchor_pairs(
        edge["anchor_file"]
    )

    upper_genes = species_genes[
        upper
    ]

    lower_genes = species_genes[
        lower
    ]

    mapped_pairs = []
    chromosome_pair_counts = Counter()

    for query_gene, subject_gene in full_pairs:
        if stored_query == upper:
            upper_gene = query_gene
            lower_gene = subject_gene
        else:
            upper_gene = subject_gene
            lower_gene = query_gene

        if (
            upper_gene not in upper_genes
            or lower_gene not in lower_genes
        ):
            continue

        upper_chromosome = (
            upper_genes[
                upper_gene
            ]["chromosome"]
        )

        lower_chromosome = (
            lower_genes[
                lower_gene
            ]["chromosome"]
        )

        mapped_pairs.append(
            (
                upper_gene,
                lower_gene,
                upper_chromosome,
                lower_chromosome,
            )
        )

        chromosome_pair_counts[
            (
                upper_chromosome,
                lower_chromosome,
            )
        ] += 1

    if not mapped_pairs:
        raise SystemExit(
            "ERROR: No anchor pairs mapped for "
            f"{upper}-{lower}."
        )

    upper_to_lower = defaultdict(
        Counter
    )

    lower_to_upper = defaultdict(
        Counter
    )

    for (
        upper_chromosome,
        lower_chromosome,
    ), count in chromosome_pair_counts.items():
        upper_to_lower[
            upper_chromosome
        ][lower_chromosome] = count

        lower_to_upper[
            lower_chromosome
        ][upper_chromosome] = count

    upper_strong = {
        chromosome: strong_partners(
            counts
        )
        for chromosome, counts
        in upper_to_lower.items()
    }

    lower_strong = {
        chromosome: strong_partners(
            counts
        )
        for chromosome, counts
        in lower_to_upper.items()
    }

    rearrangement_pairs = set()

    for chromosome, partners in (
        upper_strong.items()
    ):
        if len(partners) > 1:
            for partner in partners:
                rearrangement_pairs.add(
                    (
                        chromosome,
                        partner,
                    )
                )

    for chromosome, partners in (
        lower_strong.items()
    ):
        if len(partners) > 1:
            for partner in partners:
                rearrangement_pairs.add(
                    (
                        partner,
                        chromosome,
                    )
                )

    graph = defaultdict(set)

    for upper_chromosome, partners in (
        upper_strong.items()
    ):
        for lower_chromosome in partners:
            upper_node = (
                "U",
                upper_chromosome,
            )

            lower_node = (
                "L",
                lower_chromosome,
            )

            graph[
                upper_node
            ].add(
                lower_node
            )

            graph[
                lower_node
            ].add(
                upper_node
            )

    component_ratios = Counter()
    visited = set()

    for starting_node in graph:
        if starting_node in visited:
            continue

        queue = deque(
            [starting_node]
        )

        visited.add(
            starting_node
        )

        upper_component = set()
        lower_component = set()

        while queue:
            side, chromosome = (
                queue.popleft()
            )

            if side == "U":
                upper_component.add(
                    chromosome
                )
            else:
                lower_component.add(
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

        if (
            upper_component
            and lower_component
        ):
            ratio = (
                f"{len(upper_component)}:"
                f"{len(lower_component)}"
            )

            weight = sum(
                chromosome_pair_counts.get(
                    (
                        upper_chromosome,
                        lower_chromosome,
                    ),
                    0,
                )
                for upper_chromosome
                in upper_component
                for lower_chromosome
                in lower_component
            )

            component_ratios[
                ratio
            ] += weight

    modal_ratio = (
        component_ratios.most_common(
            1
        )[0][0]
        if component_ratios
        else "NA"
    )

    return {
        "mapped_pairs": mapped_pairs,
        "chromosome_pair_counts": (
            chromosome_pair_counts
        ),
        "rearrangement_pairs": (
            rearrangement_pairs
        ),
        "modal_ratio": modal_ratio,
    }


# ============================================================
# Chromosome track layout
# ============================================================

def build_track_layout(
    species_order,
):
    layouts = {}

    for code in species_order:
        chromosomes = (
            ordered_chromosomes[
                code
            ]
        )

        lengths = chromosome_lengths[
            code
        ]

        total_length = sum(
            lengths[
                chromosome
            ]
            for chromosome
            in chromosomes
        )

        usable_width = (
            1.0
            - CHROMOSOME_GAP
            * max(
                0,
                len(chromosomes) - 1,
            )
        )

        cursor = 0.0
        track = {}

        for chromosome in chromosomes:
            width = (
                usable_width
                * lengths[
                    chromosome
                ]
                / total_length
            )

            track[
                chromosome
            ] = {
                "start": cursor,
                "end": cursor + width,
                "midpoint": (
                    cursor
                    + width / 2.0
                ),
                "length": lengths[
                    chromosome
                ],
            }

            cursor += (
                width
                + CHROMOSOME_GAP
            )

        layouts[
            code
        ] = track

    return layouts


def gene_x(
    code,
    gene_id,
    track_layout,
):
    gene = species_genes[
        code
    ][gene_id]

    chromosome = gene[
        "chromosome"
    ]

    layout = track_layout[
        code
    ][chromosome]

    relative_position = (
        gene["midpoint"]
        / layout["length"]
    )

    relative_position = min(
        1.0,
        max(
            0.0,
            relative_position,
        ),
    )

    return (
        layout["start"]
        + relative_position
        * (
            layout["end"]
            - layout["start"]
        )
    )


def deterministic_subsample(
    records,
    maximum_records,
):
    if len(records) <= maximum_records:
        return records

    step = (
        len(records)
        / maximum_records
    )

    sampled = []
    position = 0.0

    while (
        len(sampled)
        < maximum_records
        and int(position)
        < len(records)
    ):
        sampled.append(
            records[
                int(position)
            ]
        )

        position += step

    return sampled


# ============================================================
# Main drawing function
# ============================================================

def draw_plot(
    plot_name,
    species_order,
    edges,
    output_prefix,
):
    track_layout = build_track_layout(
        species_order
    )

    y_positions = {
        code: (
            len(species_order)
            - 1
            - index
        )
        for index, code
        in enumerate(species_order)
    }

    analyses = {}

    ratio_rows = []
    rearrangement_rows = []

    for edge in edges:
        analysis = analyse_edge(
            edge
        )

        analyses[
            edge["edge_rank"]
        ] = analysis

        ratio_rows.append(
            {
                "plot_name": plot_name,
                "edge_rank": edge[
                    "edge_rank"
                ],
                "upper_species": edge[
                    "upper_species"
                ],
                "lower_species": edge[
                    "lower_species"
                ],
                "comparison_id": edge[
                    "comparison_id"
                ],
                "mapped_anchor_pairs": len(
                    analysis[
                        "mapped_pairs"
                    ]
                ),
                "operational_component_ratio": (
                    analysis[
                        "modal_ratio"
                    ]
                ),
                "multiple_strong_partner_pairs": len(
                    analysis[
                        "rearrangement_pairs"
                    ]
                ),
                "minimum_anchor_threshold": (
                    MIN_STRONG_ANCHORS
                ),
                "relative_partner_threshold": (
                    RELATIVE_STRONG_THRESHOLD
                ),
            }
        )

        for (
            upper_chromosome,
            lower_chromosome,
        ) in sorted(
            analysis[
                "rearrangement_pairs"
            ]
        ):
            rearrangement_rows.append(
                {
                    "plot_name": plot_name,
                    "edge_rank": edge[
                        "edge_rank"
                    ],
                    "upper_species": edge[
                        "upper_species"
                    ],
                    "lower_species": edge[
                        "lower_species"
                    ],
                    "upper_chromosome": (
                        upper_chromosome
                    ),
                    "lower_chromosome": (
                        lower_chromosome
                    ),
                    "anchor_pairs": (
                        analysis[
                            "chromosome_pair_counts"
                        ][
                            (
                                upper_chromosome,
                                lower_chromosome,
                            )
                        ]
                    ),
                    "classification": (
                        "MULTIPLE_STRONG_PARTNER_RELATIONSHIP"
                    ),
                }
            )

    figure_height = max(
        10.5,
        1.45 * len(
            species_order
        ),
    )

    figure, axis = plt.subplots(
        figsize=(
            19,
            figure_height,
        )
    )

    axis.set_xlim(
        -0.13,
        1.105,
    )

    axis.set_ylim(
        -0.58,
        len(species_order) - 0.42,
    )

    axis.axis(
        "off"
    )

    # --------------------------------------------------------
    # Draw synteny links first
    # --------------------------------------------------------

    for edge in edges:
        upper = edge[
            "upper_species"
        ]

        lower = edge[
            "lower_species"
        ]

        upper_y = y_positions[
            upper
        ]

        lower_y = y_positions[
            lower
        ]

        analysis = analyses[
            edge["edge_rank"]
        ]

        mapped_pairs = (
            deterministic_subsample(
                analysis[
                    "mapped_pairs"
                ],
                MAXIMUM_LINKS_PER_EDGE,
            )
        )

        for (
            upper_gene,
            lower_gene,
            upper_chromosome,
            lower_chromosome,
        ) in mapped_pairs:
            upper_x = gene_x(
                upper,
                upper_gene,
                track_layout,
            )

            lower_x = gene_x(
                lower,
                lower_gene,
                track_layout,
            )

            middle_y = (
                upper_y
                + lower_y
            ) / 2.0

            vertices = [
                (
                    upper_x,
                    upper_y
                    - CHROMOSOME_HEIGHT
                    / 2.0,
                ),
                (
                    upper_x,
                    middle_y,
                ),
                (
                    lower_x,
                    middle_y,
                ),
                (
                    lower_x,
                    lower_y
                    + CHROMOSOME_HEIGHT
                    / 2.0,
                ),
            ]

            codes = [
                MplPath.MOVETO,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CURVE4,
            ]

            is_rearrangement = (
                (
                    upper_chromosome,
                    lower_chromosome,
                )
                in analysis[
                    "rearrangement_pairs"
                ]
            )

            if is_rearrangement:
                link_colour = (
                    REARRANGEMENT_LINK_COLOUR
                )

                link_alpha = (
                    REARRANGEMENT_LINK_ALPHA
                )

                link_width = (
                    REARRANGEMENT_LINK_WIDTH
                )
            else:
                link_colour = (
                    NORMAL_LINK_COLOUR
                )

                link_alpha = (
                    NORMAL_LINK_ALPHA
                )

                link_width = (
                    NORMAL_LINK_WIDTH
                )

            axis.add_patch(
                PathPatch(
                    MplPath(
                        vertices,
                        codes,
                    ),
                    facecolor="none",
                    edgecolor=link_colour,
                    linewidth=link_width,
                    alpha=link_alpha,
                    zorder=1,
                )
            )

        ratio = analysis[
            "modal_ratio"
        ]

        axis.text(
            RATIO_LABEL_X,
            (
                upper_y
                + lower_y
            ) / 2.0,
            (
                f"{upper}-{lower}\n"
                f"ratio {ratio}"
            ),
            fontsize=8.5,
            ha="left",
            va="center",
            color="#333333",
        )

    # --------------------------------------------------------
    # Draw species labels and chromosome bars
    # --------------------------------------------------------

    for code in species_order:
        y = y_positions[
            code
        ]

        colour = SPECIES_COLOURS[
            code
        ]

        ploidy = PLOIDY_LABELS[
            code
        ]

        scientific_name = (
            species_metadata[
                code
            ]["scientific_name"]
        )

        axis.text(
            LEFT_LABEL_X,
            y + 0.045,
            f"{ploidy}  {code}",
            fontsize=14,
            fontweight="bold",
            ha="right",
            va="center",
            color=colour,
        )

        axis.text(
            LEFT_LABEL_X,
            y - 0.145,
            scientific_name,
            fontsize=8.1,
            fontstyle="italic",
            ha="right",
            va="center",
            color="#555555",
        )

        for chromosome in (
            ordered_chromosomes[
                code
            ]
        ):
            layout = track_layout[
                code
            ][chromosome]

            chromosome_patch = (
                FancyBboxPatch(
                    (
                        layout["start"],
                        y
                        - CHROMOSOME_HEIGHT
                        / 2.0,
                    ),
                    layout["end"]
                    - layout["start"],
                    CHROMOSOME_HEIGHT,
                    boxstyle=(
                        "round,pad=0.002,"
                        "rounding_size=0.010"
                    ),
                    facecolor=colour,
                    edgecolor="#222222",
                    linewidth=0.75,
                    alpha=0.96,
                    zorder=3,
                )
            )

            axis.add_patch(
                chromosome_patch
            )

            axis.text(
                layout["midpoint"],
                y + 0.135,
                chromosome,
                fontsize=10.5,
                fontweight="bold",
                ha="center",
                va="bottom",
                color="#222222",
                zorder=4,
            )

    if plot_name == "FULL_9_SPECIES":
        title = (
            "Macrosynteny overview of nine species"
        )
    else:
        title = (
            "Macrosynteny overview of eight "
            "Veronica species"
        )

    axis.set_title(
        title,
        fontsize=17,
        fontweight="bold",
        pad=20,
    )

    figure.text(
        0.012,
        0.014,
        (
            "Blue links: mapped full-anchor synteny pairs. "
            "Red links: chromosome pairs belonging to "
            "multiple strong-partner relationships. "
            "Ratios are operational chromosome-component "
            "ratios."
        ),
        fontsize=8.5,
        ha="left",
        va="bottom",
        color="#333333",
    )

    figure.tight_layout(
        rect=(
            0.035,
            0.045,
            0.94,
            0.97,
        )
    )

    pdf_file = Path(
        str(output_prefix)
        + ".pdf"
    )

    svg_file = Path(
        str(output_prefix)
        + ".svg"
    )

    png_file = Path(
        str(output_prefix)
        + ".png"
    )

    figure.savefig(
        pdf_file,
        bbox_inches="tight",
    )

    figure.savefig(
        svg_file,
        bbox_inches="tight",
    )

    figure.savefig(
        png_file,
        dpi=400,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

    return {
        "ratio_rows": ratio_rows,
        "rearrangement_rows": (
            rearrangement_rows
        ),
        "pdf": pdf_file,
        "svg": svg_file,
        "png": png_file,
    }


# ============================================================
# Generate the two figures
# ============================================================

full_outputs = draw_plot(
    "FULL_9_SPECIES",
    full_order,
    full_edges,
    PLOT_DIR
    / "full_9_species.macrosynteny.publication_style",
)

veronica_outputs = draw_plot(
    "VERONICA_ONLY",
    veronica_order,
    veronica_edges,
    PLOT_DIR
    / "veronica_only_8_species.macrosynteny.publication_style",
)


# ============================================================
# Write supporting tables
# ============================================================

ratio_fields = [
    "plot_name",
    "edge_rank",
    "upper_species",
    "lower_species",
    "comparison_id",
    "mapped_anchor_pairs",
    "operational_component_ratio",
    "multiple_strong_partner_pairs",
    "minimum_anchor_threshold",
    "relative_partner_threshold",
]

rearrangement_fields = [
    "plot_name",
    "edge_rank",
    "upper_species",
    "lower_species",
    "upper_chromosome",
    "lower_chromosome",
    "anchor_pairs",
    "classification",
]

write_tsv(
    TABLE_DIR
    / "publication_full_plot_synteny_ratios.tsv",
    full_outputs[
        "ratio_rows"
    ],
    ratio_fields,
)

write_tsv(
    TABLE_DIR
    / "publication_veronica_plot_synteny_ratios.tsv",
    veronica_outputs[
        "ratio_rows"
    ],
    ratio_fields,
)

write_tsv(
    TABLE_DIR
    / "publication_full_plot_rearrangement_pairs.tsv",
    full_outputs[
        "rearrangement_rows"
    ],
    rearrangement_fields,
)

write_tsv(
    TABLE_DIR
    / "publication_veronica_plot_rearrangement_pairs.tsv",
    veronica_outputs[
        "rearrangement_rows"
    ],
    rearrangement_fields,
)

summary_rows = [
    {
        "metric": "full_plot_species",
        "value": len(full_order),
    },
    {
        "metric": "veronica_plot_species",
        "value": len(veronica_order),
    },
    {
        "metric": "full_plot_edges",
        "value": len(full_edges),
    },
    {
        "metric": "veronica_plot_edges",
        "value": len(veronica_edges),
    },
    {
        "metric": "full_plot_rearrangement_pairs",
        "value": len(
            full_outputs[
                "rearrangement_rows"
            ]
        ),
    },
    {
        "metric": "veronica_plot_rearrangement_pairs",
        "value": len(
            veronica_outputs[
                "rearrangement_rows"
            ]
        ),
    },
    {
        "metric": "full_plot_pdf",
        "value": str(
            full_outputs["pdf"]
        ),
    },
    {
        "metric": "full_plot_svg",
        "value": str(
            full_outputs["svg"]
        ),
    },
    {
        "metric": "full_plot_png",
        "value": str(
            full_outputs["png"]
        ),
    },
    {
        "metric": "veronica_plot_pdf",
        "value": str(
            veronica_outputs["pdf"]
        ),
    },
    {
        "metric": "veronica_plot_svg",
        "value": str(
            veronica_outputs["svg"]
        ),
    },
    {
        "metric": "veronica_plot_png",
        "value": str(
            veronica_outputs["png"]
        ),
    },
    {
        "metric": "status",
        "value": "PASS",
    },
]

write_tsv(
    TABLE_DIR
    / "publication_macrosynteny_plot_summary.tsv",
    summary_rows,
    [
        "metric",
        "value",
    ],
)

print(
    "Full plot order:",
    ",".join(
        full_order
    ),
)

print(
    "Veronica-only order:",
    ",".join(
        veronica_order
    ),
)

print(
    "Full plot rearrangement pairs:",
    len(
        full_outputs[
            "rearrangement_rows"
        ]
    ),
)

print(
    "Veronica-only rearrangement pairs:",
    len(
        veronica_outputs[
            "rearrangement_rows"
        ]
    ),
)

print(
    "Step 34I plotting: PASS"
)
PY

for FILE in \
    "${PLOT_DIR}/full_9_species.macrosynteny.publication_style.pdf" \
    "${PLOT_DIR}/full_9_species.macrosynteny.publication_style.svg" \
    "${PLOT_DIR}/full_9_species.macrosynteny.publication_style.png" \
    "${PLOT_DIR}/veronica_only_8_species.macrosynteny.publication_style.pdf" \
    "${PLOT_DIR}/veronica_only_8_species.macrosynteny.publication_style.svg" \
    "${PLOT_DIR}/veronica_only_8_species.macrosynteny.publication_style.png" \
    "${TABLE_DIR}/publication_full_plot_synteny_ratios.tsv" \
    "${TABLE_DIR}/publication_veronica_plot_synteny_ratios.tsv" \
    "${TABLE_DIR}/publication_full_plot_rearrangement_pairs.tsv" \
    "${TABLE_DIR}/publication_veronica_plot_rearrangement_pairs.tsv" \
    "${TABLE_DIR}/publication_macrosynteny_plot_summary.tsv"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Expected output is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

echo
echo "============================================================"
echo "Publication-style plot summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/publication_macrosynteny_plot_summary.tsv"

echo
echo "============================================================"
echo "Full plot ratios"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/publication_full_plot_synteny_ratios.tsv"

cat > "${CHECKPOINT_DIR}/PUBLICATION_MACROSYNTENY_PLOTS_COMPLETE.txt" <<EOF2
checkpoint=publication_style_macrosynteny_plots
date=$(date --iso-8601=seconds)
status=PASS
chromosome_colouring=species_specific
species_labels=ploidy_plus_abbreviation
normal_links=blue_full_anchor_pairs
rearrangement_links=red_multiple_strong_partner_pairs
ratios=operational_chromosome_component_ratios
next_step=inspect_publication_style_figures
EOF2

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name "sha256_checksums.txt" \
    -print0 |
sort -z |
xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "Step 34I completed successfully."
echo "Checkpoint:"
cat \
    "${CHECKPOINT_DIR}/PUBLICATION_MACROSYNTENY_PLOTS_COMPLETE.txt"
