#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=macro34h
#SBATCH --output=10_synteny/logs/macro34h_%j.out
#SBATCH --error=10_synteny/logs/macro34h_%j.err

set -euo pipefail

export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

REFINED_DIR="${SYNTENY_DIR}/refined_macrosynteny"
PLOT_DIR="${REFINED_DIR}/plots"
TABLE_DIR="${REFINED_DIR}/tables"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_refined_macrosynteny_plots"

mkdir -p \
    "${SYNTENY_DIR}/logs" \
    "${PLOT_DIR}" \
    "${TABLE_DIR}" \
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
    echo "ERROR: Conda could not be found." >&2
    exit 1
fi

eval "$("${CONDA_EXE_PATH}" shell.bash hook)"
conda activate jcvi_env

export MPLBACKEND=Agg

echo "Python:"
python --version

rm -f \
    "${CHECKPOINT_DIR}"/* \
    "${TABLE_DIR}/refined_full_plot_edge_mapping_qc.tsv" \
    "${TABLE_DIR}/refined_veronica_only_edge_mapping_qc.tsv" \
    "${TABLE_DIR}/refined_full_plot_synteny_ratios.tsv" \
    "${TABLE_DIR}/refined_veronica_only_plot_synteny_ratios.tsv" \
    "${TABLE_DIR}/refined_full_plot_rearrangement_candidates.tsv" \
    "${TABLE_DIR}/refined_veronica_only_plot_rearrangement_candidates.tsv" \
    "${TABLE_DIR}/refined_macrosynteny_plot_summary.tsv"

python - <<'PY'
from __future__ import annotations

import csv
import gzip
import math
import re
from collections import Counter, defaultdict, deque
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch


PROJECT_DIR = Path(
    os.environ["PROJECT_DIR"]
).resolve()

SYNTENY_DIR = PROJECT_DIR / "10_synteny"

REFINED_DIR = (
    SYNTENY_DIR
    / "refined_macrosynteny"
)

PLOT_DIR = (
    REFINED_DIR
    / "plots"
)

TABLE_DIR = (
    REFINED_DIR
    / "tables"
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

for required_file in [
    INPUT_MANIFEST,
    STEP34_RESULTS,
    STEP34G_RESULTS,
    FULL_ORDER_FILE,
    VERONICA_ORDER_FILE,
]:
    if (
        not required_file.is_file()
        or required_file.stat().st_size == 0
    ):
        raise SystemExit(
            "ERROR: Required input is missing or empty: "
            f"{required_file}"
        )

PLOT_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

TABLE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)


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


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(
            path,
            "rt",
            encoding="utf-8",
            errors="replace",
        )

    return path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    )


def read_fasta_lengths(path):
    lengths = {}
    descriptions = {}

    current_id = None
    current_description = None
    current_length = 0

    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current_id is not None:
                    lengths[
                        current_id
                    ] = current_length

                    descriptions[
                        current_id
                    ] = current_description

                current_description = (
                    line[1:].strip()
                )

                current_id = (
                    current_description.split()[0]
                )

                current_length = 0
            else:
                current_length += len(
                    line.strip()
                )

    if current_id is not None:
        lengths[
            current_id
        ] = current_length

        descriptions[
            current_id
        ] = current_description

    return lengths, descriptions


def chromosome_number(
    seqid,
    description,
):
    search_strings = [
        description or "",
        seqid,
    ]

    patterns = [
        r"\bchromosome[_ :.-]*(\d+)\b",
        r"\bchr[_ :.-]*(\d+)\b",
        r"(?:^|[_-])(\d+)$",
        r"(\d+)$",
    ]

    for source in search_strings:
        for pattern in patterns:
            match = re.search(
                pattern,
                source,
                flags=re.IGNORECASE,
            )

            if match:
                return int(
                    match.group(1)
                )

    return None


def read_bed(path):
    genes = {}
    chromosome_maximum = defaultdict(int)

    with path.open(
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

            chromosome = fields[0]

            try:
                start = int(
                    fields[1]
                )

                end = int(
                    fields[2]
                )
            except ValueError:
                continue

            gene_id = fields[3]

            genes[gene_id] = {
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "midpoint": (
                    start + end
                ) / 2.0,
            }

            chromosome_maximum[
                chromosome
            ] = max(
                chromosome_maximum[
                    chromosome
                ],
                end,
            )

    return genes, chromosome_maximum


def read_full_anchor_pairs(path):
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


def species_track_colour(
    species_code,
    ploidy,
):
    if species_code == "PMAJ":
        return "#4d4d4d"

    if str(ploidy).startswith("4"):
        return "#7a3db8"

    return "#2166ac"


def build_numeric_chromosome_order(
    fasta_lengths,
    fasta_descriptions,
    bed_chromosomes,
):
    records = []

    for sequence_id in bed_chromosomes:
        numeric_value = chromosome_number(
            sequence_id,
            fasta_descriptions.get(
                sequence_id,
                "",
            ),
        )

        records.append(
            {
                "sequence_id": sequence_id,
                "numeric_value": (
                    numeric_value
                    if numeric_value is not None
                    else math.inf
                ),
                "numeric_label": (
                    str(numeric_value)
                    if numeric_value is not None
                    else ""
                ),
                "length": fasta_lengths.get(
                    sequence_id,
                    bed_chromosomes[
                        sequence_id
                    ],
                ),
            }
        )

    records.sort(
        key=lambda row: (
            row["numeric_value"],
            row["sequence_id"],
        )
    )

    fallback_number = 1

    used_labels = {
        row["numeric_label"]
        for row in records
        if row["numeric_label"]
    }

    for row in records:
        if row["numeric_label"]:
            continue

        while str(
            fallback_number
        ) in used_labels:
            fallback_number += 1

        row["numeric_label"] = str(
            fallback_number
        )

        used_labels.add(
            row["numeric_label"]
        )

        fallback_number += 1

    return records


species_rows = read_tsv(
    INPUT_MANIFEST
)

species_metadata = {
    row["species_code"]: row
    for row in species_rows
}

full_order = read_plot_order(
    FULL_ORDER_FILE
)

veronica_order = read_plot_order(
    VERONICA_ORDER_FILE
)

required_species = set(
    full_order
    + veronica_order
)

for species_code in required_species:
    if species_code not in species_metadata:
        raise SystemExit(
            "ERROR: Species absent from verified manifest: "
            f"{species_code}"
        )

species_assets = {}
chromosome_mapping_rows = []

for species_code in full_order:
    metadata = species_metadata[
        species_code
    ]

    genome_path = Path(
        metadata["genome_fasta"]
    )

    bed_path = (
        JCVI_BED_DIR
        / f"{species_code}.bed"
    )

    if not genome_path.is_file():
        raise SystemExit(
            "ERROR: Genome FASTA missing for "
            f"{species_code}: {genome_path}"
        )

    if not bed_path.is_file():
        raise SystemExit(
            "ERROR: BED file missing for "
            f"{species_code}: {bed_path}"
        )

    fasta_lengths, fasta_descriptions = (
        read_fasta_lengths(
            genome_path
        )
    )

    genes, bed_chromosome_maximum = (
        read_bed(
            bed_path
        )
    )

    chromosome_records = (
        build_numeric_chromosome_order(
            fasta_lengths,
            fasta_descriptions,
            bed_chromosome_maximum,
        )
    )

    original_to_numeric = {
        row["sequence_id"]: row[
            "numeric_label"
        ]
        for row in chromosome_records
    }

    chromosome_lengths = {
        row["numeric_label"]: int(
            row["length"]
        )
        for row in chromosome_records
    }

    ordered_chromosomes = [
        row["numeric_label"]
        for row in chromosome_records
    ]

    numeric_genes = {}

    for gene_id, gene in genes.items():
        original_chromosome = gene[
            "chromosome"
        ]

        if (
            original_chromosome
            not in original_to_numeric
        ):
            continue

        numeric_chromosome = (
            original_to_numeric[
                original_chromosome
            ]
        )

        numeric_genes[
            gene_id
        ] = {
            "chromosome": (
                numeric_chromosome
            ),
            "midpoint": gene[
                "midpoint"
            ],
        }

    for row in chromosome_records:
        chromosome_mapping_rows.append(
            {
                "species_code": species_code,
                "original_sequence_id": row[
                    "sequence_id"
                ],
                "numeric_chromosome": row[
                    "numeric_label"
                ],
                "length_bp": row[
                    "length"
                ],
            }
        )

    species_assets[
        species_code
    ] = {
        "genes": numeric_genes,
        "chromosome_lengths": (
            chromosome_lengths
        ),
        "ordered_chromosomes": (
            ordered_chromosomes
        ),
        "ploidy": metadata[
            "ploidy"
        ],
        "scientific_name": metadata[
            "scientific_name"
        ],
    }

write_tsv(
    TABLE_DIR
    / "refined_numeric_chromosome_mapping.tsv",
    chromosome_mapping_rows,
    [
        "species_code",
        "original_sequence_id",
        "numeric_chromosome",
        "length_bp",
    ],
)

step34_rows = read_tsv(
    STEP34_RESULTS
)

step34g_rows = read_tsv(
    STEP34G_RESULTS
)

comparison_lookup = {}

for row in (
    step34_rows
    + step34g_rows
):
    if row.get(
        "status",
        ""
    ) != "PASS":
        continue

    query_species = row[
        "query_species"
    ]

    subject_species = row[
        "subject_species"
    ]

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
                query_species,
                subject_species,
            ]
        )
    )

    comparison_lookup[key] = {
        "comparison_id": row[
            "comparison_id"
        ],
        "stored_query": query_species,
        "stored_subject": (
            subject_species
        ),
        "anchor_file": str(
            anchor_path
        ),
    }


def build_edge_rows(
    plot_name,
    species_order,
):
    edge_rows = []

    for edge_index in range(
        len(species_order) - 1
    ):
        upper_species = (
            species_order[
                edge_index
            ]
        )

        lower_species = (
            species_order[
                edge_index + 1
            ]
        )

        lookup_key = tuple(
            sorted(
                [
                    upper_species,
                    lower_species,
                ]
            )
        )

        if lookup_key not in comparison_lookup:
            raise SystemExit(
                "ERROR: No full anchor file for "
                f"{upper_species}-{lower_species}"
            )

        comparison = comparison_lookup[
            lookup_key
        ]

        edge_rows.append(
            {
                "plot_name": plot_name,
                "edge_rank": (
                    edge_index + 1
                ),
                "upper_species": (
                    upper_species
                ),
                "lower_species": (
                    lower_species
                ),
                **comparison,
            }
        )

    return edge_rows


full_edges = build_edge_rows(
    "FULL_9_SPECIES",
    full_order,
)

veronica_edges = build_edge_rows(
    "VERONICA_ONLY",
    veronica_order,
)


def build_track_layout(
    species_order,
):
    layout = {}

    for species_code in species_order:
        chromosomes = (
            species_assets[
                species_code
            ]["ordered_chromosomes"]
        )

        lengths = species_assets[
            species_code
        ]["chromosome_lengths"]

        total_length = sum(
            lengths[
                chromosome
            ]
            for chromosome
            in chromosomes
        )

        gap = 0.012

        usable_width = (
            1.0
            - gap
            * max(
                0,
                len(chromosomes) - 1,
            )
        )

        cursor = 0.0

        species_layout = {}

        for chromosome in chromosomes:
            width = (
                usable_width
                * lengths[
                    chromosome
                ]
                / total_length
            )

            species_layout[
                chromosome
            ] = {
                "start": cursor,
                "end": (
                    cursor + width
                ),
                "midpoint": (
                    cursor
                    + width / 2.0
                ),
                "length": lengths[
                    chromosome
                ],
            }

            cursor += (
                width + gap
            )

        layout[
            species_code
        ] = species_layout

    return layout


def gene_x_position(
    species_code,
    gene_id,
    track_layout,
):
    gene = species_assets[
        species_code
    ]["genes"].get(
        gene_id
    )

    if gene is None:
        return None

    chromosome = gene[
        "chromosome"
    ]

    chromosome_layout = (
        track_layout[
            species_code
        ].get(
            chromosome
        )
    )

    if chromosome_layout is None:
        return None

    chromosome_length = (
        chromosome_layout[
            "length"
        ]
    )

    if chromosome_length <= 0:
        return None

    relative_position = (
        gene["midpoint"]
        / chromosome_length
    )

    relative_position = max(
        0.0,
        min(
            1.0,
            relative_position,
        ),
    )

    return (
        chromosome_layout["start"]
        + relative_position
        * (
            chromosome_layout["end"]
            - chromosome_layout["start"]
        )
    )


MIN_STRONG_ANCHORS = 20
RELATIVE_STRONG_THRESHOLD = 0.25


def strong_partner_set(
    partner_counts,
):
    if not partner_counts:
        return set()

    maximum_count = max(
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
            * maximum_count
        )
    }


def analyse_edge(edge):
    upper_species = edge[
        "upper_species"
    ]

    lower_species = edge[
        "lower_species"
    ]

    stored_query = edge[
        "stored_query"
    ]

    stored_subject = edge[
        "stored_subject"
    ]

    full_pairs = read_full_anchor_pairs(
        Path(
            edge["anchor_file"]
        )
    )

    upper_genes = species_assets[
        upper_species
    ]["genes"]

    lower_genes = species_assets[
        lower_species
    ]["genes"]

    mapped_pairs = []
    unknown_query_genes = 0
    unknown_subject_genes = 0

    chromosome_pair_counts = Counter()

    for query_gene, subject_gene in full_pairs:
        if stored_query == upper_species:
            upper_gene = query_gene
            lower_gene = subject_gene
        elif stored_query == lower_species:
            upper_gene = subject_gene
            lower_gene = query_gene
        else:
            raise SystemExit(
                "ERROR: Stored comparison orientation "
                f"does not match {upper_species}-"
                f"{lower_species}"
            )

        if upper_gene not in upper_genes:
            unknown_query_genes += 1
            continue

        if lower_gene not in lower_genes:
            unknown_subject_genes += 1
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

        chromosome_pair_counts[
            (
                upper_chromosome,
                lower_chromosome,
            )
        ] += 1

        mapped_pairs.append(
            (
                upper_gene,
                lower_gene,
                upper_chromosome,
                lower_chromosome,
            )
        )

    total_pairs = len(
        full_pairs
    )

    mapped_pair_count = len(
        mapped_pairs
    )

    mapping_fraction = (
        mapped_pair_count
        / total_pairs
        if total_pairs
        else 0.0
    )

    if total_pairs == 0:
        raise SystemExit(
            "ERROR: Full anchor file contains no "
            f"gene pairs: {edge['anchor_file']}"
        )

    if mapped_pair_count == 0:
        raise SystemExit(
            "ERROR: Zero anchor pairs mapped for "
            f"{upper_species}-{lower_species}. "
            "This indicates an ID or orientation error."
        )

    if mapping_fraction < 0.80:
        raise SystemExit(
            "ERROR: Anchor mapping fraction below "
            f"80% for {upper_species}-{lower_species}: "
            f"{mapped_pair_count}/{total_pairs} "
            f"({mapping_fraction:.4f})."
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

    upper_strong_partners = {
        chromosome: strong_partner_set(
            counts
        )
        for chromosome, counts
        in upper_to_lower.items()
    }

    lower_strong_partners = {
        chromosome: strong_partner_set(
            counts
        )
        for chromosome, counts
        in lower_to_upper.items()
    }

    candidate_pairs = set()

    for upper_chromosome, partners in (
        upper_strong_partners.items()
    ):
        if len(partners) > 1:
            for lower_chromosome in partners:
                candidate_pairs.add(
                    (
                        upper_chromosome,
                        lower_chromosome,
                    )
                )

    for lower_chromosome, partners in (
        lower_strong_partners.items()
    ):
        if len(partners) > 1:
            for upper_chromosome in partners:
                candidate_pairs.add(
                    (
                        upper_chromosome,
                        lower_chromosome,
                    )
                )

    graph = defaultdict(set)

    for (
        upper_chromosome,
        lower_chromosome,
    ), count in chromosome_pair_counts.items():
        if (
            lower_chromosome
            in upper_strong_partners.get(
                upper_chromosome,
                set(),
            )
        ):
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

    visited = set()
    component_ratios = Counter()

    for start_node in graph:
        if start_node in visited:
            continue

        queue = deque(
            [start_node]
        )

        visited.add(
            start_node
        )

        upper_component = set()
        lower_component = set()

        while queue:
            node = queue.popleft()

            side, chromosome = node

            if side == "U":
                upper_component.add(
                    chromosome
                )
            else:
                lower_component.add(
                    chromosome
                )

            for neighbour in graph[node]:
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
            ratio_label = (
                f"{len(upper_component)}:"
                f"{len(lower_component)}"
            )

            component_weight = 0

            for upper_chromosome in (
                upper_component
            ):
                for lower_chromosome in (
                    lower_component
                ):
                    component_weight += (
                        chromosome_pair_counts.get(
                            (
                                upper_chromosome,
                                lower_chromosome,
                            ),
                            0,
                        )
                    )

            component_ratios[
                ratio_label
            ] += component_weight

    modal_ratio = (
        component_ratios.most_common(
            1
        )[0][0]
        if component_ratios
        else "NA"
    )

    return {
        "total_pairs": total_pairs,
        "mapped_pairs": mapped_pairs,
        "mapped_pair_count": (
            mapped_pair_count
        ),
        "mapping_fraction": (
            mapping_fraction
        ),
        "unknown_upper_genes": (
            unknown_query_genes
        ),
        "unknown_lower_genes": (
            unknown_subject_genes
        ),
        "chromosome_pair_counts": (
            chromosome_pair_counts
        ),
        "candidate_pairs": (
            candidate_pairs
        ),
        "modal_ratio": (
            modal_ratio
        ),
    }


def deterministic_subsample(
    pairs,
    maximum_pairs,
):
    if len(pairs) <= maximum_pairs:
        return pairs

    step = len(
        pairs
    ) / maximum_pairs

    sampled = []

    position = 0.0

    while (
        int(position) < len(pairs)
        and len(sampled) < maximum_pairs
    ):
        sampled.append(
            pairs[
                int(position)
            ]
        )

        position += step

    return sampled


def draw_macrosynteny(
    plot_name,
    species_order,
    edge_rows,
    output_stem,
):
    track_layout = (
        build_track_layout(
            species_order
        )
    )

    figure_height = max(
        10,
        1.55 * len(
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
        1.10,
    )

    axis.set_ylim(
        -0.55,
        len(species_order) - 0.45,
    )

    axis.axis(
        "off"
    )

    y_positions = {
        species_code: (
            len(species_order)
            - 1
            - index
        )
        for index, species_code
        in enumerate(species_order)
    }

    qc_rows = []
    ratio_rows = []
    rearrangement_rows = []

    edge_analyses = {}

    for edge in edge_rows:
        analysis = analyse_edge(
            edge
        )

        edge_analyses[
            edge["edge_rank"]
        ] = analysis

        qc_rows.append(
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
                "anchor_file": edge[
                    "anchor_file"
                ],
                "full_anchor_pairs": analysis[
                    "total_pairs"
                ],
                "mapped_anchor_pairs": analysis[
                    "mapped_pair_count"
                ],
                "mapping_fraction": (
                    f"{analysis['mapping_fraction']:.6f}"
                ),
                "unmapped_upper_gene_pairs": analysis[
                    "unknown_upper_genes"
                ],
                "unmapped_lower_gene_pairs": analysis[
                    "unknown_lower_genes"
                ],
                "status": "PASS",
            }
        )

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
                "full_anchor_pairs": analysis[
                    "total_pairs"
                ],
                "mapped_anchor_pairs": analysis[
                    "mapped_pair_count"
                ],
                "modal_synteny_ratio": analysis[
                    "modal_ratio"
                ],
                "minimum_strong_anchor_pairs": (
                    MIN_STRONG_ANCHORS
                ),
                "relative_strong_threshold": (
                    RELATIVE_STRONG_THRESHOLD
                ),
                "candidate_chromosome_pairs": len(
                    analysis[
                        "candidate_pairs"
                    ]
                ),
            }
        )

        for (
            upper_chromosome,
            lower_chromosome,
        ) in sorted(
            analysis[
                "candidate_pairs"
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
                    "comparison_id": edge[
                        "comparison_id"
                    ],
                    "upper_chromosome": (
                        upper_chromosome
                    ),
                    "lower_chromosome": (
                        lower_chromosome
                    ),
                    "full_anchor_pairs": (
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
                        "MULTIPLE_STRONG_CHROMOSOME_PARTNERS"
                    ),
                }
            )

    # Draw synteny ribbons before chromosome tracks.
    for edge in edge_rows:
        upper_species = edge[
            "upper_species"
        ]

        lower_species = edge[
            "lower_species"
        ]

        upper_y = y_positions[
            upper_species
        ]

        lower_y = y_positions[
            lower_species
        ]

        analysis = edge_analyses[
            edge["edge_rank"]
        ]

        plot_pairs = (
            deterministic_subsample(
                analysis[
                    "mapped_pairs"
                ],
                maximum_pairs=3500,
            )
        )

        for (
            upper_gene,
            lower_gene,
            upper_chromosome,
            lower_chromosome,
        ) in plot_pairs:
            upper_x = gene_x_position(
                upper_species,
                upper_gene,
                track_layout,
            )

            lower_x = gene_x_position(
                lower_species,
                lower_gene,
                track_layout,
            )

            if (
                upper_x is None
                or lower_x is None
            ):
                continue

            middle_y = (
                upper_y
                + lower_y
            ) / 2.0

            path_vertices = [
                (
                    upper_x,
                    upper_y - 0.035,
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
                    lower_y + 0.035,
                ),
            ]

            path_codes = [
                MplPath.MOVETO,
                MplPath.CURVE4,
                MplPath.CURVE4,
                MplPath.CURVE4,
            ]

            is_candidate = (
                (
                    upper_chromosome,
                    lower_chromosome,
                )
                in analysis[
                    "candidate_pairs"
                ]
            )

            if is_candidate:
                edge_colour = (
                    "#d73027"
                )

                edge_alpha = 0.12
                edge_width = 0.45
            else:
                edge_colour = (
                    "#4393c3"
                )

                edge_alpha = 0.075
                edge_width = 0.34

            axis.add_patch(
                PathPatch(
                    MplPath(
                        path_vertices,
                        path_codes,
                    ),
                    facecolor="none",
                    edgecolor=edge_colour,
                    linewidth=edge_width,
                    alpha=edge_alpha,
                    zorder=1,
                )
            )

        axis.text(
            1.015,
            (
                upper_y + lower_y
            ) / 2.0,
            (
                f"{upper_species}-"
                f"{lower_species}\n"
                f"ratio "
                f"{analysis['modal_ratio']}"
            ),
            fontsize=8.5,
            ha="left",
            va="center",
            color="#222222",
        )

    # Draw chromosomes and labels.
    for species_code in species_order:
        y_position = y_positions[
            species_code
        ]

        metadata = species_assets[
            species_code
        ]

        label_colour = (
            species_track_colour(
                species_code,
                metadata[
                    "ploidy"
                ],
            )
        )

        axis.text(
            -0.085,
            y_position,
            species_code,
            fontsize=14,
            fontweight="bold",
            ha="right",
            va="center",
            color=label_colour,
        )

        axis.text(
            -0.085,
            y_position - 0.20,
            metadata[
                "scientific_name"
            ],
            fontsize=7.5,
            fontstyle="italic",
            ha="right",
            va="center",
            color="#555555",
        )

        for chromosome in metadata[
            "ordered_chromosomes"
        ]:
            chromosome_layout = (
                track_layout[
                    species_code
                ][chromosome]
            )

            axis.plot(
                [
                    chromosome_layout[
                        "start"
                    ],
                    chromosome_layout[
                        "end"
                    ],
                ],
                [
                    y_position,
                    y_position,
                ],
                color="#111111",
                linewidth=5.0,
                solid_capstyle="butt",
                zorder=3,
            )

            axis.text(
                chromosome_layout[
                    "midpoint"
                ],
                y_position + 0.145,
                chromosome,
                fontsize=10.5,
                fontweight="bold",
                ha="center",
                va="bottom",
                color="#111111",
                zorder=4,
            )

    if plot_name == "FULL_9_SPECIES":
        title = (
            "Macrosynteny overview of nine species "
            "(tree-based order; PMAJ as outgroup)"
        )
    else:
        title = (
            "Macrosynteny overview of eight "
            "Veronica species (tree-based order)"
        )

    axis.set_title(
        title,
        fontsize=16,
        fontweight="bold",
        pad=20,
    )

    figure.text(
        0.012,
        0.012,
        (
            "Blue links: mapped full-anchor synteny pairs. "
            "Red links: chromosome pairs belonging to "
            "multiple strong-partner relationships. "
            "Ratios are operational chromosome-component "
            "ratios and require interpretation with Ks, "
            "self-synteny and syntenic-depth analyses."
        ),
        fontsize=8.2,
        ha="left",
        va="bottom",
    )

    figure.tight_layout(
        rect=(
            0.04,
            0.045,
            0.94,
            0.97,
        )
    )

    pdf_file = Path(
        str(output_stem)
        + ".pdf"
    )

    png_file = Path(
        str(output_stem)
        + ".png"
    )

    svg_file = Path(
        str(output_stem)
        + ".svg"
    )

    figure.savefig(
        pdf_file,
        bbox_inches="tight",
    )

    figure.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight",
    )

    figure.savefig(
        svg_file,
        bbox_inches="tight",
    )

    plt.close(
        figure
    )

    return {
        "qc_rows": qc_rows,
        "ratio_rows": ratio_rows,
        "rearrangement_rows": (
            rearrangement_rows
        ),
        "pdf": pdf_file,
        "png": png_file,
        "svg": svg_file,
    }


full_outputs = draw_macrosynteny(
    "FULL_9_SPECIES",
    full_order,
    full_edges,
    PLOT_DIR
    / "full_9_species.macrosynteny.refined",
)

veronica_outputs = draw_macrosynteny(
    "VERONICA_ONLY",
    veronica_order,
    veronica_edges,
    PLOT_DIR
    / "veronica_only_8_species.macrosynteny.refined",
)

qc_fieldnames = [
    "plot_name",
    "edge_rank",
    "upper_species",
    "lower_species",
    "comparison_id",
    "anchor_file",
    "full_anchor_pairs",
    "mapped_anchor_pairs",
    "mapping_fraction",
    "unmapped_upper_gene_pairs",
    "unmapped_lower_gene_pairs",
    "status",
]

ratio_fieldnames = [
    "plot_name",
    "edge_rank",
    "upper_species",
    "lower_species",
    "comparison_id",
    "full_anchor_pairs",
    "mapped_anchor_pairs",
    "modal_synteny_ratio",
    "minimum_strong_anchor_pairs",
    "relative_strong_threshold",
    "candidate_chromosome_pairs",
]

rearrangement_fieldnames = [
    "plot_name",
    "edge_rank",
    "upper_species",
    "lower_species",
    "comparison_id",
    "upper_chromosome",
    "lower_chromosome",
    "full_anchor_pairs",
    "classification",
]

write_tsv(
    TABLE_DIR
    / "refined_full_plot_edge_mapping_qc.tsv",
    full_outputs[
        "qc_rows"
    ],
    qc_fieldnames,
)

write_tsv(
    TABLE_DIR
    / "refined_veronica_only_edge_mapping_qc.tsv",
    veronica_outputs[
        "qc_rows"
    ],
    qc_fieldnames,
)

write_tsv(
    TABLE_DIR
    / "refined_full_plot_synteny_ratios.tsv",
    full_outputs[
        "ratio_rows"
    ],
    ratio_fieldnames,
)

write_tsv(
    TABLE_DIR
    / "refined_veronica_only_plot_synteny_ratios.tsv",
    veronica_outputs[
        "ratio_rows"
    ],
    ratio_fieldnames,
)

write_tsv(
    TABLE_DIR
    / "refined_full_plot_rearrangement_candidates.tsv",
    full_outputs[
        "rearrangement_rows"
    ],
    rearrangement_fieldnames,
)

write_tsv(
    TABLE_DIR
    / "refined_veronica_only_plot_rearrangement_candidates.tsv",
    veronica_outputs[
        "rearrangement_rows"
    ],
    rearrangement_fieldnames,
)

summary_rows = [
    {
        "metric": "full_plot_species_count",
        "value": len(
            full_order
        ),
    },
    {
        "metric": "veronica_only_species_count",
        "value": len(
            veronica_order
        ),
    },
    {
        "metric": "full_plot_edges",
        "value": len(
            full_edges
        ),
    },
    {
        "metric": "veronica_only_edges",
        "value": len(
            veronica_edges
        ),
    },
    {
        "metric": "full_plot_total_anchor_pairs",
        "value": sum(
            int(
                row[
                    "full_anchor_pairs"
                ]
            )
            for row in full_outputs[
                "qc_rows"
            ]
        ),
    },
    {
        "metric": "full_plot_mapped_anchor_pairs",
        "value": sum(
            int(
                row[
                    "mapped_anchor_pairs"
                ]
            )
            for row in full_outputs[
                "qc_rows"
            ]
        ),
    },
    {
        "metric": "veronica_plot_total_anchor_pairs",
        "value": sum(
            int(
                row[
                    "full_anchor_pairs"
                ]
            )
            for row in veronica_outputs[
                "qc_rows"
            ]
        ),
    },
    {
        "metric": "veronica_plot_mapped_anchor_pairs",
        "value": sum(
            int(
                row[
                    "mapped_anchor_pairs"
                ]
            )
            for row in veronica_outputs[
                "qc_rows"
            ]
        ),
    },
    {
        "metric": "full_plot_rearrangement_rows",
        "value": len(
            full_outputs[
                "rearrangement_rows"
            ]
        ),
    },
    {
        "metric": "veronica_only_rearrangement_rows",
        "value": len(
            veronica_outputs[
                "rearrangement_rows"
            ]
        ),
    },
    {
        "metric": "full_plot_pdf",
        "value": str(
            full_outputs[
                "pdf"
            ]
        ),
    },
    {
        "metric": "full_plot_png",
        "value": str(
            full_outputs[
                "png"
            ]
        ),
    },
    {
        "metric": "full_plot_svg",
        "value": str(
            full_outputs[
                "svg"
            ]
        ),
    },
    {
        "metric": "veronica_plot_pdf",
        "value": str(
            veronica_outputs[
                "pdf"
            ]
        ),
    },
    {
        "metric": "veronica_plot_png",
        "value": str(
            veronica_outputs[
                "png"
            ]
        ),
    },
    {
        "metric": "veronica_plot_svg",
        "value": str(
            veronica_outputs[
                "svg"
            ]
        ),
    },
    {
        "metric": "status",
        "value": "PASS",
    },
]

write_tsv(
    TABLE_DIR
    / "refined_macrosynteny_plot_summary.tsv",
    summary_rows,
    [
        "metric",
        "value",
    ],
)

print(
    "Full order:",
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
    "Full mapped anchor pairs:",
    sum(
        int(
            row[
                "mapped_anchor_pairs"
            ]
        )
        for row in full_outputs[
            "qc_rows"
        ]
    ),
)

print(
    "Veronica-only mapped anchor pairs:",
    sum(
        int(
            row[
                "mapped_anchor_pairs"
            ]
        )
        for row in veronica_outputs[
            "qc_rows"
        ]
    ),
)

print(
    "Corrected macrosynteny plotting: PASS"
)
PY

for FILE in \
    "${PLOT_DIR}/full_9_species.macrosynteny.refined.pdf" \
    "${PLOT_DIR}/full_9_species.macrosynteny.refined.png" \
    "${PLOT_DIR}/full_9_species.macrosynteny.refined.svg" \
    "${PLOT_DIR}/veronica_only_8_species.macrosynteny.refined.pdf" \
    "${PLOT_DIR}/veronica_only_8_species.macrosynteny.refined.png" \
    "${PLOT_DIR}/veronica_only_8_species.macrosynteny.refined.svg" \
    "${TABLE_DIR}/refined_full_plot_edge_mapping_qc.tsv" \
    "${TABLE_DIR}/refined_veronica_only_edge_mapping_qc.tsv" \
    "${TABLE_DIR}/refined_full_plot_synteny_ratios.tsv" \
    "${TABLE_DIR}/refined_veronica_only_plot_synteny_ratios.tsv" \
    "${TABLE_DIR}/refined_full_plot_rearrangement_candidates.tsv" \
    "${TABLE_DIR}/refined_veronica_only_plot_rearrangement_candidates.tsv" \
    "${TABLE_DIR}/refined_macrosynteny_plot_summary.tsv"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Expected output missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

echo
echo "============================================================"
echo "Full plot anchor mapping QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/refined_full_plot_edge_mapping_qc.tsv"

echo
echo "============================================================"
echo "Veronica-only anchor mapping QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/refined_veronica_only_edge_mapping_qc.tsv"

echo
echo "============================================================"
echo "Synteny-ratio results"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/refined_full_plot_synteny_ratios.tsv"

echo
echo "============================================================"
echo "Plot summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/refined_macrosynteny_plot_summary.tsv"

cat > "${CHECKPOINT_DIR}/REFINED_MACROSYNTENY_PLOTS_COMPLETE.txt" <<EOF2
checkpoint=refined_macrosynteny_plots_full_anchor_corrected
date=$(date --iso-8601=seconds)
status=PASS
anchor_source=full_jcvi_anchor_files
minimum_mapping_fraction=0.80
next_step=inspect_corrected_figures_and_interpret_synteny_depth
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
echo "Corrected Step 34H completed successfully."
echo "Checkpoint:"
cat \
    "${CHECKPOINT_DIR}/REFINED_MACROSYNTENY_PLOTS_COMPLETE.txt"
