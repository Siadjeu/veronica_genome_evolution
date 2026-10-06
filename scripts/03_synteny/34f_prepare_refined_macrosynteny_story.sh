#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=macro_refine
#SBATCH --output=10_synteny/logs/macro_refine_%j.out
#SBATCH --error=10_synteny/logs/macro_refine_%j.err

set -euo pipefail

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

TREE_FILE="${PROJECT_DIR}/09_orthology/checkpoint_duplication_node_mapping/SpeciesTree_explicitly_rooted_with_PMAJ.nwk"

STEP34_RESULTS="${SYNTENY_DIR}/jcvi_step34/tables/jcvi_step34_results.tsv"
STEP34_CHECKPOINT="${SYNTENY_DIR}/checkpoint_jcvi_step34/JCVI_STEP34_COMPLETE.txt"
INPUT_MANIFEST="${SYNTENY_DIR}/manifests/jcvi_input_manifest.tsv"

REFINED_DIR="${SYNTENY_DIR}/refined_macrosynteny"
TABLE_DIR="${REFINED_DIR}/tables"
ORDER_DIR="${REFINED_DIR}/orders"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_refined_macrosynteny_story"

mkdir -p \
    "${TABLE_DIR}" \
    "${ORDER_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${SYNTENY_DIR}/logs"

cd "${PROJECT_DIR}"

# ============================================================
# Validate required inputs
# ============================================================

for FILE in \
    "${TREE_FILE}" \
    "${STEP34_RESULTS}" \
    "${STEP34_CHECKPOINT}" \
    "${INPUT_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Required input is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

if ! grep -q '^status=PASS$' "${STEP34_CHECKPOINT}"; then
    echo "ERROR: Step 34 checkpoint is not PASS." >&2
    cat "${STEP34_CHECKPOINT}" >&2
    exit 1
fi

# ============================================================
# Activate jcvi_env
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

# ============================================================
# Clean only Step 34F outputs
# ============================================================

rm -f "${TABLE_DIR}"/*
rm -f "${ORDER_DIR}"/*
rm -f "${CHECKPOINT_DIR}"/*

# ============================================================
# Build refined storyline tables
# ============================================================

python - \
    "${TREE_FILE}" \
    "${STEP34_RESULTS}" \
    "${INPUT_MANIFEST}" \
    "${TABLE_DIR}" \
    "${ORDER_DIR}" <<'PY'
from __future__ import annotations

import csv
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

(
    tree_file_name,
    step34_results_name,
    input_manifest_name,
    table_dir_name,
    order_dir_name,
) = sys.argv[1:]

tree_file = Path(tree_file_name)
step34_results_file = Path(step34_results_name)
input_manifest_file = Path(input_manifest_name)
table_dir = Path(table_dir_name)
order_dir = Path(order_dir_name)

table_dir.mkdir(parents=True, exist_ok=True)
order_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# Species metadata
# ============================================================

species_meta = {
    "PMAJ": {
        "scientific_name": "Plantago major",
        "ploidy": "2x",
        "role": "outgroup",
    },
    "VPAN": {
        "scientific_name": "Veronica panormitana",
        "ploidy": "2x",
        "role": "Veronica",
    },
    "VSCU": {
        "scientific_name": "Veronica scutellata",
        "ploidy": "2x",
        "role": "Veronica",
    },
    "VANA": {
        "scientific_name": "Veronica anagallis-aquatica",
        "ploidy": "4x",
        "role": "Veronica",
    },
    "VARV": {
        "scientific_name": "Veronica arvensis",
        "ploidy": "2x",
        "role": "Veronica",
    },
    "VPER": {
        "scientific_name": "Veronica persica",
        "ploidy": "4x",
        "role": "Veronica",
    },
    "VSER": {
        "scientific_name": "Veronica serpyllifolia",
        "ploidy": "2x",
        "role": "Veronica",
    },
    "VTRI": {
        "scientific_name": "Veronica triphyllos",
        "ploidy": "2x",
        "role": "Veronica",
    },
    "VVER": {
        "scientific_name": "Veronica verna",
        "ploidy": "2x",
        "role": "Veronica",
    },
}

expected_codes = list(species_meta)

veronica_codes = {
    code
    for code, metadata in species_meta.items()
    if metadata["role"] == "Veronica"
}

# ============================================================
# Parse species order directly from rooted Newick text
#
# Internal OrthoFinder node names such as N3 are ignored.
# Species are ordered according to the first occurrence of
# their fixed project abbreviations in the Newick string.
# ============================================================

newick_text = tree_file.read_text(
    encoding="utf-8",
    errors="replace",
).strip()

if not newick_text:
    raise SystemExit(
        f"ERROR: Rooted species tree is empty: {tree_file}"
    )

occurrences = []

for code in expected_codes:
    matches = list(
        re.finditer(
            rf"(?<![A-Za-z0-9]){re.escape(code)}(?![A-Za-z0-9])",
            newick_text,
            flags=re.IGNORECASE,
        )
    )

    if not matches:
        matches = list(
            re.finditer(
                re.escape(code),
                newick_text,
                flags=re.IGNORECASE,
            )
        )

    if matches:
        occurrences.append(
            (
                matches[0].start(),
                code,
            )
        )

occurrences.sort(
    key=lambda item: item[0]
)

recognized_order = [
    code
    for _, code in occurrences
]

missing_species = sorted(
    set(expected_codes)
    - set(recognized_order)
)

if missing_species:
    raise SystemExit(
        "ERROR: Species missing from rooted tree: "
        + ",".join(missing_species)
        + f"; tree={tree_file}"
    )

if len(recognized_order) != 9:
    raise SystemExit(
        f"ERROR: Expected 9 species in tree order; "
        f"observed {len(recognized_order)}: "
        + ",".join(recognized_order)
    )

veronica_tree_order = [
    code
    for code in recognized_order
    if code in veronica_codes
]

if len(veronica_tree_order) != 8:
    raise SystemExit(
        f"ERROR: Expected 8 Veronica species; "
        f"observed {len(veronica_tree_order)}."
    )

# Keep the outgroup at the beginning of the full plot.
full_tree_order = [
    "PMAJ"
] + [
    code
    for code in veronica_tree_order
]

print(
    "Rooted tree:",
    tree_file,
)

print(
    "Recognized tree order:",
    ",".join(recognized_order),
)

print(
    "Full plot order:",
    ",".join(full_tree_order),
)

print(
    "Veronica-only order:",
    ",".join(veronica_tree_order),
)

# ============================================================
# Read validated JCVI inputs
# ============================================================

with input_manifest_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    input_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

inputs = {
    row["species_code"]: row
    for row in input_rows
}

for code in expected_codes:
    if code not in inputs:
        raise SystemExit(
            f"ERROR: Species {code} is absent from "
            "jcvi_input_manifest.tsv."
        )

    bed_file = Path(
        inputs[code]["bed_file"]
    )

    if not bed_file.is_file():
        raise SystemExit(
            f"ERROR: BED file is missing for {code}: "
            f"{bed_file}"
        )

# ============================================================
# Read Step 34 comparisons
# ============================================================

with step34_results_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    result_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

pass_rows = [
    row
    for row in result_rows
    if row["status"] == "PASS"
]

pair_index = {}

for row in pass_rows:
    query = row["query_species"]
    subject = row["subject_species"]

    pair_index[(query, subject)] = row


def resolve_pair(species_a, species_b):
    if (species_a, species_b) in pair_index:
        row = pair_index[
            (species_a, species_b)
        ]

        return {
            "comparison_id": row["comparison_id"],
            "stored_query": row["query_species"],
            "stored_subject": row["subject_species"],
            "orientation": "FORWARD",
            "anchor_file": row["anchor_file"],
            "simple_anchor_file": row[
                "simple_anchor_file"
            ],
            "status": "AVAILABLE",
        }

    if (species_b, species_a) in pair_index:
        row = pair_index[
            (species_b, species_a)
        ]

        return {
            "comparison_id": row["comparison_id"],
            "stored_query": row["query_species"],
            "stored_subject": row["subject_species"],
            "orientation": "REVERSE_AVAILABLE",
            "anchor_file": row["anchor_file"],
            "simple_anchor_file": row[
                "simple_anchor_file"
            ],
            "status": "AVAILABLE",
        }

    return {
        "comparison_id": (
            f"{species_a}__{species_b}"
        ),
        "stored_query": species_a,
        "stored_subject": species_b,
        "orientation": "NONE",
        "anchor_file": "",
        "simple_anchor_file": "",
        "status": "MISSING",
    }

# ============================================================
# Write plot-order tables
# ============================================================

def write_order_table(
    path,
    plot_name,
    species_order,
):
    rows = []

    for rank, code in enumerate(
        species_order,
        start=1,
    ):
        rows.append(
            {
                "plot_name": plot_name,
                "plot_rank": rank,
                "species_code": code,
                "display_label": code,
                "scientific_name": species_meta[
                    code
                ]["scientific_name"],
                "ploidy": species_meta[
                    code
                ]["ploidy"],
                "role": species_meta[
                    code
                ]["role"],
                "highlight_group": (
                    "TETRAPLOID"
                    if species_meta[code][
                        "ploidy"
                    ] == "4x"
                    else (
                        "OUTGROUP"
                        if code == "PMAJ"
                        else "DIPLOID"
                    )
                ),
            }
        )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                rows[0].keys()
            ),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


write_order_table(
    order_dir
    / "full_9_species_plot_order.tsv",
    "FULL_9_SPECIES",
    full_tree_order,
)

write_order_table(
    order_dir
    / "veronica_only_plot_order.tsv",
    "VERONICA_ONLY",
    veronica_tree_order,
)

# ============================================================
# Determine required adjacent comparisons
# ============================================================

def adjacent_comparisons(
    plot_name,
    species_order,
):
    rows = []

    for edge_rank in range(
        1,
        len(species_order),
    ):
        upper = species_order[
            edge_rank - 1
        ]

        lower = species_order[
            edge_rank
        ]

        resolved = resolve_pair(
            upper,
            lower,
        )

        rows.append(
            {
                "plot_name": plot_name,
                "edge_rank": edge_rank,
                "upper_species": upper,
                "lower_species": lower,
                "comparison_id": resolved[
                    "comparison_id"
                ],
                "stored_query": resolved[
                    "stored_query"
                ],
                "stored_subject": resolved[
                    "stored_subject"
                ],
                "orientation": resolved[
                    "orientation"
                ],
                "anchor_file": resolved[
                    "anchor_file"
                ],
                "simple_anchor_file": resolved[
                    "simple_anchor_file"
                ],
                "status": resolved[
                    "status"
                ],
            }
        )

    return rows


full_edges = adjacent_comparisons(
    "FULL_9_SPECIES",
    full_tree_order,
)

veronica_edges = adjacent_comparisons(
    "VERONICA_ONLY",
    veronica_tree_order,
)


def write_rows(path, rows, fieldnames=None):
    if fieldnames is None:
        fieldnames = list(
            rows[0].keys()
        )

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


write_rows(
    table_dir
    / "full_9_species_required_adjacent_comparisons.tsv",
    full_edges,
)

write_rows(
    table_dir
    / "veronica_only_required_adjacent_comparisons.tsv",
    veronica_edges,
)

missing_edges = []
seen_missing = set()

for row in full_edges + veronica_edges:
    if row["status"] != "MISSING":
        continue

    key = tuple(
        sorted(
            [
                row["upper_species"],
                row["lower_species"],
            ]
        )
    )

    if key in seen_missing:
        continue

    seen_missing.add(key)
    missing_edges.append(row)

edge_fields = list(
    full_edges[0].keys()
)

write_rows(
    table_dir
    / "missing_tree_order_adjacent_comparisons.tsv",
    missing_edges,
    fieldnames=edge_fields,
)

# ============================================================
# Read BED files
# ============================================================

def read_bed(path):
    gene_to_chromosome = {}
    chromosome_genes = defaultdict(set)

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
            gene_id = fields[3]

            gene_to_chromosome[
                gene_id
            ] = chromosome

            chromosome_genes[
                chromosome
            ].add(gene_id)

    return (
        gene_to_chromosome,
        chromosome_genes,
    )


bed_data = {}

for code in expected_codes:
    bed_data[code] = read_bed(
        Path(
            inputs[code]["bed_file"]
        )
    )

# ============================================================
# Parse anchor files and calculate synteny ratios
#
# A strong chromosome partner must contain:
# - at least 20 anchor pairs; and
# - at least 25% of the strongest chromosome-pair count for
#   the focal chromosome.
# ============================================================

MIN_ANCHORS = 20
RELATIVE_THRESHOLD = 0.25


def parse_anchor_file(path):
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

            if len(fields) >= 2:
                pairs.append(
                    (
                        fields[0],
                        fields[1],
                    )
                )

    return pairs


def ploidy_number(label):
    match = re.match(
        r"^(\d+)x$",
        label,
    )

    if not match:
        return 0

    return int(
        match.group(1)
    )


def strong_partners(
    partner_counts,
):
    if not partner_counts:
        return []

    maximum = max(
        partner_counts.values()
    )

    return sorted(
        [
            partner
            for partner, count
            in partner_counts.items()
            if (
                count >= MIN_ANCHORS
                and count
                >= RELATIVE_THRESHOLD
                * maximum
            )
        ]
    )


chromosome_pair_rows = []
ratio_rows = []
rearrangement_rows = []

for row in pass_rows:
    if row["comparison_type"] == "self":
        continue

    comparison_id = row[
        "comparison_id"
    ]
    query = row[
        "query_species"
    ]
    subject = row[
        "subject_species"
    ]

    anchor_file = Path(
        row["anchor_file"]
    )

    if (
        not anchor_file.is_file()
        or anchor_file.stat().st_size == 0
    ):
        raise SystemExit(
            f"ERROR: Anchor file is missing "
            f"for {comparison_id}: "
            f"{anchor_file}"
        )

    query_gene_map = bed_data[
        query
    ][0]

    subject_gene_map = bed_data[
        subject
    ][0]

    chromosome_pair_counts = Counter()

    unknown_pairs = 0

    for query_gene, subject_gene in (
        parse_anchor_file(anchor_file)
    ):
        query_chromosome = (
            query_gene_map.get(
                query_gene
            )
        )

        subject_chromosome = (
            subject_gene_map.get(
                subject_gene
            )
        )

        if (
            query_chromosome is None
            or subject_chromosome is None
        ):
            unknown_pairs += 1
            continue

        chromosome_pair_counts[
            (
                query_chromosome,
                subject_chromosome,
            )
        ] += 1

    query_to_subject = defaultdict(
        Counter
    )

    subject_to_query = defaultdict(
        Counter
    )

    for (
        query_chromosome,
        subject_chromosome,
    ), count in chromosome_pair_counts.items():
        query_to_subject[
            query_chromosome
        ][subject_chromosome] = count

        subject_to_query[
            subject_chromosome
        ][query_chromosome] = count

        chromosome_pair_rows.append(
            {
                "comparison_id": comparison_id,
                "query_species": query,
                "subject_species": subject,
                "query_chromosome": (
                    query_chromosome
                ),
                "subject_chromosome": (
                    subject_chromosome
                ),
                "anchor_pairs": count,
            }
        )

    query_partner_depths = {
        chromosome: len(
            strong_partners(
                counts
            )
        )
        for chromosome, counts
        in query_to_subject.items()
    }

    subject_partner_depths = {
        chromosome: len(
            strong_partners(
                counts
            )
        )
        for chromosome, counts
        in subject_to_query.items()
    }

    query_ploidy = species_meta[
        query
    ]["ploidy"]

    subject_ploidy = species_meta[
        subject
    ]["ploidy"]

    query_ploidy_number = (
        ploidy_number(
            query_ploidy
        )
    )

    subject_ploidy_number = (
        ploidy_number(
            subject_ploidy
        )
    )

    if query_partner_depths:
        query_median_depth = (
            statistics.median(
                query_partner_depths.values()
            )
        )
    else:
        query_median_depth = 0

    if subject_partner_depths:
        subject_median_depth = (
            statistics.median(
                subject_partner_depths.values()
            )
        )
    else:
        subject_median_depth = 0

    if (
        query_ploidy_number
        < subject_ploidy_number
    ):
        smaller_species = query
        larger_species = subject
        smaller_to_larger = (
            query_median_depth
        )
        larger_to_smaller = (
            subject_median_depth
        )

    elif (
        subject_ploidy_number
        < query_ploidy_number
    ):
        smaller_species = subject
        larger_species = query
        smaller_to_larger = (
            subject_median_depth
        )
        larger_to_smaller = (
            query_median_depth
        )

    else:
        smaller_species = ""
        larger_species = ""
        smaller_to_larger = ""
        larger_to_smaller = ""

    if (
        smaller_to_larger != ""
        and larger_to_smaller != ""
    ):
        interpreted_ratio = (
            f"{smaller_to_larger}:"
            f"{larger_to_smaller}"
        )
    else:
        interpreted_ratio = ""

    ratio_rows.append(
        {
            "comparison_id": comparison_id,
            "query_species": query,
            "subject_species": subject,
            "query_ploidy": query_ploidy,
            "subject_ploidy": subject_ploidy,
            "query_median_strong_partner_depth": (
                query_median_depth
            ),
            "subject_median_strong_partner_depth": (
                subject_median_depth
            ),
            "smaller_ploidy_species": (
                smaller_species
            ),
            "larger_ploidy_species": (
                larger_species
            ),
            "median_larger_chromosomes_per_smaller_chromosome": (
                smaller_to_larger
            ),
            "median_smaller_chromosomes_per_larger_chromosome": (
                larger_to_smaller
            ),
            "interpreted_synteny_ratio": (
                interpreted_ratio
            ),
            "minimum_anchor_pairs": (
                MIN_ANCHORS
            ),
            "relative_strength_threshold": (
                RELATIVE_THRESHOLD
            ),
            "unknown_anchor_pairs": (
                unknown_pairs
            ),
        }
    )

    for focal_chromosome, counts in (
        query_to_subject.items()
    ):
        partners = strong_partners(
            counts
        )

        if len(partners) > 1:
            rearrangement_rows.append(
                {
                    "comparison_id": (
                        comparison_id
                    ),
                    "focal_species": query,
                    "partner_species": (
                        subject
                    ),
                    "focal_chromosome": (
                        focal_chromosome
                    ),
                    "strong_partner_chromosomes": (
                        ",".join(partners)
                    ),
                    "strong_partner_count": (
                        len(partners)
                    ),
                    "interpretation": (
                        "POTENTIAL_SPLIT_FUSION_TRANSLOCATION_OR_DUPLICATION"
                    ),
                }
            )

    for focal_chromosome, counts in (
        subject_to_query.items()
    ):
        partners = strong_partners(
            counts
        )

        if len(partners) > 1:
            rearrangement_rows.append(
                {
                    "comparison_id": (
                        comparison_id
                    ),
                    "focal_species": subject,
                    "partner_species": query,
                    "focal_chromosome": (
                        focal_chromosome
                    ),
                    "strong_partner_chromosomes": (
                        ",".join(partners)
                    ),
                    "strong_partner_count": (
                        len(partners)
                    ),
                    "interpretation": (
                        "POTENTIAL_SPLIT_FUSION_TRANSLOCATION_OR_DUPLICATION"
                    ),
                }
            )

write_rows(
    table_dir
    / "chromosome_pair_anchor_counts.tsv",
    chromosome_pair_rows,
    fieldnames=[
        "comparison_id",
        "query_species",
        "subject_species",
        "query_chromosome",
        "subject_chromosome",
        "anchor_pairs",
    ],
)

write_rows(
    table_dir
    / "synteny_ratio_summary.tsv",
    ratio_rows,
    fieldnames=[
        "comparison_id",
        "query_species",
        "subject_species",
        "query_ploidy",
        "subject_ploidy",
        "query_median_strong_partner_depth",
        "subject_median_strong_partner_depth",
        "smaller_ploidy_species",
        "larger_ploidy_species",
        "median_larger_chromosomes_per_smaller_chromosome",
        "median_smaller_chromosomes_per_larger_chromosome",
        "interpreted_synteny_ratio",
        "minimum_anchor_pairs",
        "relative_strength_threshold",
        "unknown_anchor_pairs",
    ],
)

write_rows(
    table_dir
    / "chromosome_rearrangement_candidates.tsv",
    rearrangement_rows,
    fieldnames=[
        "comparison_id",
        "focal_species",
        "partner_species",
        "focal_chromosome",
        "strong_partner_chromosomes",
        "strong_partner_count",
        "interpretation",
    ],
)

# ============================================================
# Write compact summary
# ============================================================

summary_rows = [
    {
        "metric": "rooted_species_tree",
        "value": str(tree_file),
    },
    {
        "metric": "recognized_tree_order",
        "value": ",".join(
            recognized_order
        ),
    },
    {
        "metric": "full_plot_order",
        "value": ",".join(
            full_tree_order
        ),
    },
    {
        "metric": "veronica_only_order",
        "value": ",".join(
            veronica_tree_order
        ),
    },
    {
        "metric": "full_plot_missing_adjacent_comparisons",
        "value": str(
            sum(
                row["status"] == "MISSING"
                for row in full_edges
            )
        ),
    },
    {
        "metric": "veronica_plot_missing_adjacent_comparisons",
        "value": str(
            sum(
                row["status"] == "MISSING"
                for row in veronica_edges
            )
        ),
    },
    {
        "metric": "unique_missing_comparisons",
        "value": str(
            len(missing_edges)
        ),
    },
    {
        "metric": "nonself_synteny_ratio_comparisons",
        "value": str(
            len(ratio_rows)
        ),
    },
    {
        "metric": "rearrangement_candidate_rows",
        "value": str(
            len(rearrangement_rows)
        ),
    },
    {
        "metric": "status",
        "value": "PASS",
    },
]

write_rows(
    table_dir
    / "refined_storyline_summary.tsv",
    summary_rows,
    fieldnames=[
        "metric",
        "value",
    ],
)

print(
    "Missing adjacent comparisons:",
    len(missing_edges),
)

print(
    "Synteny-ratio comparisons:",
    len(ratio_rows),
)

print(
    "Rearrangement candidates:",
    len(rearrangement_rows),
)

print(
    "Step 34F analysis: PASS"
)
PY

# ============================================================
# Validate outputs
# ============================================================

for FILE in \
    "${ORDER_DIR}/full_9_species_plot_order.tsv" \
    "${ORDER_DIR}/veronica_only_plot_order.tsv" \
    "${TABLE_DIR}/full_9_species_required_adjacent_comparisons.tsv" \
    "${TABLE_DIR}/veronica_only_required_adjacent_comparisons.tsv" \
    "${TABLE_DIR}/missing_tree_order_adjacent_comparisons.tsv" \
    "${TABLE_DIR}/chromosome_pair_anchor_counts.tsv" \
    "${TABLE_DIR}/synteny_ratio_summary.tsv" \
    "${TABLE_DIR}/chromosome_rearrangement_candidates.tsv" \
    "${TABLE_DIR}/refined_storyline_summary.tsv"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Expected output is missing or empty:" >&2
        echo "${FILE}" >&2
        exit 1
    fi
done

# ============================================================
# Display outputs
# ============================================================

echo
echo "============================================================"
echo "Refined storyline summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/refined_storyline_summary.tsv"

echo
echo "============================================================"
echo "Full 9-species phylogenetic plot order"
echo "============================================================"

column -t -s $'\t' \
    "${ORDER_DIR}/full_9_species_plot_order.tsv"

echo
echo "============================================================"
echo "Veronica-only phylogenetic plot order"
echo "============================================================"

column -t -s $'\t' \
    "${ORDER_DIR}/veronica_only_plot_order.tsv"

echo
echo "============================================================"
echo "Missing adjacent comparisons"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/missing_tree_order_adjacent_comparisons.tsv"

echo
echo "============================================================"
echo "Synteny-ratio summary"
echo "============================================================"

column -t -s $'\t' \
    "${TABLE_DIR}/synteny_ratio_summary.tsv"

echo
echo "============================================================"
echo "First 50 rearrangement candidates"
echo "============================================================"

head -n 51 \
    "${TABLE_DIR}/chromosome_rearrangement_candidates.tsv" \
    | column -t -s $'\t'

# ============================================================
# Checkpoint
# ============================================================

cp -f \
    "${TABLE_DIR}/refined_storyline_summary.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${ORDER_DIR}/full_9_species_plot_order.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${ORDER_DIR}/veronica_only_plot_order.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/missing_tree_order_adjacent_comparisons.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/synteny_ratio_summary.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLE_DIR}/chromosome_rearrangement_candidates.tsv" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/REFINED_MACROSYNTENY_STORY_COMPLETE.txt" <<EOF2
checkpoint=refined_macrosynteny_story
date=$(date --iso-8601=seconds)
status=PASS
next_step=run_missing_tree_order_comparisons_or_generate_refined_plots
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
echo "============================================================"
echo "Step 34F completed successfully"
echo "============================================================"
echo "Checkpoint:"
echo "${CHECKPOINT_DIR}/REFINED_MACROSYNTENY_STORY_COMPLETE.txt"
