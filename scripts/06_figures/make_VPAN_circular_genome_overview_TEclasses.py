#!/usr/bin/env python3

"""
VPAN circular genome overview.

Tracks from outside to inside:
    1. Gene count
    2. LTR/Gypsy proportion
    3. LTR/Copia proportion
    4. DNA-transposon proportion
    5. GC content

All statistics are calculated in non-overlapping 1-Mb windows.

RepeatMasker classes:
    Gypsy:
        LTR/Gypsy

    Copia:
        LTR/Copia

    DNA transposons:
        DNA
        DNA/*

RC/Helitron is intentionally not included in the DNA track.

Outputs:
    <prefix>.window_statistics.tsv
    <prefix>.pdf
    <prefix>.svg
    <prefix>.png
"""

from __future__ import annotations

import argparse
import gzip
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


###############################################################################
# Colours
###############################################################################

GENE_COLOR = "#2E86C1"
GYPSY_COLOR = "#8E44AD"
COPIA_COLOR = "#F39C12"
DNA_COLOR = "#27AE60"
GC_COLOR = "#D73027"

CHR_COLOR = "#F2F2F2"
CHR_EDGE = "#222222"


###############################################################################
# Radial layout
###############################################################################

CHR_BOTTOM = 0.91
CHR_HEIGHT = 0.075

GENE_BOTTOM = 0.74
GENE_HEIGHT = 0.14

GYPSY_BOTTOM = 0.59
GYPSY_HEIGHT = 0.11

COPIA_BOTTOM = 0.45
COPIA_HEIGHT = 0.11

DNA_BOTTOM = 0.31
DNA_HEIGHT = 0.11

GC_BOTTOM = 0.16
GC_HEIGHT = 0.11


###############################################################################
# Text input
###############################################################################

def open_text(path: Path):

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


###############################################################################
# Chromosome mapping
###############################################################################

def read_chromosome_map(path: Path):

    repeat_to_genome = {}
    genome_to_display = {}

    with path.open(
        "r",
        encoding="utf-8",
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

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 2:

                raise SystemExit(
                    f"ERROR: chromosome-map line {line_number} "
                    "does not contain exactly two tab-separated columns."
                )

            genome_id = fields[0].strip()
            repeat_id = fields[1].strip()

            if (
                genome_id == "genome_id"
                and repeat_id == "repeatmasker_id"
            ):
                continue

            repeat_to_genome[
                repeat_id
            ] = genome_id

            genome_to_display[
                genome_id
            ] = repeat_id.lower()

    return repeat_to_genome, genome_to_display


###############################################################################
# FASTA
###############################################################################

def fasta_lengths(path: Path):

    lengths = {}

    current = None
    length = 0

    with open_text(path) as handle:

        for line in handle:

            if line.startswith(">"):

                if current is not None:
                    lengths[current] = length

                current = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

                length = 0

            else:

                length += len(
                    line.strip()
                )

    if current is not None:
        lengths[current] = length

    return lengths


def load_selected_sequences(
    path: Path,
    selected: set[str],
):

    chunks = defaultdict(list)
    current = None

    with open_text(path) as handle:

        for line in handle:

            if line.startswith(">"):

                current = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

            elif current in selected:

                chunks[current].append(
                    line.strip().upper()
                )

    sequences = {
        seqid: "".join(parts)
        for seqid, parts
        in chunks.items()
    }

    missing = (
        selected
        - set(sequences)
    )

    if missing:

        raise SystemExit(
            "ERROR: Missing FASTA sequences:\n"
            + "\n".join(
                sorted(missing)
            )
        )

    return sequences


###############################################################################
# Genes
###############################################################################

def parse_gene_midpoints(
    gff_path: Path,
    selected: set[str],
):

    genes = defaultdict(list)

    with open_text(gff_path) as handle:

        for line in handle:

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                continue

            seqid = fields[0]
            feature = fields[2]

            if seqid not in selected:
                continue

            if feature.lower() != "gene":
                continue

            try:
                start = int(fields[3]) - 1
                end = int(fields[4])

            except ValueError:
                continue

            genes[
                seqid
            ].append(
                (
                    start + end
                ) // 2
            )

    return genes


###############################################################################
# RepeatMasker .out parser
###############################################################################

def classify_repeat(
    repeat_class: str,
):

    if repeat_class == "LTR/Gypsy":
        return "gypsy"

    if repeat_class == "LTR/Copia":
        return "copia"

    if (
        repeat_class == "DNA"
        or repeat_class.startswith("DNA/")
    ):
        return "dna"

    return None


def parse_repeatmasker_out(
    repeat_out: Path,
    selected: set[str],
    repeat_to_genome: dict[str, str],
):

    intervals = {
        "gypsy": defaultdict(list),
        "copia": defaultdict(list),
        "dna": defaultdict(list),
    }

    class_counts = defaultdict(int)

    with open_text(repeat_out) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            fields = line.split()

            # RepeatMasker data rows normally begin with SW score.
            if len(fields) < 11:
                continue

            try:
                int(fields[0])
            except ValueError:
                continue

            repeat_chr = fields[4]

            genome_chr = repeat_to_genome.get(
                repeat_chr
            )

            if genome_chr is None:
                continue

            if genome_chr not in selected:
                continue

            try:
                start = int(fields[5]) - 1
                end = int(fields[6])

            except ValueError:
                continue

            repeat_class = fields[10]

            category = classify_repeat(
                repeat_class
            )

            if category is None:
                continue

            if end <= start:
                continue

            intervals[
                category
            ][
                genome_chr
            ].append(
                (
                    start,
                    end,
                )
            )

            class_counts[
                category
            ] += 1

    print()
    print("============================================================")
    print("RepeatMasker class parsing")
    print("============================================================")

    print(
        f"LTR/Gypsy records: "
        f"{class_counts['gypsy']:,}"
    )

    print(
        f"LTR/Copia records: "
        f"{class_counts['copia']:,}"
    )

    print(
        f"DNA-transposon records: "
        f"{class_counts['dna']:,}"
    )

    for category in [
        "gypsy",
        "copia",
        "dna",
    ]:

        if class_counts[
            category
        ] == 0:

            raise SystemExit(
                f"ERROR: Zero records found for {category}."
            )

    return intervals


###############################################################################
# Interval helpers
###############################################################################

def merge_intervals(intervals):

    if not intervals:
        return []

    intervals = sorted(
        intervals
    )

    merged = [
        [
            intervals[0][0],
            intervals[0][1],
        ]
    ]

    for start, end in intervals[1:]:

        previous = merged[-1]

        if start <= previous[1]:

            previous[1] = max(
                previous[1],
                end,
            )

        else:

            merged.append(
                [
                    start,
                    end,
                ]
            )

    return [
        (
            start,
            end,
        )
        for start, end
        in merged
    ]


def overlap_bp(
    intervals,
    window_start,
    window_end,
):

    total = 0

    for start, end in intervals:

        if end <= window_start:
            continue

        if start >= window_end:
            break

        overlap = (
            min(
                end,
                window_end,
            )
            - max(
                start,
                window_start,
            )
        )

        if overlap > 0:
            total += overlap

    return total


###############################################################################
# GC
###############################################################################

def gc_fraction(sequence):

    sequence = sequence.upper()

    a = sequence.count("A")
    c = sequence.count("C")
    g = sequence.count("G")
    t = sequence.count("T")

    denominator = (
        a + c + g + t
    )

    if denominator == 0:
        return float("nan")

    return (
        g + c
    ) / denominator


###############################################################################
# Window calculations
###############################################################################

def calculate_windows(
    chromosomes,
    lengths,
    sequences,
    genes,
    repeats,
    window_size,
):

    rows = []

    merged_repeats = {}

    for category in [
        "gypsy",
        "copia",
        "dna",
    ]:

        merged_repeats[
            category
        ] = {}

        for chromosome in chromosomes:

            merged_repeats[
                category
            ][
                chromosome
            ] = merge_intervals(
                repeats[
                    category
                ].get(
                    chromosome,
                    [],
                )
            )

    for chromosome in chromosomes:

        chromosome_length = lengths[
            chromosome
        ]

        chromosome_genes = sorted(
            genes.get(
                chromosome,
                [],
            )
        )

        sequence = sequences[
            chromosome
        ]

        n_windows = math.ceil(
            chromosome_length
            / window_size
        )

        for index in range(
            n_windows
        ):

            start = (
                index
                * window_size
            )

            end = min(
                start + window_size,
                chromosome_length,
            )

            window_bp = (
                end - start
            )

            gene_count = sum(
                start <= position < end
                for position
                in chromosome_genes
            )

            gypsy_bp = overlap_bp(
                merged_repeats[
                    "gypsy"
                ][
                    chromosome
                ],
                start,
                end,
            )

            copia_bp = overlap_bp(
                merged_repeats[
                    "copia"
                ][
                    chromosome
                ],
                start,
                end,
            )

            dna_bp = overlap_bp(
                merged_repeats[
                    "dna"
                ][
                    chromosome
                ],
                start,
                end,
            )

            gc = gc_fraction(
                sequence[
                    start:end
                ]
            )

            rows.append(
                {
                    "chromosome": chromosome,
                    "window_index": index + 1,
                    "start": start,
                    "end": end,
                    "window_bp": window_bp,
                    "gene_count": gene_count,

                    "gypsy_bp": gypsy_bp,
                    "gypsy_fraction": (
                        gypsy_bp / window_bp
                    ),

                    "copia_bp": copia_bp,
                    "copia_fraction": (
                        copia_bp / window_bp
                    ),

                    "dna_bp": dna_bp,
                    "dna_fraction": (
                        dna_bp / window_bp
                    ),

                    "gc_fraction": gc,
                }
            )

    return rows


###############################################################################
# Circular geometry
###############################################################################

def chromosome_layout(
    chromosomes,
    lengths,
    gap_degrees,
):

    gap = math.radians(
        gap_degrees
    )

    total_gap = (
        len(chromosomes)
        * gap
    )

    usable = (
        2 * math.pi
        - total_gap
    )

    total_bp = sum(
        lengths[c]
        for c
        in chromosomes
    )

    layout = {}

    current = 0.0

    for chromosome in chromosomes:

        width = (
            usable
            * lengths[
                chromosome
            ]
            / total_bp
        )

        start = current
        end = start + width

        layout[
            chromosome
        ] = {
            "start": start,
            "end": end,
            "center": (
                start + end
            ) / 2,
            "width": width,
        }

        current = (
            end + gap
        )

    return layout


def window_angle(
    row,
    lengths,
    layout,
):

    chromosome = row[
        "chromosome"
    ]

    chromosome_length = lengths[
        chromosome
    ]

    data = layout[
        chromosome
    ]

    frac_start = (
        row["start"]
        / chromosome_length
    )

    frac_end = (
        row["end"]
        / chromosome_length
    )

    theta_start = (
        data["start"]
        + data["width"]
        * frac_start
    )

    theta_end = (
        data["start"]
        + data["width"]
        * frac_end
    )

    return (
        (
            theta_start
            + theta_end
        ) / 2,
        (
            theta_end
            - theta_start
        ),
    )


def label_rotation(theta):

    degrees = math.degrees(
        theta
    )

    rotation = -degrees

    while rotation <= -180:
        rotation += 360

    while rotation > 180:
        rotation -= 360

    if rotation < -90:
        rotation += 180

    elif rotation > 90:
        rotation -= 180

    return rotation


###############################################################################
# Plot
###############################################################################

def make_plot(
    rows,
    chromosomes,
    lengths,
    display_names,
    output_prefix,
    gap_degrees,
    dpi,
):

    layout = chromosome_layout(
        chromosomes,
        lengths,
        gap_degrees,
    )

    gene_max = max(
        row[
            "gene_count"
        ]
        for row
        in rows
    )

    gc_values = [
        row[
            "gc_fraction"
        ]
        for row
        in rows
        if np.isfinite(
            row[
                "gc_fraction"
            ]
        )
    ]

    gc_min = min(
        gc_values
    )

    gc_max = max(
        gc_values
    )

    gypsy_max = max(
        row[
            "gypsy_fraction"
        ]
        for row
        in rows
    )

    copia_max = max(
        row[
            "copia_fraction"
        ]
        for row
        in rows
    )

    dna_max = max(
        row[
            "dna_fraction"
        ]
        for row
        in rows
    )

    print()
    print("============================================================")
    print("Track ranges")
    print("============================================================")

    print(
        f"Gene max:       {gene_max}/Mb"
    )

    print(
        f"Gypsy max:      {gypsy_max:.4f}"
    )

    print(
        f"Copia max:      {copia_max:.4f}"
    )

    print(
        f"DNA max:        {dna_max:.4f}"
    )

    print(
        f"GC range:       "
        f"{gc_min:.4f}-{gc_max:.4f}"
    )

    fig = plt.figure(
        figsize=(11, 11),
        facecolor="white",
    )

    ax = fig.add_subplot(
        111,
        projection="polar",
    )

    ax.set_theta_offset(
        math.pi / 2
    )

    ax.set_theta_direction(
        -1
    )

    ax.set_ylim(
        0,
        1.01,
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)

    ax.spines[
        "polar"
    ].set_visible(
        False
    )

    ###########################################################################
    # Chromosome bars
    ###########################################################################

    for chromosome in chromosomes:

        data = layout[
            chromosome
        ]

        ax.bar(
            data[
                "center"
            ],
            CHR_HEIGHT,
            width=data[
                "width"
            ],
            bottom=CHR_BOTTOM,
            color=CHR_COLOR,
            edgecolor=CHR_EDGE,
            linewidth=0.9,
            align="center",
            zorder=6,
        )

        ax.text(
            data[
                "center"
            ],
            CHR_BOTTOM
            + CHR_HEIGHT / 2,
            display_names[
                chromosome
            ],
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            rotation=label_rotation(
                data[
                    "center"
                ]
            ),
            rotation_mode="anchor",
            zorder=10,
        )

    ###########################################################################
    # Tracks
    ###########################################################################

    for row in rows:

        theta, width = window_angle(
            row,
            lengths,
            layout,
        )

        width *= 0.96

        # Gene count
        ax.bar(
            theta,
            GENE_HEIGHT
            * row[
                "gene_count"
            ]
            / gene_max,
            width=width,
            bottom=GENE_BOTTOM,
            color=GENE_COLOR,
            edgecolor="none",
            zorder=3,
        )

        # Gypsy
        ax.bar(
            theta,
            GYPSY_HEIGHT
            * row[
                "gypsy_fraction"
            ],
            width=width,
            bottom=GYPSY_BOTTOM,
            color=GYPSY_COLOR,
            edgecolor="none",
            zorder=3,
        )

        # Copia
        ax.bar(
            theta,
            COPIA_HEIGHT
            * row[
                "copia_fraction"
            ],
            width=width,
            bottom=COPIA_BOTTOM,
            color=COPIA_COLOR,
            edgecolor="none",
            zorder=3,
        )

        # DNA transposons
        ax.bar(
            theta,
            DNA_HEIGHT
            * row[
                "dna_fraction"
            ],
            width=width,
            bottom=DNA_BOTTOM,
            color=DNA_COLOR,
            edgecolor="none",
            zorder=3,
        )

        # GC
        if np.isfinite(
            row[
                "gc_fraction"
            ]
        ):

            normalized_gc = (
                row[
                    "gc_fraction"
                ]
                - gc_min
            ) / (
                gc_max
                - gc_min
            )

            ax.bar(
                theta,
                GC_HEIGHT
                * normalized_gc,
                width=width,
                bottom=GC_BOTTOM,
                color=GC_COLOR,
                edgecolor="none",
                zorder=3,
            )

    ###########################################################################
    # Baselines
    ###########################################################################

    theta_values = np.linspace(
        0,
        2 * math.pi,
        2000,
    )

    for radius in [
        GENE_BOTTOM,
        GYPSY_BOTTOM,
        COPIA_BOTTOM,
        DNA_BOTTOM,
        GC_BOTTOM,
    ]:

        ax.plot(
            theta_values,
            np.full(
                len(theta_values),
                radius,
            ),
            color="#B8B8B8",
            linewidth=0.4,
            zorder=1,
        )

    ###########################################################################
    # Centre
    ###########################################################################

    ax.text(
        0.5,
        0.515,
        "VPAN",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=17,
        fontweight="bold",
    )

    ax.text(
        0.5,
        0.482,
        "1-Mb windows",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=9,
    )

    ###########################################################################
    # Legend
    ###########################################################################

    legend = [
        mpatches.Patch(
            color=GENE_COLOR,
            label="Gene count",
        ),
        mpatches.Patch(
            color=GYPSY_COLOR,
            label="LTR/Gypsy",
        ),
        mpatches.Patch(
            color=COPIA_COLOR,
            label="LTR/Copia",
        ),
        mpatches.Patch(
            color=DNA_COLOR,
            label="DNA transposons",
        ),
        mpatches.Patch(
            color=GC_COLOR,
            label="GC content",
        ),
    ]

    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            0.008,
        ),
        ncol=5,
        frameon=False,
        fontsize=9,
    )

    fig.text(
        0.02,
        0.012,
        (
            f"Gene max: {gene_max}/Mb; "
            f"Gypsy max: {gypsy_max:.3f}; "
            f"Copia max: {copia_max:.3f}; "
            f"DNA max: {dna_max:.3f}; "
            f"GC: {gc_min:.3f}-{gc_max:.3f}"
        ),
        fontsize=6.5,
    )

    plt.subplots_adjust(
        left=0.025,
        right=0.975,
        top=0.975,
        bottom=0.075,
    )

    pdf = Path(
        str(
            output_prefix
        )
        + ".pdf"
    )

    svg = Path(
        str(
            output_prefix
        )
        + ".svg"
    )

    png = Path(
        str(
            output_prefix
        )
        + ".png"
    )

    fig.savefig(
        pdf,
        bbox_inches="tight",
    )

    fig.savefig(
        svg,
        bbox_inches="tight",
    )

    fig.savefig(
        png,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(fig)


###############################################################################
# Main
###############################################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genome",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--gff",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--repeat-out",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--chromosome-map",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--scaffolds",
        required=True,
    )

    parser.add_argument(
        "--window-size",
        type=int,
        default=1_000_000,
    )

    parser.add_argument(
        "--gap-degrees",
        type=float,
        default=2.5,
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
    )

    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
    )

    args = parser.parse_args()

    chromosomes = [
        item.strip()
        for item
        in args.scaffolds.split(",")
        if item.strip()
    ]

    selected = set(
        chromosomes
    )

    (
        repeat_to_genome,
        display_names,
    ) = read_chromosome_map(
        args.chromosome_map
    )

    lengths = fasta_lengths(
        args.genome
    )

    sequences = load_selected_sequences(
        args.genome,
        selected,
    )

    genes = parse_gene_midpoints(
        args.gff,
        selected,
    )

    repeats = parse_repeatmasker_out(
        args.repeat_out,
        selected,
        repeat_to_genome,
    )

    rows = calculate_windows(
        chromosomes,
        lengths,
        sequences,
        genes,
        repeats,
        args.window_size,
    )

    ###########################################################################
    # TSV
    ###########################################################################

    stats_file = Path(
        str(
            args.output_prefix
        )
        + ".window_statistics.tsv"
    )

    stats_file.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    columns = [
        "chromosome",
        "display_chromosome",
        "window_index",
        "start",
        "end",
        "window_bp",
        "gene_count",
        "gypsy_bp",
        "gypsy_fraction",
        "copia_bp",
        "copia_fraction",
        "dna_bp",
        "dna_fraction",
        "gc_fraction",
    ]

    with stats_file.open(
        "w",
        encoding="utf-8",
    ) as handle:

        handle.write(
            "\t".join(
                columns
            )
            + "\n"
        )

        for row in rows:

            values = [
                row["chromosome"],
                display_names[
                    row[
                        "chromosome"
                    ]
                ],
                str(
                    row[
                        "window_index"
                    ]
                ),
                str(
                    row[
                        "start"
                    ]
                ),
                str(
                    row[
                        "end"
                    ]
                ),
                str(
                    row[
                        "window_bp"
                    ]
                ),
                str(
                    row[
                        "gene_count"
                    ]
                ),
                str(
                    row[
                        "gypsy_bp"
                    ]
                ),
                f"{row['gypsy_fraction']:.8f}",
                str(
                    row[
                        "copia_bp"
                    ]
                ),
                f"{row['copia_fraction']:.8f}",
                str(
                    row[
                        "dna_bp"
                    ]
                ),
                f"{row['dna_fraction']:.8f}",
                f"{row['gc_fraction']:.8f}",
            ]

            handle.write(
                "\t".join(
                    values
                )
                + "\n"
            )

    make_plot(
        rows,
        chromosomes,
        lengths,
        display_names,
        args.output_prefix,
        args.gap_degrees,
        args.dpi,
    )

    print()
    print("Completed successfully.")
    print(
        f"Statistics: {stats_file}"
    )
    print(
        f"PDF: {args.output_prefix}.pdf"
    )
    print(
        f"SVG: {args.output_prefix}.svg"
    )
    print(
        f"PNG: {args.output_prefix}.png"
    )


if __name__ == "__main__":
    main()
