#!/usr/bin/env python3

"""
Replot VPAN circular genome overview from precomputed window statistics.

Tracks from outside to inside:
    1. Gene count
    2. LTR/Gypsy
    3. LTR/Copia
    4. DNA transposons
    5. GC content

Chromosome ideograms:
    - labels: chr1, chr2, ...
    - labels centered inside each chromosome bar
    - 1-Mb minor ticks by default
    - 5-Mb major ticks by default
    - optional major tick labels

The TE-class tracks are scaled independently to their own observed maxima
for visual clarity. The underlying values remain unchanged in the input TSV.

Input TSV columns expected:
    chromosome
    display_chromosome
    window_index
    start
    end
    window_bp
    gene_count
    gypsy_bp
    gypsy_fraction
    copia_bp
    copia_fraction
    dna_bp
    dna_fraction
    gc_fraction

Outputs:
    <prefix>.pdf
    <prefix>.svg
    <prefix>.png

Example:
    python 05_scripts/replot_VPAN_circular_TEclasses_with_ticks.py \
        --stats figures/VPAN.circular_genome_TEclasses.window_statistics.tsv \
        --output-prefix figures/VPAN.circular_genome_TEclasses_final \
        --major-tick-mb 5 \
        --minor-tick-mb 1 \
        --show-tick-labels
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


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

GUIDE_COLOR = "#B8B8B8"
TICK_COLOR = "#222222"


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
# Required columns
###############################################################################

REQUIRED_COLUMNS = [
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


###############################################################################
# Geometry helpers
###############################################################################

def label_rotation(theta: float) -> float:
    """
    Return a readable tangential rotation for circular labels.

    The polar axis is configured with:
        0 degrees at 12 o'clock
        clockwise direction
    """

    degrees = math.degrees(theta)
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


def chromosome_layout(
    df: pd.DataFrame,
    chromosomes: list[str],
    gap_degrees: float = 2.5,
):
    """
    Allocate angular space to chromosomes proportional to their lengths.
    """

    lengths = (
        df.groupby("chromosome")["end"]
        .max()
        .to_dict()
    )

    gap = math.radians(gap_degrees)

    total_gap = (
        len(chromosomes)
        * gap
    )

    usable_angle = (
        2 * math.pi
        - total_gap
    )

    total_bp = sum(
        lengths[chrom]
        for chrom
        in chromosomes
    )

    layout = {}

    current = 0.0

    for chrom in chromosomes:

        angular_width = (
            usable_angle
            * lengths[chrom]
            / total_bp
        )

        start = current
        end = start + angular_width

        layout[chrom] = {
            "start": start,
            "end": end,
            "center": (
                start + end
            ) / 2,
            "width": angular_width,
            "length": lengths[chrom],
        }

        current = (
            end
            + gap
        )

    return layout


def window_angle(
    row,
    layout,
):
    """
    Convert a genomic window to angular center and angular width.
    """

    chrom = row["chromosome"]
    data = layout[chrom]

    start_fraction = (
        row["start"]
        / data["length"]
    )

    end_fraction = (
        row["end"]
        / data["length"]
    )

    theta_start = (
        data["start"]
        + data["width"]
        * start_fraction
    )

    theta_end = (
        data["start"]
        + data["width"]
        * end_fraction
    )

    theta_center = (
        theta_start
        + theta_end
    ) / 2

    theta_width = (
        theta_end
        - theta_start
    )

    return (
        theta_center,
        theta_width,
    )


###############################################################################
# Chromosome scale graduations
###############################################################################

def draw_chromosome_ticks(
    ax,
    chromosomes,
    layout,
    chr_bottom,
    chr_height,
    major_tick_mb=5.0,
    minor_tick_mb=1.0,
    show_tick_labels=False,
):
    """
    Draw chromosome-length graduations along each chromosome ideogram.

    Minor ticks:
        default every 1 Mb

    Major ticks:
        default every 5 Mb

    Major tick labels:
        optional, controlled by --show-tick-labels

    Labels are shown as:
        5, 10, 15, ...
    representing Mb from the start of each chromosome.
    """

    if major_tick_mb <= 0:
        raise ValueError(
            "major_tick_mb must be greater than zero."
        )

    if minor_tick_mb <= 0:
        raise ValueError(
            "minor_tick_mb must be greater than zero."
        )

    major_bp = int(
        major_tick_mb
        * 1_000_000
    )

    minor_bp = int(
        minor_tick_mb
        * 1_000_000
    )

    minor_tick_length = 0.010
    major_tick_length = 0.018

    tick_label_offset = 0.017

    minor_linewidth = 0.45
    major_linewidth = 0.80

    label_fontsize = 6.5

    for chrom in chromosomes:

        data = layout[chrom]
        chrom_length = int(
            data["length"]
        )

        # Tick positions start after 0 Mb.
        for bp in range(
            minor_bp,
            chrom_length,
            minor_bp,
        ):

            fraction = (
                bp
                / chrom_length
            )

            theta = (
                data["start"]
                + data["width"]
                * fraction
            )

            is_major = (
                bp % major_bp == 0
            )

            if is_major:

                tick_length = (
                    major_tick_length
                )

                linewidth = (
                    major_linewidth
                )

            else:

                tick_length = (
                    minor_tick_length
                )

                linewidth = (
                    minor_linewidth
                )

            # Draw ticks outward from chromosome bar.
            radius_start = (
                chr_bottom
                + chr_height
            )

            radius_end = (
                radius_start
                + tick_length
            )

            ax.plot(
                [
                    theta,
                    theta,
                ],
                [
                    radius_start,
                    radius_end,
                ],
                color=TICK_COLOR,
                linewidth=linewidth,
                solid_capstyle="butt",
                zorder=12,
            )

            if (
                is_major
                and show_tick_labels
            ):

                tick_label = str(
                    int(
                        bp
                        / 1_000_000
                    )
                )

                ax.text(
                    theta,
                    radius_end
                    + tick_label_offset,
                    tick_label,
                    ha="center",
                    va="center",
                    fontsize=label_fontsize,
                    rotation=label_rotation(
                        theta
                    ),
                    rotation_mode="anchor",
                    color="black",
                    zorder=13,
                )


###############################################################################
# Plot
###############################################################################

def make_plot(
    df: pd.DataFrame,
    chromosomes: list[str],
    display_names: dict[str, str],
    output_prefix: Path,
    gap_degrees: float,
    major_tick_mb: float,
    minor_tick_mb: float,
    show_tick_labels: bool,
    dpi: int,
):
    """
    Generate final circular figure.
    """

    layout = chromosome_layout(
        df,
        chromosomes,
        gap_degrees=gap_degrees,
    )

    ###########################################################################
    # Track scaling
    ###########################################################################

    gene_max = float(
        df["gene_count"].max()
    )

    gypsy_max = float(
        df["gypsy_fraction"].max()
    )

    copia_max = float(
        df["copia_fraction"].max()
    )

    dna_max = float(
        df["dna_fraction"].max()
    )

    gc_min = float(
        df["gc_fraction"].min()
    )

    gc_max = float(
        df["gc_fraction"].max()
    )

    if gene_max <= 0:
        raise SystemExit(
            "ERROR: gene_count maximum is zero."
        )

    if gypsy_max <= 0:
        raise SystemExit(
            "ERROR: gypsy_fraction maximum is zero."
        )

    if copia_max <= 0:
        raise SystemExit(
            "ERROR: copia_fraction maximum is zero."
        )

    if dna_max <= 0:
        raise SystemExit(
            "ERROR: dna_fraction maximum is zero."
        )

    if gc_max <= gc_min:
        raise SystemExit(
            "ERROR: invalid GC-content range."
        )

    print()
    print("============================================================")
    print("Track ranges used for plotting")
    print("============================================================")

    print(
        f"Gene count maximum: {gene_max:.0f}/Mb"
    )

    print(
        f"LTR/Gypsy maximum fraction: {gypsy_max:.6f}"
    )

    print(
        f"LTR/Copia maximum fraction: {copia_max:.6f}"
    )

    print(
        f"DNA-transposon maximum fraction: {dna_max:.6f}"
    )

    print(
        f"GC-content range: {gc_min:.6f}-{gc_max:.6f}"
    )

    ###########################################################################
    # Figure
    ###########################################################################

    fig = plt.figure(
        figsize=(11, 11),
        facecolor="white",
    )

    ax = fig.add_subplot(
        111,
        projection="polar",
    )

    # 12 o'clock start.
    ax.set_theta_offset(
        math.pi / 2
    )

    # Clockwise.
    ax.set_theta_direction(
        -1
    )

    # Extended slightly above chromosome bar for ticks/labels.
    if show_tick_labels:

        radial_max = 1.08

    else:

        radial_max = 1.03

    ax.set_ylim(
        0,
        radial_max,
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
    # Chromosome ideograms
    ###########################################################################

    for chrom in chromosomes:

        data = layout[
            chrom
        ]

        ax.bar(
            data["center"],
            CHR_HEIGHT,
            width=data["width"],
            bottom=CHR_BOTTOM,
            color=CHR_COLOR,
            edgecolor=CHR_EDGE,
            linewidth=0.9,
            align="center",
            zorder=6,
        )

        ax.text(
            data["center"],
            CHR_BOTTOM
            + CHR_HEIGHT / 2,
            display_names[
                chrom
            ],
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            rotation=label_rotation(
                data["center"]
            ),
            rotation_mode="anchor",
            color="black",
            zorder=10,
        )

    ###########################################################################
    # Chromosome graduations
    ###########################################################################

    draw_chromosome_ticks(
        ax=ax,
        chromosomes=chromosomes,
        layout=layout,
        chr_bottom=CHR_BOTTOM,
        chr_height=CHR_HEIGHT,
        major_tick_mb=major_tick_mb,
        minor_tick_mb=minor_tick_mb,
        show_tick_labels=show_tick_labels,
    )

    ###########################################################################
    # Tracks
    ###########################################################################

    for _, row in df.iterrows():

        theta, width = window_angle(
            row,
            layout,
        )

        # Small gap between adjacent genomic windows.
        width *= 0.96

        # ------------------------------------------------------
        # Gene count
        # ------------------------------------------------------

        gene_height = (
            GENE_HEIGHT
            * row["gene_count"]
            / gene_max
        )

        ax.bar(
            theta,
            gene_height,
            width=width,
            bottom=GENE_BOTTOM,
            color=GENE_COLOR,
            edgecolor="none",
            align="center",
            zorder=3,
        )

        # ------------------------------------------------------
        # LTR/Gypsy
        # independently normalized
        # ------------------------------------------------------

        gypsy_height = (
            GYPSY_HEIGHT
            * row["gypsy_fraction"]
            / gypsy_max
        )

        ax.bar(
            theta,
            gypsy_height,
            width=width,
            bottom=GYPSY_BOTTOM,
            color=GYPSY_COLOR,
            edgecolor="none",
            align="center",
            zorder=3,
        )

        # ------------------------------------------------------
        # LTR/Copia
        # independently normalized
        # ------------------------------------------------------

        copia_height = (
            COPIA_HEIGHT
            * row["copia_fraction"]
            / copia_max
        )

        ax.bar(
            theta,
            copia_height,
            width=width,
            bottom=COPIA_BOTTOM,
            color=COPIA_COLOR,
            edgecolor="none",
            align="center",
            zorder=3,
        )

        # ------------------------------------------------------
        # DNA transposons
        # independently normalized
        # ------------------------------------------------------

        dna_height = (
            DNA_HEIGHT
            * row["dna_fraction"]
            / dna_max
        )

        ax.bar(
            theta,
            dna_height,
            width=width,
            bottom=DNA_BOTTOM,
            color=DNA_COLOR,
            edgecolor="none",
            align="center",
            zorder=3,
        )

        # ------------------------------------------------------
        # GC content
        # ------------------------------------------------------

        gc_value = row[
            "gc_fraction"
        ]

        if np.isfinite(
            gc_value
        ):

            gc_normalized = (
                gc_value
                - gc_min
            ) / (
                gc_max
                - gc_min
            )

            gc_height = (
                GC_HEIGHT
                * gc_normalized
            )

            ax.bar(
                theta,
                gc_height,
                width=width,
                bottom=GC_BOTTOM,
                color=GC_COLOR,
                edgecolor="none",
                align="center",
                zorder=3,
            )

    ###########################################################################
    # Track guide circles
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
            color=GUIDE_COLOR,
            linewidth=0.4,
            zorder=1,
        )

    ###########################################################################
    # Centre labels
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

    legend_handles = [
        mpatches.Patch(
            facecolor=GENE_COLOR,
            edgecolor="none",
            label="Gene count",
        ),
        mpatches.Patch(
            facecolor=GYPSY_COLOR,
            edgecolor="none",
            label="LTR/Gypsy",
        ),
        mpatches.Patch(
            facecolor=COPIA_COLOR,
            edgecolor="none",
            label="LTR/Copia",
        ),
        mpatches.Patch(
            facecolor=DNA_COLOR,
            edgecolor="none",
            label="DNA transposons",
        ),
        mpatches.Patch(
            facecolor=GC_COLOR,
            edgecolor="none",
            label="GC content",
        ),
    ]

    fig.legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(
            0.5,
            0.018,
        ),
        ncol=5,
        frameon=False,
        fontsize=9,
    )

    ###########################################################################
    # Layout
    ###########################################################################

    if show_tick_labels:

        top_margin = 0.96

    else:

        top_margin = 0.975

    plt.subplots_adjust(
        left=0.025,
        right=0.975,
        top=top_margin,
        bottom=0.075,
    )

    ###########################################################################
    # Save
    ###########################################################################

    output_prefix.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    pdf_path = Path(
        str(
            output_prefix
        )
        + ".pdf"
    )

    svg_path = Path(
        str(
            output_prefix
        )
        + ".svg"
    )

    png_path = Path(
        str(
            output_prefix
        )
        + ".png"
    )

    fig.savefig(
        pdf_path,
        bbox_inches="tight",
    )

    fig.savefig(
        svg_path,
        bbox_inches="tight",
    )

    fig.savefig(
        png_path,
        dpi=dpi,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )

    print()
    print("============================================================")
    print("Created figure files")
    print("============================================================")

    print(
        pdf_path
    )

    print(
        svg_path
    )

    print(
        png_path
    )


###############################################################################
# Main
###############################################################################

def main():

    parser = argparse.ArgumentParser(
        description=(
            "Replot VPAN circular TE-class genome overview "
            "with chromosome-length graduations."
        )
    )

    parser.add_argument(
        "--stats",
        required=True,
        type=Path,
        help=(
            "Window-statistics TSV generated by the "
            "VPAN TE-class analysis."
        ),
    )

    parser.add_argument(
        "--output-prefix",
        required=True,
        type=Path,
        help=(
            "Output filename prefix."
        ),
    )

    parser.add_argument(
        "--gap-degrees",
        type=float,
        default=2.5,
        help=(
            "Angular gap between chromosomes."
        ),
    )

    parser.add_argument(
        "--major-tick-mb",
        type=float,
        default=5.0,
        help=(
            "Major chromosome tick interval in Mb. "
            "Default: 5."
        ),
    )

    parser.add_argument(
        "--minor-tick-mb",
        type=float,
        default=1.0,
        help=(
            "Minor chromosome tick interval in Mb. "
            "Default: 1."
        ),
    )

    parser.add_argument(
        "--show-tick-labels",
        action="store_true",
        help=(
            "Show numeric labels on major chromosome ticks."
        ),
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help=(
            "PNG resolution. Default: 600 dpi."
        ),
    )

    args = parser.parse_args()

    ###########################################################################
    # Validate input
    ###########################################################################

    if (
        not args.stats.is_file()
        or args.stats.stat().st_size == 0
    ):

        raise SystemExit(
            f"ERROR: Missing or empty statistics TSV:\n{args.stats}"
        )

    ###########################################################################
    # Read statistics
    ###########################################################################

    df = pd.read_csv(
        args.stats,
        sep="\t",
    )

    missing_columns = [
        column
        for column
        in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing_columns:

        raise SystemExit(
            "ERROR: Required columns missing from statistics TSV:\n"
            + "\n".join(
                missing_columns
            )
        )

    ###########################################################################
    # Preserve chromosome order from TSV
    ###########################################################################

    chromosome_table = (
        df[
            [
                "chromosome",
                "display_chromosome",
            ]
        ]
        .drop_duplicates()
    )

    chromosomes = (
        chromosome_table[
            "chromosome"
        ]
        .tolist()
    )

    display_names = (
        chromosome_table
        .set_index(
            "chromosome"
        )[
            "display_chromosome"
        ]
        .to_dict()
    )

    print()
    print("============================================================")
    print("Chromosomes")
    print("============================================================")

    for chrom in chromosomes:

        chromosome_length = int(
            df.loc[
                df[
                    "chromosome"
                ]
                == chrom,
                "end",
            ].max()
        )

        print(
            f"{display_names[chrom]}\t"
            f"{chromosome_length:,} bp\t"
            f"{chromosome_length / 1_000_000:.2f} Mb"
        )

    print()
    print(
        f"Minor chromosome ticks: "
        f"{args.minor_tick_mb:g} Mb"
    )

    print(
        f"Major chromosome ticks: "
        f"{args.major_tick_mb:g} Mb"
    )

    print(
        f"Major tick labels: "
        f"{'yes' if args.show_tick_labels else 'no'}"
    )

    ###########################################################################
    # Plot
    ###########################################################################

    make_plot(
        df=df,
        chromosomes=chromosomes,
        display_names=display_names,
        output_prefix=args.output_prefix,
        gap_degrees=args.gap_degrees,
        major_tick_mb=args.major_tick_mb,
        minor_tick_mb=args.minor_tick_mb,
        show_tick_labels=args.show_tick_labels,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
