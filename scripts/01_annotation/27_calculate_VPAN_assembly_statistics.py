#!/usr/bin/env python3

"""
Calculate assembly statistics from a scaffold/chromosome FASTA.

Reports:
    assembled genome size
    number of scaffolds
    scaffold N50
    longest scaffold
    number of pseudochromosomes
    total N bases
    assembly gap percentage
    number of N-runs

Additionally reconstructs gap-delimited contigs by splitting scaffold
sequences at runs of N/n characters.

Important:
The reconstructed contig statistics represent sequence blocks separated by
assembly gaps. They should not automatically be interpreted as the original
pre-scaffolding contigs unless the assembly workflow preserved those
boundaries exactly.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path


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


def read_fasta(path: Path):

    sequences = {}

    name = None
    chunks = []

    with open_text(path) as handle:

        for line in handle:

            if line.startswith(">"):

                if name is not None:
                    sequences[name] = "".join(chunks)

                name = (
                    line[1:]
                    .strip()
                    .split()[0]
                )

                chunks = []

            else:

                chunks.append(
                    line.strip()
                )

    if name is not None:
        sequences[name] = "".join(chunks)

    return sequences


def n50(lengths):

    if not lengths:
        return 0

    total = sum(lengths)
    half = total / 2

    cumulative = 0

    for length in sorted(
        lengths,
        reverse=True,
    ):

        cumulative += length

        if cumulative >= half:
            return length

    return 0


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--min-pseudochromosome-bp",
        type=int,
        default=1_000_000,
        help=(
            "Minimum sequence length counted as a pseudochromosome. "
            "Default: 1 Mb."
        ),
    )

    args = parser.parse_args()

    sequences = read_fasta(
        args.fasta
    )

    if not sequences:

        raise SystemExit(
            "ERROR: no FASTA sequences found."
        )

    scaffold_lengths = [
        len(seq)
        for seq
        in sequences.values()
    ]

    total_bp = sum(
        scaffold_lengths
    )

    # ----------------------------------------------------------
    # N/gap statistics
    # ----------------------------------------------------------

    total_n_bases = 0
    n_runs = []

    inferred_contig_lengths = []

    for seq in sequences.values():

        runs = re.findall(
            r"[Nn]+",
            seq,
        )

        total_n_bases += sum(
            len(run)
            for run
            in runs
        )

        n_runs.extend(
            len(run)
            for run
            in runs
        )

        # Split into non-N sequence blocks.
        pieces = re.split(
            r"[Nn]+",
            seq,
        )

        inferred_contig_lengths.extend(
            len(piece)
            for piece
            in pieces
            if len(piece) > 0
        )

    gap_percent = (
        100
        * total_n_bases
        / total_bp
    )

    pseudochromosomes = [
        name
        for name, seq
        in sequences.items()
        if len(seq)
        >= args.min_pseudochromosome_bp
    ]

    # ----------------------------------------------------------
    # Output
    # ----------------------------------------------------------

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    with args.output.open(
        "w",
        encoding="utf-8",
    ) as out:

        out.write(
            "metric\tvalue\n"
        )

        out.write(
            f"assembled_genome_size_bp\t"
            f"{total_bp}\n"
        )

        out.write(
            f"assembled_genome_size_Mb\t"
            f"{total_bp / 1_000_000:.3f}\n"
        )

        out.write(
            f"number_scaffolds\t"
            f"{len(scaffold_lengths)}\n"
        )

        out.write(
            f"scaffold_N50_bp\t"
            f"{n50(scaffold_lengths)}\n"
        )

        out.write(
            f"scaffold_N50_Mb\t"
            f"{n50(scaffold_lengths) / 1_000_000:.3f}\n"
        )

        out.write(
            f"longest_scaffold_bp\t"
            f"{max(scaffold_lengths)}\n"
        )

        out.write(
            f"longest_scaffold_Mb\t"
            f"{max(scaffold_lengths) / 1_000_000:.3f}\n"
        )

        out.write(
            f"number_pseudochromosomes\t"
            f"{len(pseudochromosomes)}\n"
        )

        out.write(
            f"total_N_bases\t"
            f"{total_n_bases}\n"
        )

        out.write(
            f"assembly_gap_percent\t"
            f"{gap_percent:.6f}\n"
        )

        out.write(
            f"number_N_runs\t"
            f"{len(n_runs)}\n"
        )

        out.write(
            f"inferred_contig_count\t"
            f"{len(inferred_contig_lengths)}\n"
        )

        out.write(
            f"inferred_contig_N50_bp\t"
            f"{n50(inferred_contig_lengths)}\n"
        )

        out.write(
            f"inferred_contig_N50_Mb\t"
            f"{n50(inferred_contig_lengths) / 1_000_000:.3f}\n"
        )

    # Chromosome/scaffold lengths
    length_table = Path(
        str(args.output)
        + ".scaffold_lengths.tsv"
    )

    with length_table.open(
        "w",
        encoding="utf-8",
    ) as out:

        out.write(
            "sequence\tlength_bp\tlength_Mb\n"
        )

        for name, seq in sequences.items():

            out.write(
                f"{name}\t"
                f"{len(seq)}\t"
                f"{len(seq) / 1_000_000:.3f}\n"
            )

    print("Completed successfully.")
    print()
    print(f"Statistics: {args.output}")
    print(f"Lengths:    {length_table}")

    print()
    print("Assembly size:")
    print(
        f"  {total_bp:,} bp "
        f"({total_bp / 1_000_000:.3f} Mb)"
    )

    print(
        f"Scaffolds: {len(scaffold_lengths)}"
    )

    print(
        "Scaffold N50: "
        f"{n50(scaffold_lengths):,} bp"
    )

    print(
        f"N-runs: {len(n_runs):,}"
    )

    print(
        f"Inferred contigs: "
        f"{len(inferred_contig_lengths):,}"
    )

    print(
        "Inferred contig N50: "
        f"{n50(inferred_contig_lengths):,} bp"
    )


if __name__ == "__main__":
    main()
