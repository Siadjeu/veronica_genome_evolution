#!/usr/bin/env python3

from pathlib import Path
import csv
import sys


TOTAL = 28313

DIAMOND = Path(
    "results/annotation/functional/VPAN_swissprot.diamond.tsv"
)

INTERPRO = Path(
    "results/annotation/functional/interproscan/"
    "VPAN.interproscan.selected_apps.tsv"
)

OUTDIR = Path(
    "results/annotation/functional"
)

SUMMARY = OUTDIR / "VPAN_functional_annotation_summary.tsv"

EVIDENCE = OUTDIR / "VPAN_functional_annotation_evidence.tsv"


def pct(n):
    return 100.0 * n / TOTAL


def main():

    OUTDIR.mkdir(parents=True, exist_ok=True)

    if not DIAMOND.is_file():
        sys.exit(
            f"ERROR: DIAMOND result not found:\n{DIAMOND}"
        )

    if not INTERPRO.is_file():
        sys.exit(
            f"ERROR: InterProScan result not found:\n{INTERPRO}"
        )

    # ---------------------------------------------------------
    # Swiss-Prot strict support
    #
    # Expected DIAMOND columns:
    # 1  qseqid
    # 2  sseqid
    # 3  pident
    # 4  length
    # 5  qlen
    # 6  slen
    # 7  evalue
    # 8  bitscore
    # 9  qcovhsp
    # 10 scovhsp
    # 11 stitle
    #
    # Strict support:
    # pident >= 30%
    # qcovhsp >= 50%
    # ---------------------------------------------------------

    swiss = set()

    with DIAMOND.open() as fh:
        for lineno, line in enumerate(fh, start=1):

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 9:
                sys.exit(
                    "ERROR: malformed DIAMOND row "
                    f"at line {lineno}"
                )

            qseqid = fields[0]
            pident = float(fields[2])
            qcov = float(fields[8])

            if pident >= 30.0 and qcov >= 50.0:
                swiss.add(qseqid)

    # ---------------------------------------------------------
    # InterProScan support
    #
    # Any InterProScan row counts as domain/family support.
    # Protein accession is column 1.
    # ---------------------------------------------------------

    interpro = set()

    application_counts = {}

    with INTERPRO.open() as fh:
        for lineno, line in enumerate(fh, start=1):

            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) < 4:
                sys.exit(
                    "ERROR: malformed InterProScan row "
                    f"at line {lineno}"
                )

            protein_id = fields[0]
            application = fields[3]

            interpro.add(protein_id)

            application_counts[application] = (
                application_counts.get(application, 0) + 1
            )

    # ---------------------------------------------------------
    # Union and overlap
    # ---------------------------------------------------------

    both = swiss & interpro
    swiss_only = swiss - interpro
    interpro_only = interpro - swiss
    union = swiss | interpro
    unsupported = TOTAL - len(union)

    # ---------------------------------------------------------
    # Validation
    # ---------------------------------------------------------

    if len(union) > TOTAL:
        sys.exit(
            "ERROR: union contains more proteins than expected "
            f"({len(union)} > {TOTAL})"
        )

    # ---------------------------------------------------------
    # Terminal summary
    # ---------------------------------------------------------

    print()
    print("====================================================")
    print("VPAN FINAL FUNCTIONAL ANNOTATION")
    print("====================================================")
    print(f"Representative proteins : {TOTAL}")
    print()

    print(
        f"Swiss-Prot strict       : "
        f"{len(swiss):6d} "
        f"({pct(len(swiss)):.2f}%)"
    )

    print(
        f"InterProScan            : "
        f"{len(interpro):6d} "
        f"({pct(len(interpro)):.2f}%)"
    )

    print(
        f"Both                    : "
        f"{len(both):6d} "
        f"({pct(len(both)):.2f}%)"
    )

    print(
        f"Swiss-Prot only         : "
        f"{len(swiss_only):6d} "
        f"({pct(len(swiss_only)):.2f}%)"
    )

    print(
        f"InterProScan only       : "
        f"{len(interpro_only):6d} "
        f"({pct(len(interpro_only)):.2f}%)"
    )

    print("----------------------------------------------------")

    print(
        f"FUNCTIONALLY ANNOTATED  : "
        f"{len(union):6d} "
        f"({pct(len(union)):.2f}%)"
    )

    print(
        f"Unsupported             : "
        f"{unsupported:6d} "
        f"({pct(unsupported):.2f}%)"
    )

    print()
    print("InterProScan rows by application:")

    for application, count in sorted(
        application_counts.items(),
        key=lambda x: (-x[1], x[0])
    ):
        print(
            f"{application:15s} {count:8d}"
        )

    # ---------------------------------------------------------
    # Write summary table
    # ---------------------------------------------------------

    with SUMMARY.open("w", newline="") as fh:

        writer = csv.writer(
            fh,
            delimiter="\t"
        )

        writer.writerow([
            "category",
            "protein_count",
            "percentage"
        ])

        rows = [
            (
                "representative_proteins",
                TOTAL,
                100.0
            ),
            (
                "swissprot_strict",
                len(swiss),
                pct(len(swiss))
            ),
            (
                "interproscan",
                len(interpro),
                pct(len(interpro))
            ),
            (
                "both",
                len(both),
                pct(len(both))
            ),
            (
                "swissprot_only",
                len(swiss_only),
                pct(len(swiss_only))
            ),
            (
                "interproscan_only",
                len(interpro_only),
                pct(len(interpro_only))
            ),
            (
                "functional_annotation_union",
                len(union),
                pct(len(union))
            ),
            (
                "no_functional_support",
                unsupported,
                pct(unsupported)
            ),
        ]

        for category, count, percentage in rows:

            writer.writerow([
                category,
                count,
                f"{percentage:.2f}"
            ])

    # ---------------------------------------------------------
    # Write per-protein evidence table
    #
    # We include all proteins observed in either evidence source.
    # ---------------------------------------------------------

    all_supported = sorted(union)

    with EVIDENCE.open("w", newline="") as fh:

        writer = csv.writer(
            fh,
            delimiter="\t"
        )

        writer.writerow([
            "protein_id",
            "swissprot_strict",
            "interproscan",
            "functional_annotation_supported"
        ])

        for protein_id in all_supported:

            writer.writerow([
                protein_id,
                "YES" if protein_id in swiss else "NO",
                "YES" if protein_id in interpro else "NO",
                "YES"
            ])

    print()
    print("Saved:")
    print(SUMMARY)
    print(EVIDENCE)


if __name__ == "__main__":
    main()
