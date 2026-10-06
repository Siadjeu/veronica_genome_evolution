#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=6000
#SBATCH --job-name=macro_all9
#SBATCH --output=10_synteny/logs/macrosynteny_all9_%j.out
#SBATCH --error=10_synteny/logs/macrosynteny_all9_%j.err

set -euo pipefail

# ============================================================
# Project paths
# ============================================================

PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"
SYNTENY_DIR="${PROJECT_DIR}/10_synteny"

STEP34_CHECKPOINT="${SYNTENY_DIR}/checkpoint_jcvi_step34/JCVI_STEP34_COMPLETE.txt"
STEP34_RESULTS="${SYNTENY_DIR}/jcvi_step34/tables/jcvi_step34_results.tsv"
JCVI_MANIFEST="${SYNTENY_DIR}/manifests/jcvi_input_manifest.tsv"

OUTDIR="${SYNTENY_DIR}/jcvi_step34/macrosynteny_all_species"
BEDDIR="${OUTDIR}/numeric_beds"
ANCHORDIR="${OUTDIR}/anchors"
TABLEDIR="${OUTDIR}/tables"
CONFIGDIR="${OUTDIR}/config"
PLOTDIR="${OUTDIR}/plots"
LOGDIR="${OUTDIR}/logs"
CHECKPOINT_DIR="${SYNTENY_DIR}/checkpoint_macrosynteny_all_species"

mkdir -p \
    "${BEDDIR}" \
    "${ANCHORDIR}" \
    "${TABLEDIR}" \
    "${CONFIGDIR}" \
    "${PLOTDIR}" \
    "${LOGDIR}" \
    "${CHECKPOINT_DIR}" \
    "${SYNTENY_DIR}/logs"

cd "${PROJECT_DIR}"

# ============================================================
# Validate inputs
# ============================================================

for FILE in \
    "${STEP34_CHECKPOINT}" \
    "${STEP34_RESULTS}" \
    "${JCVI_MANIFEST}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Missing or empty required input:" >&2
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

export MPLBACKEND=Agg

echo "Python:"
python --version

echo "JCVI:"
python - <<'PY'
import jcvi
print(getattr(jcvi, "__version__", "unknown"))
PY

# ============================================================
# Clean old Step 34E outputs
# ============================================================

rm -f "${BEDDIR}"/*
rm -f "${ANCHORDIR}"/*
rm -f "${TABLEDIR}"/*
rm -f "${CONFIGDIR}"/*
rm -f "${PLOTDIR}"/*
rm -f "${LOGDIR}"/*
rm -f "${CHECKPOINT_DIR}"/*

# ============================================================
# Build numeric chromosome BED files and configurations
# ============================================================

python - \
    "${JCVI_MANIFEST}" \
    "${STEP34_RESULTS}" \
    "${BEDDIR}" \
    "${ANCHORDIR}" \
    "${TABLEDIR}" \
    "${CONFIGDIR}" <<'PY'
from __future__ import annotations

import csv
import gzip
import re
import shutil
import sys
from pathlib import Path

(
    manifest_name,
    results_name,
    bed_dir_name,
    anchor_dir_name,
    table_dir_name,
    config_dir_name,
) = sys.argv[1:]

manifest_file = Path(manifest_name)
results_file = Path(results_name)

bed_dir = Path(bed_dir_name)
anchor_dir = Path(anchor_dir_name)
table_dir = Path(table_dir_name)
config_dir = Path(config_dir_name)

for directory in [
    bed_dir,
    anchor_dir,
    table_dir,
    config_dir,
]:
    directory.mkdir(
        parents=True,
        exist_ok=True,
    )

# ============================================================
# Plotting order
#
# This order guarantees an available Step 34 anchor comparison
# between every adjacent pair.
# ============================================================

species_order = [
    "PMAJ",
    "VPAN",
    "VSCU",
    "VANA",
    "VARV",
    "VVER",
    "VPER",
    "VSER",
    "VTRI",
]

scientific_labels = {
    "PMAJ": "Plantago major",
    "VPAN": "Veronica panormitana",
    "VSCU": "Veronica scutellata",
    "VANA": "Veronica anagallis-aquatica",
    "VARV": "Veronica arvensis",
    "VVER": "Veronica verna",
    "VPER": "Veronica persica",
    "VSER": "Veronica serpyllifolia",
    "VTRI": "Veronica triphyllos",
}

# Adjacent track comparison and orientation.
#
# stored_query/stored_subject correspond to the gene columns
# in the Step 34 anchor file. Gene IDs are species-prefixed,
# so JCVI can connect them even when the visual order is the
# reverse of the original query-subject order.
edge_comparisons = [
    ("PMAJ", "VPAN", "PMAJ__VPAN"),
    ("VPAN", "VSCU", "VPAN__VSCU"),
    ("VSCU", "VANA", "VANA__VSCU"),
    ("VANA", "VARV", "VANA__VARV"),
    ("VARV", "VVER", "VARV__VVER"),
    ("VVER", "VPER", "VPER__VVER"),
    ("VPER", "VSER", "VPER__VSER"),
    ("VSER", "VTRI", "VSER__VTRI"),
]

with manifest_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    manifest_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

with results_file.open(
    newline="",
    encoding="utf-8-sig",
) as handle:
    result_rows = list(
        csv.DictReader(
            handle,
            delimiter="\t",
        )
    )

manifest = {
    row["species_code"]: row
    for row in manifest_rows
}

results = {
    row["comparison_id"]: row
    for row in result_rows
}

missing_species = [
    code
    for code in species_order
    if code not in manifest
]

if missing_species:
    raise SystemExit(
        "ERROR: Species missing from JCVI manifest: "
        + ",".join(missing_species)
    )


def open_text(path: Path):
    if path.name.lower().endswith(".gz"):
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


def fasta_headers(path: Path):
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def natural_key(value: str):
    parts = re.split(r"(\d+)", value)

    return [
        int(part)
        if part.isdigit()
        else part.lower()
        for part in parts
    ]


def chromosome_number_from_header(
    primary_id: str,
    full_header: str,
):
    patterns = [
        r"chromosome\s*[:=_-]?\s*([0-9]+)",
        r"\bchr(?:omosome)?[_-]?([0-9]+)\b",
        r"^Chr([0-9]+)$",
        r"^chr([0-9]+)$",
    ]

    searchable = (
        primary_id
        + " "
        + full_header
    )

    for pattern in patterns:
        match = re.search(
            pattern,
            searchable,
            flags=re.IGNORECASE,
        )

        if match:
            return int(match.group(1))

    return None


mapping_rows = []
species_seqids = {}
bed_qc_rows = []

for code in species_order:
    metadata = manifest[code]

    genome_file = Path(
        metadata["genome_fasta"]
    )
    bed_file = Path(
        metadata["bed_file"]
    )

    if not genome_file.is_file():
        raise SystemExit(
            f"ERROR: Missing genome for {code}: "
            f"{genome_file}"
        )

    if not bed_file.is_file():
        raise SystemExit(
            f"ERROR: Missing BED for {code}: "
            f"{bed_file}"
        )

    headers = list(
        fasta_headers(genome_file)
    )

    header_records = []

    for order_index, full_header in enumerate(
        headers,
        start=1,
    ):
        primary_id = full_header.split()[0]

        parsed_number = (
            chromosome_number_from_header(
                primary_id,
                full_header,
            )
        )

        header_records.append(
            {
                "primary_id": primary_id,
                "full_header": full_header,
                "fasta_order": order_index,
                "parsed_number": parsed_number,
            }
        )

    parsed_numbers = [
        record["parsed_number"]
        for record in header_records
    ]

    parsed_numbers_are_complete = (
        all(
            number is not None
            for number in parsed_numbers
        )
        and len(set(parsed_numbers))
        == len(parsed_numbers)
        and set(parsed_numbers)
        == set(
            range(
                1,
                len(parsed_numbers) + 1,
            )
        )
    )

    if parsed_numbers_are_complete:
        ordered_records = sorted(
            header_records,
            key=lambda record:
                record["parsed_number"],
        )

        mapping_method = (
            "PARSED_FROM_FASTA_HEADER"
        )
    else:
        # Preserve original FASTA pseudomolecule order.
        ordered_records = sorted(
            header_records,
            key=lambda record:
                record["fasta_order"],
        )

        mapping_method = (
            "FASTA_SEQUENCE_ORDER"
        )

    original_to_numeric = {}

    for numeric_number, record in enumerate(
        ordered_records,
        start=1,
    ):
        numeric_id = str(
            numeric_number
        )

        original_to_numeric[
            record["primary_id"]
        ] = numeric_id

        mapping_rows.append(
            {
                "species_code": code,
                "scientific_name": (
                    scientific_labels[code]
                ),
                "original_chromosome_id": (
                    record["primary_id"]
                ),
                "original_full_header": (
                    record["full_header"]
                ),
                "original_fasta_order": (
                    record["fasta_order"]
                ),
                "parsed_chromosome_number": (
                    record["parsed_number"]
                    if record[
                        "parsed_number"
                    ] is not None
                    else ""
                ),
                "numeric_plot_chromosome_id": (
                    numeric_id
                ),
                "mapping_method": (
                    mapping_method
                ),
            }
        )

    numeric_bed_file = (
        bed_dir / f"{code}.numeric.bed"
    )

    input_gene_count = 0
    output_gene_count = 0
    unknown_bed_chromosomes = set()

    with bed_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as source, numeric_bed_file.open(
        "w",
        encoding="utf-8",
    ) as destination:
        for line_number, line in enumerate(
            source,
            start=1,
        ):
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split(
                "\t"
            )

            if len(fields) < 6:
                raise SystemExit(
                    f"ERROR: Malformed BED line "
                    f"{line_number} in {bed_file}"
                )

            input_gene_count += 1

            chromosome = fields[0]

            numeric_chromosome = (
                original_to_numeric.get(
                    chromosome
                )
            )

            if numeric_chromosome is None:
                unknown_bed_chromosomes.add(
                    chromosome
                )
                continue

            fields[0] = numeric_chromosome

            destination.write(
                "\t".join(fields)
                + "\n"
            )

            output_gene_count += 1

    if unknown_bed_chromosomes:
        raise SystemExit(
            f"ERROR: BED chromosomes for {code} "
            f"were not found in genome mapping: "
            + ",".join(
                sorted(
                    unknown_bed_chromosomes,
                    key=natural_key,
                )
            )
        )

    if input_gene_count != output_gene_count:
        raise SystemExit(
            f"ERROR: BED gene count changed for "
            f"{code}: input={input_gene_count}; "
            f"output={output_gene_count}"
        )

    species_seqids[code] = [
        str(number)
        for number in range(
            1,
            len(ordered_records) + 1,
        )
    ]

    bed_qc_rows.append(
        {
            "species_code": code,
            "chromosome_count": len(
                ordered_records
            ),
            "input_bed_genes": (
                input_gene_count
            ),
            "output_bed_genes": (
                output_gene_count
            ),
            "mapping_method": (
                mapping_method
            ),
            "status": "PASS",
        }
    )

# ============================================================
# Copy the simple anchors required for adjacent tracks
# ============================================================

edge_rows = []

for upper_species, lower_species, comparison_id in (
    edge_comparisons
):
    if comparison_id not in results:
        raise SystemExit(
            f"ERROR: Step 34 comparison missing: "
            f"{comparison_id}"
        )

    result = results[comparison_id]

    if result["status"] != "PASS":
        raise SystemExit(
            f"ERROR: Step 34 comparison is not "
            f"PASS: {comparison_id}"
        )

    simple_anchor = Path(
        result["simple_anchor_file"]
    )

    if (
        not simple_anchor.is_file()
        or simple_anchor.stat().st_size == 0
    ):
        raise SystemExit(
            f"ERROR: Missing simple anchor file "
            f"for {comparison_id}: "
            f"{simple_anchor}"
        )

    copied_anchor = (
        anchor_dir
        / f"{comparison_id}.anchors.simple"
    )

    shutil.copyfile(
        simple_anchor,
        copied_anchor,
    )

    edge_rows.append(
        {
            "upper_species": upper_species,
            "lower_species": lower_species,
            "comparison_id": comparison_id,
            "source_anchor_file": str(
                simple_anchor
            ),
            "plot_anchor_file": str(
                copied_anchor
            ),
            "status": "PASS",
        }
    )

# ============================================================
# Write chromosome mapping table
# ============================================================

mapping_output = (
    table_dir
    / "numeric_chromosome_mapping.tsv"
)

with mapping_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "scientific_name",
            "original_chromosome_id",
            "original_full_header",
            "original_fasta_order",
            "parsed_chromosome_number",
            "numeric_plot_chromosome_id",
            "mapping_method",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(mapping_rows)

bed_qc_output = (
    table_dir
    / "numeric_bed_qc.tsv"
)

with bed_qc_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "species_code",
            "chromosome_count",
            "input_bed_genes",
            "output_bed_genes",
            "mapping_method",
            "status",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(bed_qc_rows)

edge_output = (
    table_dir
    / "macrosynteny_edges.tsv"
)

with edge_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "upper_species",
            "lower_species",
            "comparison_id",
            "source_anchor_file",
            "plot_anchor_file",
            "status",
        ],
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writeheader()
    writer.writerows(edge_rows)

# ============================================================
# Write JCVI seqids file
# ============================================================

seqids_file = (
    config_dir
    / "all_species.numeric.seqids"
)

with seqids_file.open(
    "w",
    encoding="utf-8",
) as handle:
    for code in species_order:
        handle.write(
            ",".join(
                species_seqids[code]
            )
            + "\n"
        )

# ============================================================
# Write JCVI layout file
#
# Layout columns:
# y, xstart, xend, rotation, color, label, va, bed
# ============================================================

track_y = {
    "PMAJ": 0.92,
    "VPAN": 0.82,
    "VSCU": 0.72,
    "VANA": 0.62,
    "VARV": 0.52,
    "VVER": 0.42,
    "VPER": 0.32,
    "VSER": 0.22,
    "VTRI": 0.12,
}

track_colors = {
    "PMAJ": "#6B7280",
    "VPAN": "#2563EB",
    "VSCU": "#2563EB",
    "VANA": "#D97706",
    "VARV": "#2563EB",
    "VVER": "#2563EB",
    "VPER": "#D97706",
    "VSER": "#2563EB",
    "VTRI": "#2563EB",
}

layout_file = (
    config_dir
    / "all_species.numeric.layout"
)

with layout_file.open(
    "w",
    encoding="utf-8",
) as handle:
    handle.write(
        "# y, xstart, xend, rotation, "
        "color, label, va, bed\n"
    )

    for code in species_order:
        handle.write(
            ",".join(
                [
                    f"{track_y[code]:.2f}",
                    "0.08",
                    "0.92",
                    "0",
                    track_colors[code],
                    scientific_labels[code],
                    "top",
                    str(
                        bed_dir
                        / f"{code}.numeric.bed"
                    ),
                ]
            )
            + "\n"
        )

    handle.write(
        "# edges\n"
    )

    for upper_index in range(
        len(species_order) - 1
    ):
        lower_index = upper_index + 1

        upper_species = (
            species_order[upper_index]
        )
        lower_species = (
            species_order[lower_index]
        )

        matching_edges = [
            row
            for row in edge_rows
            if row["upper_species"]
            == upper_species
            and row["lower_species"]
            == lower_species
        ]

        if len(matching_edges) != 1:
            raise SystemExit(
                "ERROR: Expected exactly one edge "
                f"for {upper_species} -> "
                f"{lower_species}; found "
                f"{len(matching_edges)}"
            )

        anchor_file = matching_edges[0][
            "plot_anchor_file"
        ]

        handle.write(
            f"e,{upper_index},{lower_index},"
            f"{anchor_file}\n"
        )

# ============================================================
# Write plotting order
# ============================================================

order_output = (
    table_dir
    / "macrosynteny_species_order.tsv"
)

with order_output.open(
    "w",
    newline="",
    encoding="utf-8",
) as handle:
    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(
        [
            "track_index",
            "species_code",
            "scientific_name",
            "ploidy",
            "chromosome_count",
        ]
    )

    for index, code in enumerate(
        species_order
    ):
        writer.writerow(
            [
                index,
                code,
                scientific_labels[code],
                manifest[code]["ploidy"],
                len(species_seqids[code]),
            ]
        )

print("Numeric chromosome mapping completed.")
print("Species tracks:", len(species_order))
print("Adjacent macrosynteny edges:", len(edge_rows))
PY

# ============================================================
# Display mapping before plotting
# ============================================================

echo
echo "============================================================"
echo "Numeric chromosome BED QC"
echo "============================================================"

column -t -s $'\t' \
    "${TABLEDIR}/numeric_bed_qc.tsv"

echo
echo "============================================================"
echo "Species plotting order"
echo "============================================================"

column -t -s $'\t' \
    "${TABLEDIR}/macrosynteny_species_order.tsv"

echo
echo "============================================================"
echo "Adjacent macrosynteny edges"
echo "============================================================"

column -t -s $'\t' \
    "${TABLEDIR}/macrosynteny_edges.tsv"

# ============================================================
# Validate numeric BED outputs
# ============================================================

for CODE in \
    PMAJ \
    VPAN \
    VSCU \
    VANA \
    VARV \
    VVER \
    VPER \
    VSER \
    VTRI
do
    BED_FILE="${BEDDIR}/${CODE}.numeric.bed"

    if [[ ! -s "${BED_FILE}" ]]; then
        echo "ERROR: Numeric BED missing or empty:" >&2
        echo "${BED_FILE}" >&2
        exit 1
    fi

    NON_NUMERIC_CHROMOSOMES=$(
        awk -F'\t' '
            $1 !~ /^[0-9]+$/ {
                count++
            }
            END {
                print count + 0
            }
        ' "${BED_FILE}"
    )

    if [[ "${NON_NUMERIC_CHROMOSOMES}" -ne 0 ]]; then
        echo "ERROR: Non-numeric chromosome IDs remain in ${CODE}." >&2
        exit 1
    fi
done

# ============================================================
# Detect supported JCVI karyotype options
# ============================================================

HELP_FILE="${LOGDIR}/jcvi_karyotype_help.txt"

python -m jcvi.graphics.karyotype --help \
    > "${HELP_FILE}" \
    2>&1 || true

EXTRA_OPTIONS=()

if grep -q -- '--notex' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--notex")
fi

if grep -q -- '--shadestyle' "${HELP_FILE}"; then
    EXTRA_OPTIONS+=("--shadestyle=curve")
fi

# ============================================================
# Generate vector PDF
# ============================================================

SEQIDS_FILE="${CONFIGDIR}/all_species.numeric.seqids"
LAYOUT_FILE="${CONFIGDIR}/all_species.numeric.layout"

PDF_FILE="${PLOTDIR}/all_9_species.macrosynteny.numeric_chromosomes.pdf"
PNG_FILE="${PLOTDIR}/all_9_species.macrosynteny.numeric_chromosomes.png"
SVG_FILE="${PLOTDIR}/all_9_species.macrosynteny.numeric_chromosomes.svg"

set +e

python -m jcvi.graphics.karyotype \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    --format=pdf \
    --outfile="${PDF_FILE}" \
    --figsize=16x14 \
    --dpi=300 \
    "${EXTRA_OPTIONS[@]}" \
    > "${LOGDIR}/macrosynteny_pdf.log" \
    2>&1

PDF_EXIT=$?

python -m jcvi.graphics.karyotype \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    --format=png \
    --outfile="${PNG_FILE}" \
    --figsize=16x14 \
    --dpi=400 \
    "${EXTRA_OPTIONS[@]}" \
    > "${LOGDIR}/macrosynteny_png.log" \
    2>&1

PNG_EXIT=$?

python -m jcvi.graphics.karyotype \
    "${SEQIDS_FILE}" \
    "${LAYOUT_FILE}" \
    --format=svg \
    --outfile="${SVG_FILE}" \
    --figsize=16x14 \
    "${EXTRA_OPTIONS[@]}" \
    > "${LOGDIR}/macrosynteny_svg.log" \
    2>&1

SVG_EXIT=$?

set -e

if [[ "${PDF_EXIT}" -ne 0 || ! -s "${PDF_FILE}" ]]; then
    echo "ERROR: PDF macrosynteny plot failed." >&2
    cat "${LOGDIR}/macrosynteny_pdf.log" >&2
    exit 1
fi

if [[ "${PNG_EXIT}" -ne 0 || ! -s "${PNG_FILE}" ]]; then
    echo "ERROR: PNG macrosynteny plot failed." >&2
    cat "${LOGDIR}/macrosynteny_png.log" >&2
    exit 1
fi

# SVG is useful for Inkscape but is optional if the installed
# Matplotlib backend does not support it through JCVI.
if [[ "${SVG_EXIT}" -ne 0 || ! -s "${SVG_FILE}" ]]; then
    echo "WARNING: Direct SVG output failed."
    echo "The vector PDF remains available for Inkscape."
    rm -f "${SVG_FILE}"
    SVG_STATUS="OPTIONAL_SVG_FAILED"
else
    SVG_STATUS="PASS"
fi

# ============================================================
# Final summary
# ============================================================

MAPPING_ROWS=$(
    awk -F'\t' '
        NR > 1 {
            count++
        }
        END {
            print count + 0
        }
    ' "${TABLEDIR}/numeric_chromosome_mapping.tsv"
)

EXPECTED_CHROMOSOMES=86

if [[ "${MAPPING_ROWS}" -ne "${EXPECTED_CHROMOSOMES}" ]]; then
    echo "ERROR: Expected ${EXPECTED_CHROMOSOMES} chromosome mappings;" >&2
    echo "observed ${MAPPING_ROWS}." >&2
    exit 1
fi

PDF_SIZE=$(
    stat -c '%s' "${PDF_FILE}"
)

PNG_SIZE=$(
    stat -c '%s' "${PNG_FILE}"
)

if [[ -s "${SVG_FILE}" ]]; then
    SVG_SIZE=$(
        stat -c '%s' "${SVG_FILE}"
    )
else
    SVG_SIZE=0
fi

SUMMARY="${TABLEDIR}/macrosynteny_plot_summary.tsv"

cat > "${SUMMARY}" <<EOF2
metricvalue
species_tracks9
adjacent_synteny_edges8
numeric_chromosome_mappings${MAPPING_ROWS}
pdf_file${PDF_FILE}
pdf_size_bytes${PDF_SIZE}
png_file${PNG_FILE}
png_size_bytes${PNG_SIZE}
svg_file${SVG_FILE}
svg_size_bytes${SVG_SIZE}
svg_status${SVG_STATUS}
statusPASS
EOF2

echo
echo "============================================================"
echo "All-species macrosynteny summary"
echo "============================================================"

column -t -s $'\t' \
    "${SUMMARY}"

# ============================================================
# Checkpoint
# ============================================================

cp -f \
    "${SUMMARY}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/numeric_chromosome_mapping.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/numeric_bed_qc.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/macrosynteny_species_order.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${TABLEDIR}/macrosynteny_edges.tsv" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${SEQIDS_FILE}" \
    "${CHECKPOINT_DIR}/"

cp -f \
    "${LAYOUT_FILE}" \
    "${CHECKPOINT_DIR}/"

cat > "${CHECKPOINT_DIR}/ALL_SPECIES_MACROSYNTENY_COMPLETE.txt" <<EOF2
checkpoint=all_species_macrosynteny_numeric_chromosomes
date=$(date --iso-8601=seconds)
species_tracks=9
adjacent_synteny_edges=8
numeric_chromosome_mappings=${MAPPING_ROWS}
pdf_status=PASS
png_status=PASS
svg_status=${SVG_STATUS}
status=PASS
next_step=inspect_and_refine_macrosynteny_figure
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
echo "Step 34E completed successfully"
echo "============================================================"
echo "Vector PDF:"
echo "${PDF_FILE}"
echo
echo "High-resolution PNG:"
echo "${PNG_FILE}"

if [[ -s "${SVG_FILE}" ]]; then
    echo
    echo "SVG:"
    echo "${SVG_FILE}"
fi
