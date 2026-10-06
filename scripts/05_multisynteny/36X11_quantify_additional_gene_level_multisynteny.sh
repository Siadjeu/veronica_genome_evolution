#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=wgdi36X11
#SBATCH --output=11_wgdi/logs/step36X11_%j.out
#SBATCH --error=11_wgdi/logs/step36X11_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="${PROJECT_ROOT}/11_wgdi"

BLOCKINFO_DIR="${WGDI_ROOT}/07_blockinfo/02_results"
LENS_DIR="${WGDI_ROOT}/01_inputs/lens"

OUT_ROOT="${WGDI_ROOT}/12_multisynteny_depth_additional"
PAIR_DIR="${OUT_ROOT}/01_unique_gene_pairs"
GENE_DIR="${OUT_ROOT}/02_gene_partner_depth"
CHR_DIR="${OUT_ROOT}/03_chromosome_partner_depth"
SUMMARY_DIR="${OUT_ROOT}/04_summary_tables"

QC_DIR="${WGDI_ROOT}/02_qc/additional_multisynteny_depth"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X11"
LOG_DIR="${WGDI_ROOT}/logs"

SUMMARY="${SUMMARY_DIR}/step36X11_gene_level_multisynteny_summary.tsv"
CHR_SUMMARY="${SUMMARY_DIR}/step36X11_chromosome_partner_summary.tsv"
MANIFEST="${WGDI_ROOT}/00_admin/step36X11_multisynteny_manifest.tsv"
QC_TABLE="${QC_DIR}/step36X11_input_qc.tsv"

mkdir -p \
    "${PAIR_DIR}" \
    "${GENE_DIR}" \
    "${CHR_DIR}" \
    "${SUMMARY_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${LOG_DIR}"

###############################################################################
# Targeted cleanup only
###############################################################################

rm -f \
    "${SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${MANIFEST}" \
    "${QC_TABLE}" \
    "${CHECKPOINT_DIR}/STEP36X11_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

for comparison in VSCU_VSER VSCU_VPAN VPAN_VPER
do
    rm -f \
        "${PAIR_DIR}/${comparison}.unique_gene_pairs.tsv" \
        "${GENE_DIR}/${comparison}.gene_partner_depth.tsv" \
        "${CHR_DIR}/${comparison}.chromosome_partner_depth.tsv"
done

###############################################################################
# Main analysis
###############################################################################

python - \
    "${PROJECT_ROOT}" \
    "${BLOCKINFO_DIR}" \
    "${LENS_DIR}" \
    "${PAIR_DIR}" \
    "${GENE_DIR}" \
    "${CHR_DIR}" \
    "${SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${MANIFEST}" \
    "${QC_TABLE}" <<'PY'
from __future__ import annotations

import ast
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

project_root = Path(sys.argv[1])
blockinfo_dir = Path(sys.argv[2])
lens_dir = Path(sys.argv[3])
pair_dir = Path(sys.argv[4])
gene_dir = Path(sys.argv[5])
chr_dir = Path(sys.argv[6])
summary_path = Path(sys.argv[7])
chr_summary_path = Path(sys.argv[8])
manifest_path = Path(sys.argv[9])
qc_path = Path(sys.argv[10])

COMPARISONS = [
    {
        "comparison": "VSCU_VSER",
        "reference": "VSCU",
        "target": "VSER",
        "comparison_class": "diploid_diploid",
        "purpose": "independent_diploid_control",
    },
    {
        "comparison": "VSCU_VPAN",
        "reference": "VSCU",
        "target": "VPAN",
        "comparison_class": "diploid_diploid",
        "purpose": "independent_diploid_control",
    },
    {
        "comparison": "VPAN_VPER",
        "reference": "VPAN",
        "target": "VPER",
        "comparison_class": "diploid_tetraploid",
        "purpose": "independent_polyploid_depth_test",
    },
]


def relative(path: Path) -> str:
    return str(path.relative_to(project_root))


def read_lens(path: Path):
    """
    WGDI lens:
      col1 = chromosome
      col2 = chromosome bp length
      col3 = ordered-gene count
    """
    chromosome_gene_counts = {}
    chromosome_order = []

    with path.open(encoding="utf-8-sig") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split()

            if len(fields) < 3:
                raise SystemExit(
                    f"ERROR: lens row <3 columns: "
                    f"{path}, line {line_number}"
                )

            chromosome = fields[0]

            try:
                gene_count = int(round(float(fields[2])))
            except ValueError as error:
                raise SystemExit(
                    f"ERROR: invalid ordered-gene count in "
                    f"{path}, line {line_number}"
                ) from error

            if gene_count <= 0:
                raise SystemExit(
                    f"ERROR: nonpositive ordered-gene count "
                    f"for {chromosome} in {path}"
                )

            if chromosome in chromosome_gene_counts:
                raise SystemExit(
                    f"ERROR: duplicate chromosome {chromosome} "
                    f"in {path}"
                )

            chromosome_gene_counts[chromosome] = gene_count
            chromosome_order.append(chromosome)

    if not chromosome_gene_counts:
        raise SystemExit(
            f"ERROR: no chromosomes parsed from {path}"
        )

    return chromosome_gene_counts, chromosome_order


def parse_gene_list(raw: str) -> list[str]:
    """
    Parse WGDI block1/block2 gene lists robustly.

    First try Python-list syntax. If that fails, use a conservative
    delimiter fallback.
    """
    raw = (raw or "").strip()

    if not raw:
        return []

    try:
        value = ast.literal_eval(raw)

        if isinstance(value, (list, tuple)):
            result = [
                str(item).strip()
                for item in value
                if str(item).strip()
            ]

            if result:
                return result

    except (ValueError, SyntaxError):
        pass

    cleaned = raw.strip("[](){}")

    parts = re.split(
        r"\s*[,;]\s*|\s+",
        cleaned,
    )

    result = []

    for item in parts:
        item = item.strip().strip("'\"")

        if item:
            result.append(item)

    return result


summary_rows = []
chr_summary_rows_all = []
manifest_rows = []
qc_rows = []

for cfg in COMPARISONS:

    comparison = cfg["comparison"]
    reference = cfg["reference"]
    target = cfg["target"]

    blockinfo = (
        blockinfo_dir
        / f"{comparison}.blockinfo.csv"
    )

    ref_lens = (
        lens_dir
        / f"{reference}.wgdi.lens"
    )

    target_lens = (
        lens_dir
        / f"{target}.wgdi.lens"
    )

    for required in [
        blockinfo,
        ref_lens,
        target_lens,
    ]:
        if (
            not required.is_file()
            or required.stat().st_size == 0
        ):
            raise SystemExit(
                f"ERROR: missing X11 input: {required}"
            )

    ref_gene_counts, ref_chr_order = read_lens(
        ref_lens
    )

    target_gene_counts, target_chr_order = read_lens(
        target_lens
    )

    reference_ordered_genes = sum(
        ref_gene_counts.values()
    )

    target_ordered_genes = sum(
        target_gene_counts.values()
    )

    with blockinfo.open(
        newline="",
        encoding="utf-8-sig",
    ) as handle:

        reader = csv.DictReader(handle)

        fieldnames = reader.fieldnames or []

        required_columns = {
            "chr1",
            "chr2",
            "block1",
            "block2",
        }

        missing_columns = (
            required_columns - set(fieldnames)
        )

        if missing_columns:
            raise SystemExit(
                f"ERROR: {comparison} blockinfo missing "
                f"columns: {sorted(missing_columns)}"
            )

        blocks = list(reader)

    if not blocks:
        raise SystemExit(
            f"ERROR: no blockinfo rows for {comparison}"
        )

    unique_pairs = set()

    reference_gene_chr = {}
    target_gene_chr = {}

    raw_pair_assignments = 0
    unequal_block_gene_lists = 0

    for block_index, block in enumerate(
        blocks,
        start=1,
    ):

        chr1 = block["chr1"]
        chr2 = block["chr2"]

        if chr1 not in ref_gene_counts:
            raise SystemExit(
                f"ERROR: {chr1} absent from "
                f"{reference} lens in {comparison}"
            )

        if chr2 not in target_gene_counts:
            raise SystemExit(
                f"ERROR: {chr2} absent from "
                f"{target} lens in {comparison}"
            )

        genes1 = parse_gene_list(
            block["block1"]
        )

        genes2 = parse_gene_list(
            block["block2"]
        )

        if len(genes1) != len(genes2):
            unequal_block_gene_lists += 1
            raise SystemExit(
                f"ERROR: unequal block1/block2 gene counts "
                f"in {comparison}, block row {block_index}: "
                f"{len(genes1)} versus {len(genes2)}"
            )

        if not genes1:
            continue

        for ref_gene, target_gene in zip(
            genes1,
            genes2,
        ):
            raw_pair_assignments += 1

            previous = reference_gene_chr.get(
                ref_gene
            )

            if (
                previous is not None
                and previous != chr1
            ):
                raise SystemExit(
                    f"ERROR: reference gene {ref_gene} "
                    f"mapped to both {previous} and {chr1}"
                )

            reference_gene_chr[ref_gene] = chr1

            previous = target_gene_chr.get(
                target_gene
            )

            if (
                previous is not None
                and previous != chr2
            ):
                raise SystemExit(
                    f"ERROR: target gene {target_gene} "
                    f"mapped to both {previous} and {chr2}"
                )

            target_gene_chr[target_gene] = chr2

            unique_pairs.add(
                (
                    ref_gene,
                    chr1,
                    target_gene,
                    chr2,
                )
            )

    if not unique_pairs:
        raise SystemExit(
            f"ERROR: no unique gene pairs for {comparison}"
        )

    pair_path = (
        pair_dir
        / f"{comparison}.unique_gene_pairs.tsv"
    )

    with pair_path.open(
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
                "comparison",
                "reference_species",
                "target_species",
                "reference_gene",
                "reference_chromosome",
                "target_gene",
                "target_chromosome",
            ]
        )

        for (
            ref_gene,
            ref_chr,
            target_gene,
            target_chr,
        ) in sorted(unique_pairs):

            writer.writerow(
                [
                    comparison,
                    reference,
                    target,
                    ref_gene,
                    ref_chr,
                    target_gene,
                    target_chr,
                ]
            )

    partners_by_reference = defaultdict(set)
    target_chromosomes_by_reference = defaultdict(set)

    ref_genes_by_chr_target_chr = defaultdict(set)

    for (
        ref_gene,
        ref_chr,
        target_gene,
        target_chr,
    ) in unique_pairs:

        partners_by_reference[ref_gene].add(
            target_gene
        )

        target_chromosomes_by_reference[
            ref_gene
        ].add(target_chr)

        ref_genes_by_chr_target_chr[
            (ref_chr, target_chr)
        ].add(ref_gene)

    gene_depth_path = (
        gene_dir
        / f"{comparison}.gene_partner_depth.tsv"
    )

    gene_depth_rows = []

    for ref_gene in sorted(
        partners_by_reference
    ):
        target_genes = sorted(
            partners_by_reference[ref_gene]
        )

        target_chromosomes = sorted(
            target_chromosomes_by_reference[
                ref_gene
            ]
        )

        gene_depth_rows.append(
            {
                "comparison": comparison,
                "reference_species": reference,
                "target_species": target,
                "reference_gene": ref_gene,
                "reference_chromosome": (
                    reference_gene_chr[
                        ref_gene
                    ]
                ),
                "unique_target_gene_partners": str(
                    len(target_genes)
                ),
                "distinct_target_chromosomes": str(
                    len(target_chromosomes)
                ),
                "target_genes": ";".join(
                    target_genes
                ),
                "target_chromosomes": ";".join(
                    target_chromosomes
                ),
            }
        )

    with gene_depth_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(
                gene_depth_rows[0]
            ),
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(gene_depth_rows)

    chr_depth_path = (
        chr_dir
        / f"{comparison}.chromosome_partner_depth.tsv"
    )

    chr_rows = []

    for ref_chr in ref_chr_order:

        for target_chr in target_chr_order:

            genes = ref_genes_by_chr_target_chr.get(
                (ref_chr, target_chr),
                set(),
            )

            if not genes:
                continue

            chr_rows.append(
                {
                    "comparison": comparison,
                    "reference_species": reference,
                    "target_species": target,
                    "reference_chromosome": ref_chr,
                    "target_chromosome": target_chr,
                    "reference_genes_with_partner": str(
                        len(genes)
                    ),
                    "reference_chromosome_ordered_genes": str(
                        ref_gene_counts[
                            ref_chr
                        ]
                    ),
                    "fraction_reference_chromosome": (
                        f"{len(genes) / ref_gene_counts[ref_chr]:.8f}"
                    ),
                }
            )

            chr_summary_rows_all.append(
                chr_rows[-1]
            )

    with chr_depth_path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        if chr_rows:
            writer = csv.DictWriter(
                handle,
                fieldnames=list(
                    chr_rows[0]
                ),
                delimiter="\t",
                lineterminator="\n",
            )

            writer.writeheader()
            writer.writerows(chr_rows)

    any_partner = len(
        partners_by_reference
    )

    one_partner = sum(
        len(partners) == 1
        for partners in partners_by_reference.values()
    )

    two_partners = sum(
        len(partners) == 2
        for partners in partners_by_reference.values()
    )

    more_than_two = sum(
        len(partners) > 2
        for partners in partners_by_reference.values()
    )

    exactly_two_target_chr = sum(
        len(chromosomes) == 2
        for chromosomes
        in target_chromosomes_by_reference.values()
    )

    at_least_two_target_chr = sum(
        len(chromosomes) >= 2
        for chromosomes
        in target_chromosomes_by_reference.values()
    )

    inferred_zero = (
        reference_ordered_genes
        - any_partner
    )

    if inferred_zero < 0:
        raise SystemExit(
            f"ERROR: syntenic reference genes exceed "
            f"reference lens gene count for {comparison}"
        )

    summary_rows.append(
        {
            "comparison": comparison,
            "comparison_class": (
                cfg["comparison_class"]
            ),
            "purpose": cfg["purpose"],
            "reference_species": reference,
            "target_species": target,
            "reference_ordered_genes": str(
                reference_ordered_genes
            ),
            "target_ordered_genes": str(
                target_ordered_genes
            ),
            "strict_blockinfo_rows": str(
                len(blocks)
            ),
            "raw_gene_pair_assignments": str(
                raw_pair_assignments
            ),
            "unique_gene_pairs": str(
                len(unique_pairs)
            ),
            "reference_genes_with_any_partner": str(
                any_partner
            ),
            "reference_genes_inferred_zero_partner": str(
                inferred_zero
            ),
            "reference_genes_with_1_partner": str(
                one_partner
            ),
            "reference_genes_with_2_partners": str(
                two_partners
            ),
            "reference_genes_with_gt2_partners": str(
                more_than_two
            ),
            "reference_genes_exactly_2_target_chromosomes": str(
                exactly_two_target_chr
            ),
            "reference_genes_ge2_target_chromosomes": str(
                at_least_two_target_chr
            ),
            "fraction_reference_with_any_partner": (
                f"{any_partner / reference_ordered_genes:.8f}"
            ),
            "fraction_1_partner_among_syntenic": (
                f"{one_partner / any_partner:.8f}"
            ),
            "fraction_2_partners_among_syntenic": (
                f"{two_partners / any_partner:.8f}"
            ),
            "fraction_gt2_partners_among_syntenic": (
                f"{more_than_two / any_partner:.8f}"
            ),
            "fraction_exactly_2_target_chr_among_syntenic": (
                f"{exactly_two_target_chr / any_partner:.8f}"
            ),
            "fraction_ge2_target_chr_among_syntenic": (
                f"{at_least_two_target_chr / any_partner:.8f}"
            ),
            "status": "PASS",
        }
    )

    qc_rows.append(
        {
            "comparison": comparison,
            "blockinfo": relative(blockinfo),
            "reference_lens": relative(
                ref_lens
            ),
            "target_lens": relative(
                target_lens
            ),
            "blockinfo_rows": str(
                len(blocks)
            ),
            "unequal_block_gene_lists": str(
                unequal_block_gene_lists
            ),
            "unique_gene_pairs": str(
                len(unique_pairs)
            ),
            "reference_gene_chr_conflicts": "0",
            "target_gene_chr_conflicts": "0",
            "status": "PASS",
        }
    )

    manifest_rows.append(
        {
            "comparison": comparison,
            "comparison_class": (
                cfg["comparison_class"]
            ),
            "purpose": cfg["purpose"],
            "reference_species": reference,
            "target_species": target,
            "source_blockinfo": relative(
                blockinfo
            ),
            "reference_lens": relative(
                ref_lens
            ),
            "target_lens": relative(
                target_lens
            ),
            "unique_gene_pairs": relative(
                pair_path
            ),
            "gene_partner_depth": relative(
                gene_depth_path
            ),
            "chromosome_partner_depth": relative(
                chr_depth_path
            ),
            "status": "PASS",
        }
    )

    print(
        f"{comparison}: "
        f"{len(blocks)} blocks; "
        f"{len(unique_pairs)} unique gene pairs; "
        f"{any_partner} reference genes with partners; "
        f"{at_least_two_target_chr} with >=2 target chromosomes.",
        flush=True,
    )


def write_tsv(path, rows):
    if not rows:
        raise SystemExit(
            f"ERROR: no rows generated for {path}"
        )

    with path.open(
        "w",
        newline="",
        encoding="utf-8",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(rows)


write_tsv(
    summary_path,
    summary_rows,
)

write_tsv(
    chr_summary_path,
    chr_summary_rows_all,
)

write_tsv(
    manifest_path,
    manifest_rows,
)

write_tsv(
    qc_path,
    qc_rows,
)

if len(summary_rows) != 3:
    raise SystemExit(
        f"ERROR: expected 3 comparisons; "
        f"found {len(summary_rows)}"
    )

if any(
    row["status"] != "PASS"
    for row in summary_rows
):
    raise SystemExit(
        "ERROR: at least one X11 comparison failed."
    )

print(
    "Step 36X11 gene-level multisynteny analysis: PASS"
)
PY

###############################################################################
# Final shell QC
###############################################################################

PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            n++
        }
        END {
            print n + 0
        }
    ' "${SUMMARY}"
)"

if [[ "${PASS_COUNT}" -ne 3 ]]; then
    echo "ERROR: ${PASS_COUNT}/3 X11 comparisons passed." >&2
    exit 1
fi

cat > "${CHECKPOINT_DIR}/STEP36X11_COMPLETE.txt" <<EOF2
checkpoint=step36X11_additional_gene_level_multisynteny_depth
date=$(date --iso-8601=seconds)
comparisons_expected=3
comparisons_complete=${PASS_COUNT}
comparisons=VSCU_VSER,VSCU_VPAN,VPAN_VPER
reference_orientation=VSCU_to_VSER,VSCU_to_VPAN,VPAN_to_VPER
source=full_strict_WGDI_blockinfo
gene_pair_deduplication=exact_reference_target_gene_pair
reference_denominator=WGDI_lens_column_3_ordered_gene_count
partner_depth=unique_target_gene_partners_per_reference_gene
chromosome_depth=distinct_target_chromosomes_per_reference_gene
zero_partner_genes=inferred_from_reference_lens_total
summary=11_wgdi/12_multisynteny_depth_additional/04_summary_tables/step36X11_gene_level_multisynteny_summary.tsv
chromosome_summary=11_wgdi/12_multisynteny_depth_additional/04_summary_tables/step36X11_chromosome_partner_summary.tsv
manifest=11_wgdi/00_admin/step36X11_multisynteny_manifest.tsv
status=PASS
next_step=step36X12_chromosome_consistent_multisynteny
EOF2

cp -f \
    "${SUMMARY}" \
    "${CHR_SUMMARY}" \
    "${MANIFEST}" \
    "${QC_TABLE}" \
    "${CHECKPOINT_DIR}/"

find \
    "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

echo
echo "============================================================"
echo "Step 36X11 gene-level multisynteny summary"
echo "============================================================"

column -t -s $'\t' "${SUMMARY}"

echo
echo "============================================================"
echo "Step 36X11 checkpoint"
echo "============================================================"

cat "${CHECKPOINT_DIR}/STEP36X11_COMPLETE.txt"
