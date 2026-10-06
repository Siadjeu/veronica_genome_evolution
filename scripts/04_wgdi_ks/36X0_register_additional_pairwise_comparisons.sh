#!/bin/bash
#SBATCH --partition=all_cpu.p
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=2000
#SBATCH --job-name=wgdi36X0
#SBATCH --output=11_wgdi/logs/step36X0_%j.out
#SBATCH --error=11_wgdi/logs/step36X0_%j.err

set -euo pipefail

PROJECT_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "${PROJECT_ROOT}"

WGDI_ROOT="11_wgdi"

ADMIN_DIR="${WGDI_ROOT}/00_admin"
QC_DIR="${WGDI_ROOT}/02_qc/additional_pairwise_comparisons"
CHECKPOINT_DIR="${WGDI_ROOT}/checkpoints/step36X0"

MANIFEST="${ADMIN_DIR}/step36X_additional_pairwise_comparisons.tsv"
INPUT_QC="${QC_DIR}/step36X0_input_qc.tsv"
SCRIPT_INVENTORY="${QC_DIR}/step36X0_existing_script_inventory.tsv"
OUTPUT_INVENTORY="${QC_DIR}/step36X0_existing_output_inventory.tsv"

mkdir -p \
    "${ADMIN_DIR}" \
    "${QC_DIR}" \
    "${CHECKPOINT_DIR}" \
    "${WGDI_ROOT}/logs"

###############################################################################
# TARGETED CLEANUP
###############################################################################

rm -f \
    "${MANIFEST}" \
    "${INPUT_QC}" \
    "${SCRIPT_INVENTORY}" \
    "${OUTPUT_INVENTORY}" \
    "${CHECKPOINT_DIR}/STEP36X0_COMPLETE.txt" \
    "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# NEW COMPARISONS
###############################################################################

cat > "${MANIFEST}" <<'TSV'
comparisonspecies1species2comparison_classscientific_rolestatus
VSCU_VSERVSCUVSERdiploid_diploidindependent_diploid_controlREGISTERED
VSCU_VPANVSCUVPANdiploid_diploidindependent_diploid_controlREGISTERED
VPAN_VPERVPANVPERdiploid_tetraploidindependent_polyploid_depth_testREGISTERED
TSV

###############################################################################
# VERIFY EXISTING PREPARED INPUTS
###############################################################################

printf \
    'species\tlens\tgff\tpep\tcds\tstatus\n' \
    > "${INPUT_QC}"

SPECIES=(
    VSCU
    VSER
    VPAN
    VPER
)

FAIL=0

for SP in "${SPECIES[@]}"
do
    LENS="${WGDI_ROOT}/01_inputs/lens/${SP}.wgdi.lens"

    GFF_CANDIDATE="$(
        find "${WGDI_ROOT}/01_inputs" \
            -type f \
            \( -iname "${SP}*.gff" \
               -o -iname "${SP}*.gff3" \
               -o -iname "${SP}*.wgdi.gff" \) \
            2>/dev/null \
            | sort \
            | head -n 1
    )"

    PEP_CANDIDATE="$(
        find "${WGDI_ROOT}/01_inputs" \
            -type f \
            \( -iname "${SP}*.pep*" \
               -o -iname "${SP}*.aa*" \
               -o -iname "${SP}*.protein*" \) \
            2>/dev/null \
            | sort \
            | head -n 1
    )"

    CDS_CANDIDATE="$(
        find "${WGDI_ROOT}/01_inputs" \
            -type f \
            -iname "${SP}*cds*" \
            2>/dev/null \
            | sort \
            | head -n 1
    )"

    STATUS="PASS"

    if [[ ! -s "${LENS}" ]]; then
        STATUS="FAIL"
        FAIL=1
    fi

    if [[ -z "${GFF_CANDIDATE}" || ! -s "${GFF_CANDIDATE}" ]]; then
        STATUS="FAIL"
        FAIL=1
    fi

    if [[ -z "${PEP_CANDIDATE}" || ! -s "${PEP_CANDIDATE}" ]]; then
        STATUS="FAIL"
        FAIL=1
    fi

    if [[ -z "${CDS_CANDIDATE}" || ! -s "${CDS_CANDIDATE}" ]]; then
        STATUS="FAIL"
        FAIL=1
    fi

    printf \
        '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${SP}" \
        "${LENS}" \
        "${GFF_CANDIDATE:-MISSING}" \
        "${PEP_CANDIDATE:-MISSING}" \
        "${CDS_CANDIDATE:-MISSING}" \
        "${STATUS}" \
        >> "${INPUT_QC}"
done

###############################################################################
# INVENTORY EXACT SCRIPTS USED IN PREVIOUS SUCCESSFUL WORKFLOW
###############################################################################

printf \
    'script\tbytes\n' \
    > "${SCRIPT_INVENTORY}"

find 05_scripts \
    -maxdepth 1 \
    -type f \
    \( \
        -name '36B*' \
        -o -name '36C*' \
        -o -name '36D*' \
        -o -name '36E*' \
        -o -name '36F*' \
        -o -name '36G*' \
        -o -name '36H*' \
    \) \
    -print0 \
    | sort -z \
    | while IFS= read -r -d '' FILE
do
    printf \
        '%s\t%s\n' \
        "${FILE}" \
        "$(stat -c '%s' "${FILE}")"
done \
    >> "${SCRIPT_INVENTORY}"

###############################################################################
# INVENTORY EXISTING OUTPUTS FOR THE THREE NEW COMPARISONS
#
# This is observational only. Nothing is deleted.
###############################################################################

printf \
    'comparison\tpath\tbytes\n' \
    > "${OUTPUT_INVENTORY}"

while IFS=$'\t' read -r \
    COMPARISON \
    SP1 \
    SP2 \
    CLASS \
    ROLE \
    STATUS
do
    [[ "${COMPARISON}" == "comparison" ]] && continue

    find "${WGDI_ROOT}" \
        -type f \
        -iname "*${COMPARISON}*" \
        -print0 \
        2>/dev/null \
        | sort -z \
        | while IFS= read -r -d '' FILE
    do
        printf \
            '%s\t%s\t%s\n' \
            "${COMPARISON}" \
            "${FILE}" \
            "$(stat -c '%s' "${FILE}")"
    done \
        >> "${OUTPUT_INVENTORY}"

done < "${MANIFEST}"

###############################################################################
# VERIFY ALL FOUR SPECIES PASSED INPUT QC
###############################################################################

PASS_COUNT="$(
    awk -F $'\t' '
        NR > 1 && $NF == "PASS" {
            n++
        }
        END {
            print n + 0
        }
    ' "${INPUT_QC}"
)"

if [[ "${PASS_COUNT}" -ne 4 ]]; then
    echo "ERROR: Only ${PASS_COUNT}/4 species passed prepared-input QC." >&2
    echo "Inspect: ${INPUT_QC}" >&2
    exit 1
fi

###############################################################################
# CHECKPOINT
###############################################################################

cat > "${CHECKPOINT_DIR}/STEP36X0_COMPLETE.txt" <<EOF2
checkpoint=step36X0_register_additional_pairwise_comparisons
date=$(date --iso-8601=seconds)
new_comparisons=3
comparison_1=VSCU_VSER
comparison_2=VSCU_VPAN
comparison_3=VPAN_VPER
VSCU_VSER_role=diploid_diploid_control
VSCU_VPAN_role=diploid_diploid_control
VPAN_VPER_role=diploid_tetraploid_test
prepared_species=4
prepared_input_qc_pass=4
existing_outputs_deleted=false
existing_workflow_scripts_modified=false
manifest=11_wgdi/00_admin/step36X_additional_pairwise_comparisons.tsv
input_qc=11_wgdi/02_qc/additional_pairwise_comparisons/step36X0_input_qc.tsv
script_inventory=11_wgdi/02_qc/additional_pairwise_comparisons/step36X0_existing_script_inventory.tsv
output_inventory=11_wgdi/02_qc/additional_pairwise_comparisons/step36X0_existing_output_inventory.tsv
status=PASS
next_step=step36X1_pairwise_DIAMOND_using_original_step36B_parameters
EOF2

cp -f \
    "${MANIFEST}" \
    "${INPUT_QC}" \
    "${SCRIPT_INVENTORY}" \
    "${OUTPUT_INVENTORY}" \
    "${CHECKPOINT_DIR}/"

find "${CHECKPOINT_DIR}" \
    -maxdepth 1 \
    -type f \
    ! -name 'sha256_checksums.txt' \
    -print0 \
    | sort -z \
    | xargs -0 sha256sum \
    > "${CHECKPOINT_DIR}/sha256_checksums.txt"

###############################################################################
# DISPLAY
###############################################################################

echo
echo "============================================================"
echo "NEW COMPARISONS"
echo "============================================================"

column -t -s $'\t' \
    "${MANIFEST}"

echo
echo "============================================================"
echo "INPUT QC"
echo "============================================================"

column -t -s $'\t' \
    "${INPUT_QC}"

echo
echo "============================================================"
echo "EXISTING RELEVANT SCRIPTS"
echo "============================================================"

column -t -s $'\t' \
    "${SCRIPT_INVENTORY}"

echo
echo "============================================================"
echo "ALREADY EXISTING FILES FOR NEW COMPARISONS"
echo "============================================================"

column -t -s $'\t' \
    "${OUTPUT_INVENTORY}"

echo
echo "============================================================"
echo "CHECKPOINT"
echo "============================================================"

cat \
    "${CHECKPOINT_DIR}/STEP36X0_COMPLETE.txt"

echo
echo "Step 36X0 completed successfully."
