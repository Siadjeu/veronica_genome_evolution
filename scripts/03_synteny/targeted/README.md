# Figure 3 — Targeted chromosome-scale synteny

This directory contains the final targeted JCVI-based workflow used to compare
chromosome relationships among diploid and tetraploid *Veronica* genomes.

## Figure 3 comparisons

- A: VSCU–VSER, diploid–diploid
- B: VSCU–VPAN, diploid–diploid
- C: VSCU–VPER, diploid–tetraploid
- D: VPAN–VPER, diploid–tetraploid

## Link interpretation

Blue links represent mapped full-anchor synteny pairs.

Red links represent chromosome pairs participating in multiple
strong-partner relationships.

## Final scripts

```text
35Y0_generate_PMAJ_VSCU_comparison.sh
35Y_generate_synteny_optimized_pairwise_plots.sh
35Y1_generate_VSCU_VSER_comparison.sh
35Y2_generate_additional_synteny_optimized_pairwise_plots.sh
35Z_create_final_jcvi_synteny_checkpoint.sh
## Notes
The optimized chromosome orders were used for visualization only and do not
alter the underlying synteny calls.
The Figure 3 comparisons were selected to contrast simpler diploid chromosome
relationships with the increased multisyntenic structure observed in
comparisons involving the tetraploid VPER genome.
