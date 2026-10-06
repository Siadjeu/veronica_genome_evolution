# Figure 4 — Ks distributions of syntenic blocks

This directory contains the final WGDI-based workflow used to calculate,
validate and visualize synonymous substitution-rate (Ks) distributions of
syntenic blocks across the comparative *Veronica* genome dataset.

## Main workflow

The analysis includes:

1. WGDI input validation
2. DIAMOND protein homology searches
3. Improved collinearity detection
4. Preparation of unique syntenic gene pairs
5. Chunked Ks estimation
6. WGDI block-information construction
7. Curation of block-level Ks datasets
8. Corrected multi-bandwidth Ks peak validation
9. Authentic WGDI KsFigure generation
10. Independent additional pairwise controls

## Figure 4 comparisons

Panel A contains the nine self-comparisons together with PMAJ–VSCU,
VSCU–VANA and VSCU–VPER.

Panel B contains the additional between-genome comparisons VSCU–VSER,
VSCU–VPAN and VPAN–VPER.

## Interpretation

The lower-Ks components observed in independent diploid–diploid comparisons
provide controls showing that these components largely represent species
divergence rather than whole-genome duplication.

An older broad Ks component is shared across multiple *Veronica* comparisons,
but its overlap with the PMAJ–VSCU orthologous-divergence distribution means
that Ks evidence alone is not treated as proof of a *Veronica*-specific
ancestral whole-genome duplication.

## Primary Ks metric

YN00 Ks estimates were used as the primary values for block-level analyses.
Only positive Ks values were used for peak analysis.

## Peak validation

The corrected peak-validation workflow used log-transformed positive
block-median Ks values between 0.1 and 3.0, multiple bandwidth factors and
bootstrap support.

The corrected `36F4R` workflow supersedes earlier peak-validation attempts.

## Software

- WGDI 0.75
- DIAMOND 2.2.4
- MAFFT 7.526
- PAML/yn00 4.10.10
