# Figure 2 — Genome-wide macrosynteny

This directory contains the final JCVI-based workflow used to generate the
chromosome-scale macrosynteny analysis across eight *Veronica* species and the
*Plantago major* outgroup.

## Species

| Code | Species | Ploidy | Chromosomes |
|------|---------|--------|-------------|
| PMAJ | *Plantago major* | 2x | 6 |
| VPAN | *Veronica panormitana* | 2x | 9 |
| VSCU | *Veronica scutellata* | 2x | 9 |
| VANA | *Veronica anagallis-aquatica* | 4x | 18 |
| VARV | *Veronica arvensis* | 2x | 8 |
| VPER | *Veronica persica* | 4x | 14 |
| VSER | *Veronica serpyllifolia* | 2x | 7 |
| VTRI | *Veronica triphyllos* | 2x | 7 |
| VVER | *Veronica verna* | 2x | 8 |

## Software

The final analysis used:

- JCVI 1.6.5
- DIAMOND 2.2.4

One representative protein-coding transcript per locus was retained before
comparative synteny analysis to avoid redundancy caused by alternative
transcripts.

## Workflow

The scripts are ordered according to the final analysis workflow.

```text
Representative proteomes
        |
        v
Input discovery and validation
        |
        v
JCVI BED/protein preparation
        |
        v
Pairwise synteny
        |
        v
Tree-order comparison completion
        |
        v
Refined macrosynteny
        |
        v
Publication macrosynteny plot
