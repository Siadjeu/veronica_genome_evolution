# Veronica panormitana genome assembly and functional annotation

## Assembly statistics

Assembly statistics were calculated from the final chromosome-scale
Veronica panormitana assembly.

QUAST was used to summarize:

- total assembly size
- number of chromosome-scale scaffolds
- scaffold N50
- GC content
- ambiguous bases

Final assembly statistics:

| Statistic | Value |
|---|---:|
| Genome size (Mb) | 237.93 |
| Scaffolds | 9 |
| Scaffold N50 (Mb) | 26.61 |
| Pseudochromosomes | 9 |
| GC content (%) | 35.26 |

## Repeat annotation

Repeats were characterized using RepeatMasker v4.1.6 with
RMBlast v2.14.1+ and the species-specific repeat library
`vpanor-families_TEclass.fa`.

The final assembly contained 38.51% repetitive sequence.

## Structural annotation

Structural annotation was performed using the BRAKER4 workflow.

For V. panormitana, annotation incorporated:

- RNA-seq accession SRR3491905
- HISAT2 v2.2.1
- Viridiplantae OrthoDB v11 protein evidence
- BRAKER3
- AUGUSTUS
- GeneMark-ETP
- TSEBRA

The final annotation contained 28,325 predicted protein-coding genes.

One representative protein per locus was selected for downstream
comparative and functional analyses, yielding 28,313 representative
protein sequences.

## Swiss-Prot annotation

Representative proteins were searched against reviewed
UniProtKB/Swiss-Prot using DIAMOND v2.2.2.

Parameters included:

- sensitive mode
- E-value <= 1e-5
- maximum one target per query

A Swiss-Prot hit was considered supported when:

- amino-acid identity >= 30%
- query coverage >= 50%

Strict Swiss-Prot support:

17,798 / 28,313 proteins = 62.86%.

## InterProScan annotation

InterProScan v5.78-109.0 was run using:

- Gene3D v4.3.0
- PANTHER v19.0
- Pfam v38.2
- SMART v9.0
- SUPERFAMILY v1.75

InterProScan support:

25,045 / 28,313 proteins = 88.46%.

## Final functional annotation

Functional annotation coverage was calculated as the non-redundant
union of proteins supported by either stringent Swiss-Prot similarity
or at least one InterProScan assignment.

Results:

| Category | Proteins | Percentage |
|---|---:|---:|
| Swiss-Prot strict | 17,798 | 62.86 |
| InterProScan | 25,045 | 88.46 |
| Both | 17,698 | 62.51 |
| Swiss-Prot only | 100 | 0.35 |
| InterProScan only | 7,347 | 25.95 |
| Functional annotation union | 25,145 | 88.81 |
| Unsupported | 3,168 | 11.19 |

Final functional annotation coverage:

25,145 / 28,313 = 88.81%.
