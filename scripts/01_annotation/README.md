# Standardized genome annotation

This directory contains the final BRAKER4-based annotation workflow used for
eight *Veronica* genomes and the *Plantago major* outgroup.

## Annotation strategy

All genomes were reannotated using a standardized workflow to minimize bias
caused by heterogeneous source annotations.

### Repeat identification and masking

- RepeatModeler2 v2.0.9
- RepeatMasker v4.2.4
- Tandem Repeats Finder v4.09.1

Genome sequences shorter than 10 kb were excluded from structural gene
prediction.

### Protein evidence

Viridiplantae proteins from OrthoDB v11 were used for all species.

### Veronica panormitana

VPAN was annotated using both transcript and protein evidence (ETP mode).
RNA-seq accession SRR3491905 was aligned with HISAT2 v2.2.1 using --dta and
processed with SAMtools v1.13.

The workflow included GeneMark-ETP, GeneMarkS-T and DIAMOND.
StringTie2 v2.2.1 was used for transcript reconstruction/UTR support.

### Remaining genomes

The other seven *Veronica* genomes and *P. major* were annotated using the
protein-supported EP/BRAKER2 workflow.

The workflow included:

- GeneMark-ES/EP+
- ProtHint v2.6.0
- DIAMOND v2.0.15
- Spaln v2.3.3f
- AUGUSTUS v3.5.0
- TSEBRA

### Annotation completeness

BUSCO v6.0.0 with embryophyta_odb12 (2,026 orthologs) was used to assess
genome and predicted-protein completeness.

## Final scripts

- `19_prepare_braker4_genomes.sh`
- `21_run_braker4_array.sh`
- `21a_resume_VPAN_braker4.sh`
- `22_run_remaining_braker4_array.sh`
- `23_validate_all_braker4_annotations.sh`
- `25_make_final_manuscript_annotation_table.sh`
- `26_create_annotation_checkpoint.sh`

Representative proteomes used for comparative genomics were prepared
downstream using `scripts/03_synteny/24_prepare_representative_proteomes.sh`.
