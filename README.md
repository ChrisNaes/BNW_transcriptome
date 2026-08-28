# Multi-organ transcriptomic responses to severe sarcoptic mange in bare-nosed wombats

This repository contains the analysis code and supporting annotation files for the manuscript:

**Extended multi-organ pathogenesis results from severe *Sarcoptes scabiei* infestation of wombats**

Authors: Christina Naesborg-Nielsen, T. Fraser, K. Mounsey, and S. Carver.

## Overview

Sarcoptic mange, caused by the ectoparasitic mite *Sarcoptes scabiei*, is primarily considered a disease of the skin. However, severe disease in bare-nosed wombats (*Vombatus ursinus*) is associated with substantial weight loss, physiological decline, and mortality.

This study used multi-tissue RNA sequencing to investigate transcriptional responses associated with advanced crusted mange across six tissues:

- Skin
- Lung
- Liver
- Lymph node
- Skeletal muscle
- Thyroid

RNA-seq data from clinically healthy and advanced mange-affected wombats were analyzed using a tissue-specific differential expression framework, followed by Gene Ontology annotation, cross-tissue comparisons, and targeted functional pathway analyses.

## Repository contents

### `scripts/`

Analysis scripts used for:

1. Raw read quality control
2. Adapter and quality trimming
3. HISAT2 alignment
4. BAM processing and quality control
5. Gene-level quantification with featureCounts
6. Differential expression analysis using DESeq2
7. Gene Ontology annotation
8. Cross-tissue transcriptional analyses
9. Figure generation

### `annotation/`

Contains the custom bare-nosed wombat gene annotation used for gene-level quantification.

RNA-seq reads were aligned to the bare-nosed wombat reference genome:

**vu-2k**  
NCBI accession: `GCA_028626985.1`

Because this genome lacks a comprehensive annotation, gene models were transferred from the annotated *Vombatus ursinus* assembly:

`GCF_900497805.2`

and supplemented with manually curated immune-associated loci.

See `annotation/README.md` for additional information.

### `metadata/`

Contains sample metadata required for the RNA-seq analyses.

### `results/`

Large analysis outputs are not stored directly in this repository.

Processed datasets supporting the manuscript are provided as Supplementary Information associated with the publication.

## RNA-seq data

Raw RNA-sequencing reads generated for this study are available through the NCBI Sequence Read Archive under:

**BioProject:** `[ADD BIOPROJECT ACCESSION]`

Publicly available bare-nosed wombat RNA-seq data used to support transcript reconstruction were obtained from:

**NCBI BioProject:** `PRJNA484583`

## Analysis workflow

The primary analysis workflow was:

```text
Raw FASTQ files
      |
      v
FastQC / MultiQC
      |
      v
TrimGalore / cutadapt
      |
      v
HISAT2 alignment
      |
      v
SAMtools BAM processing
      |
      v
featureCounts
      |
      v
DESeq2
      |
      v
Gene Ontology annotation
      |
      v
Cross-tissue and functional analyses
      |
      v
Figures and supplementary tables

## Differential expression analysis

Differential expression was analyzed independently within each tissue using DESeq2.

Genes were retained if they had at least 10 reads in at least two samples within the corresponding tissue.

Disease status was modeled as the explanatory variable, with healthy wombats used as the reference group.

Genes were classified as differentially expressed when:

Benjamini-Hochberg FDR < 0.05
and
|log2 fold-change| >= 2

Positive log2 fold-change values therefore indicate higher expression in wombats with advanced crusted mange relative to healthy wombats.

## Functional annotation

Gene-level differential expression results were linked to NCBI Gene Identifiers and Gene Ontology annotations.

GO-term transcriptional summaries were calculated using all genes for which both a log2 fold-change estimate and GO annotation were available.

These summaries describe the direction and magnitude of transcriptional regulation among genes associated with each GO term and should not be interpreted as conventional over-representation or gene-set enrichment tests.

Targeted functional modules were subsequently constructed using predefined GO-term keyword searches.

## Software

Major software used in the analysis included:

FastQC v0.11.9
MultiQC v1.14
TrimGalore v0.6.7
cutadapt v4.5
HISAT2 v2.2.2
SAMtools v1.18
StringTie v3.0.0
Subread / featureCounts v2.0.6
DESeq2 v1.52.0

Additional R package versions are recorded in the analysis environment/session information.

## Reproducibility

The complete analysis was rerun on a high-performance computing cluster using the scripts provided in this repository.

Where possible, intermediate and final outputs, software versions, command-line logs, and reference files were retained to provide a reproducible record of the analysis.

Large sequencing and alignment files are not stored on GitHub.

## Data availability

Raw RNA-sequencing data generated in this study have been deposited in the NCBI Sequence Read Archive under BioProject accession [ADD BIOPROJECT ACCESSION].

Public RNA-seq data used for transcript reconstruction are available under BioProject PRJNA484583.

The bare-nosed wombat reference genome is available under accession GCA_028626985.1.

Processed data supporting the study are provided with the manuscript as Supplementary Information.

## Citation

If you use code or resources from this repository, please cite:

Naesborg-Nielsen C., Fraser T., Mounsey K., Carver S.
Extended multi-organ pathogenesis results from severe Sarcoptes scabiei infestation of wombats.
[Journal and DOI to be added following publication]

## Contact

For questions regarding this repository or analysis, please contact:

Christina Naesborg-Nielsen
christina.naesborgnielsen@uga.edu
