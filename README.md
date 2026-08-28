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
