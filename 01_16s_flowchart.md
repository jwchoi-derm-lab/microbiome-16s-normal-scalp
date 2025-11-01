# Microbiome 16S Normal Scalp Analysis Workflow

## Overview
Automated analysis pipeline for scalp 16S rRNA amplicon data from NCBI SRA through functional prediction.

**Dataset**: 1,190 scalp samples from 1,112 subjects from 17 BioProjects, 6,455 runs

---

## Methods

### 1. Data Acquisition
- **Source**: NCBI SRA (keywords: "scalp," "skin microbiome," "16S")
- **Processing**: SRA Toolkit → FASTQ (4,700 PE + 652 SE files)
- **Filtering**: Excluded lesional, treated, non-baseline, and fungal datasets

### 2. Quality Control & Preprocessing
- **QC**: FastQC v0.12.1 + MultiQC v1.31 (Phred ≥25)
- **Trimming**: fastp + Cutadapt v5.1
- **Primer detection**: 
  - Metadata/publication review (primary)
  - Computational inference: 1,500 reads/BioProject, ≥80% match in terminal 30 bp
- **Primer removal**: Cutadapt (-g/-u flags)

### 3. ASV Inference
- **Tool**: DADA2 v1.30.0 (per-BioProject)
- **Truncation**: PE (≥20 bp overlap), SE (≥75% amplicon)
- **Chimera removal**: consensus method
- **Taxonomy**: SILVA v138.2 (genus/species, naïve Bayesian classifier)

### 4. Analysis
- **Filtering**: ≥50 reads, ≥1% prevalence, >0.1% relative abundance within a BioProject
- **Alpha diversity**: Rarefied to 500 reads, Kruskal-Wallis test
- **Beta diversity**: Bray-Curtis dissimilarity, PCoA, PERMANOVA
- **Normalization**: BioProject-specific, integrated into phyloseq
- **Transformation**: Centered log-ratio (CLR)
- **Functional prediction**: PICRUSt2 v2.6.2 (EC, KEGG, MetaCyc)

---


