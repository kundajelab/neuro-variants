# Marderstein*, Kundu* et al.

This repository contains all scripts for the manuscript **"Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart"** by **Marderstein, Kundu et al.**

For any questions, please contact:  

Andrew Marderstein & Soumya Kundu
📧 **mardera1@mskcc.org**, **soumyak@stanford.edu**  

---

# Pre-processing

These sets of scripts in the **`preprocess/`** directory generate all of the outputs used for the downstream analyses.

- **`0_process_data/`** – Scripts for processing the scATAC-seq data.
- **`1_train_chrombpnet/`** – Scripts for training the ChromBPNet models.
- **`2_score_variants/`** – Scripts for scoring the rare, common, and ASD variants.
- **`3_shap_variants/`** – Scripts for generating DeepLIFT / DeepSHAP contribution scores for both alleles of each variant.
- **`4_shap_peaks/`** - Scripts for generating DeepLIFT / DeepSHAP contribution scores for scATAC-seq peaks.
- **`5_run_modisco/`** - Scripts for running TF-MoDISco to identify the motif patterns learned by each model.
- **`6_cluster_motifs/`** - Scripts for running MotifCompendium to cluster the motif patterns from TF-MoDISco.
- **`7_run_finemo/`** - Scripts for running Fi-NeMo for identifying motif instances in the genome.

---

# Analysis

These sets of scripts in the **`analysis/`** directory generate all of the results presented in the manuscript.

---

## 1_ChromBPNet_PreprocessingAnalysis

These scripts process **ChromBPNet** variant scoring outputs and compile annotation tables for downstream analyses.

### Steps:
1. **Extract outputs**: Run `1_pull_scores.sh` to extract relevant ChromBPNet outputs.
2. **Analyze model performance**:  
   - `2a_model_performance.R` evaluates performance metrics.  
   - `2b_model_performance_plot.R` identifies model outliers.
3. **Annotate variants**:  
   - Use `3a_bed2vcf.Rare.CADD_VEP.R` to run **CADD** and **VEP**.
   - Process outputs with `3c_Process_CADD_VEP.R`.
4. **Merge results**: Run `4_mergeData.R` to integrate annotations, merge scores, and remove outliers.

---

## 2_Effects_Across_Contexts

These scripts correspond to the manuscript section **"Variants effects are shaped by genomic context and TF binding"**.  
They analyze:
- **Genomic context** – A variant’s proximity to transcribed regions.
- **Cell-type specificity** – How constrained or widespread variant effects are.
- **Regulatory magnitude** – The extent of chromatin accessibility and TF binding changes.

---

## 3_Application_to_GWAS_and_QTL_studies

These scripts support the manuscript sections:
- **"Context-specific models reveal regulatory effects of fine-mapped eQTLs"**
- **"Pinpointing disease-relevant variants using cell-type-specific chromatin models"**
- **"Microglia-driven mechanisms of Alzheimer’s disease risk"**

We use **ChromBPNet** to identify candidate causal variants affecting gene regulation and disease risk.

---

## 4_Rare_vs_Common

These scripts correspond to:
- **"Ultra-rare variants show larger and more shared regulatory effects than common variants"**
- **"Specific motifs influence constraint of fetal neuron regulation"**

They compare rare and common variant effects to understand the selective pressures that influence allele frequency distributions across human populations.

---

## 5_FLARE

These scripts correspond to:
- **"FLARE: a functional genomic model of constraint"**
- **"FLARE prioritizes de novo non-coding mutations in autism"**

Since **PhyloP** scores are not context-specific, **FLARE** models the relationship between genomic context, regulatory effects, and evolutionary conservation within **cell-type-specific contexts**. FLARE:
1. **Disentangles** accessibility and regulatory effects from conservation.
2. **Integrates** multiple functional genomic features into a unified model.
3. **Captures** regulatory potential across multiple cell types.

[We set up a FLARE repository for the FLARE method, which can be found by clicking here.](https://github.com/drewmard/FLARE)

---

## 6_FLARE_extended

These show additional FLARE applications, with scripts corresponding to:
- **"FLARE prioritizes de novo non-coding mutations in congenital heart disease"**
- **"FLARE prioritizes variants underlying expression outliers in the adult brain"**
- **"FLARE captures common variant heritability in schizophrenia"**

---

# Cite

Marderstein^, Kundu^, et al. Mapping the regulatory effects of common and rare non-coding variants across cellular and developmental contexts in the brain and heart.
