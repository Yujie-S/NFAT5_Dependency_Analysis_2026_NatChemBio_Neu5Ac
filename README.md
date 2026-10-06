# NFAT5_Dependency_Analysis_2026_NatChemBio_Neu5Ac
R scripts for Neu5Ac RNA-seq DEG and NFAT5-dependency index analysis in manuscript Sialic acid augments therapeutic T cell fitness through cellular stress response (Nat Chem Bio, 2026) by Yujie Shi et al. 

Uploaded and last edited by Yujie Shi (yujieshi@scripps.edu) Oct-06-2026.

## Introduction

To analyze the dependency on NFAT5 of Neu5Ac-induced changes, we performed RNA-seq on three conditions: untreated scrambled control (Un_Scr), Neu5Ac-treated scrambled control (Neu_Scr), and Neu5Ac-treated Nfat5 knockdown (Neu_KD). Differential expression analysis was performed to calculate the log2 fold changes (LFC) and adjusted p values (padj) of all genes between Neu_Scr vs Un_Scr (Comparison 1, LFC_Scr, padj_Scr) and Neu_KD vs Un_Scr (Comparison 2, LFC_KD, padj_KD). Genes whose Neu5Ac-induced change is strongly NFAT5-dependent should show a clear reduction or loss of differential expression upon NFAT5 knockdown. In contrast, NFAT5-independent genes should maintain similar magnitude and significance in both comparisons, or can be in different directions due to other factors. With this rationale, we defined an Nfat5-Dependency Index for each gene as DI_gene = LFC_KD / LFC_Scr. This index quantifies how Nfat5 KD alters the magnitude and/or direction of Neu5Ac-induced changes. Differentially expressed genes (DEGs) from Comparison 1 (padj_Scr < 0.05 and |log2FC| > 0.25; n = 4,891) were classified as: Strongly NFAT5-dependent (0 ≤ DI < 0.68, or padj_KD ≥ 0.05), partially NFAT5-dependent (0.68 ≤ DI < 1), or NFAT5-independent (DI ≥ 1, or DI < 0 with padj_KD < 0.05). The threshold for strong or weak dependency (DI = 0.68) was empirically set based on the calculated DI of Nfat5 itself. 

## Requirements

R packages: DESeq2, apeglm, dplyr, readr, tibble, tidyr, ggplot2, ggrepel, pheatmap, RColorBrewer.

## Scripts

- `1_pre_process_DEG.R`: DESeq2 contrasts and apeglm-shrunken log2 fold changes.
  
- `2_DI_analysis.R`: dependency index, classification, and plots. Run after `1_pre_process_DEG.R`.

Log2 fold changes used for DI are apeglm-shrunken. Adjusted p values are from the unshrunken Wald test.

## Data files

- Raw RNA-seq data is deposited in GEO under accession [GSE347039](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE347039) where experimental details can be found.
  
- Input file:
  
    `gene_count.xls`: Raw count matrix of experimental groups (three samples per group):
    - `NC_Scr` T cells cultured in normal medium, shScramble control
    - `SA_Scr` T cells cultured in Neu5Ac-supplemented medium, shScramble control
    - `SA_KD` T cells cultured in Neu5Ac-supplemented medium, sh*Nfat5*

  
  
- Output files:
  
    `DEG_SA_KD_vs_NC_Scr_full.csv`: DEG list of SA_KD group vs NC_Scr group
  
    `DEG_SA_Scr_vs_NC_Scr_full.csv`: DEG list of SA_Scr group vs NC_Scr group

    `NFAT5_Dependency_of_DEGs_induced_by_Neu5Ac.tsv`: Dependency analysis results



  
