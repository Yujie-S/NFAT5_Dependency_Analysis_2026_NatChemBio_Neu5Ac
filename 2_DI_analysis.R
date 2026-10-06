
library(readr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(tibble)

# ------------------ 1. Load the data ------------------
df_KD <- read_csv("DEG_SA_KD_vs_NC_Scr_full.csv", show_col_types = FALSE)
df_Scr <- read_csv("DEG_SA_Scr_vs_NC_Scr_full.csv", show_col_types = FALSE)
df_KD$gene_label <- paste0(df_KD$gene_name, "_", df_KD$gene_id)
df_Scr$gene_label <- paste0(df_Scr$gene_name, "_", df_Scr$gene_id)

colnames(df_KD) <- trimws(colnames(df_KD))
colnames(df_Scr) <- trimws(colnames(df_Scr))

df_KD  <- df_KD  %>% select(gene_label, log2FoldChange, padj) %>% as.data.frame()
df_Scr <- df_Scr %>% select(gene_label, log2FoldChange, padj) %>% as.data.frame()

df_KD <- df_KD %>% dplyr::rename('lfc_KD' = 'log2FoldChange',
                                 'padj_KD' = 'padj')
df_Scr <- df_Scr %>% dplyr::rename('lfc_Scr' = 'log2FoldChange',
                                   'padj_Scr' = 'padj')

common <- df_KD%>%
  select(gene_label, lfc_KD, padj_KD) %>%
  inner_join(
    df_Scr %>% select(gene_label, lfc_Scr, padj_Scr),
    by = c("gene_label")
  )

# ------------------ 2. Dependency Index Calculation & Categorize ------------------
# Filter
sig_in_scr <- common %>%
  filter(padj_Scr < 0.05 & abs(lfc_Scr) > 0.25)
cat("Number of genes significantly changed by SA in Scrambled control:", nrow(sig_in_scr), "\n\n")

# analyze dependency
# Add dependency metrics

sig_in_scr <- sig_in_scr %>%
  mutate(
    same_sign = sign(lfc_KD) == sign(lfc_Scr) | (lfc_KD == 0 & lfc_Scr == 0),
    ratio     = ifelse(lfc_Scr == 0, NA, lfc_KD / lfc_Scr),
    abs_ratio = abs(ratio),
    delta     = lfc_Scr - lfc_KD,
    dependency = case_when(
      !same_sign & padj_KD < 0.05         ~ "NFAT5-independent (Opposite direction)",
      padj_KD >= 0.05                     ~ "Strongly NFAT5-dependent",
      abs_ratio < 0.68                    ~ "Strongly NFAT5-dependent",
      abs_ratio >= 1                      ~ "NFAT5-independent",
      TRUE                                ~ "Partially dependent"
    )
  )

sig_in_scr <- sig_in_scr %>% separate(col = gene_label, into = c("gene_name", "gene_id"), sep = "_")

# Quick summary
cat("Dependency classification:\n")
print(table(sig_in_scr$dependency))

# ------------------ 3. Scatter Plot ------------------

highlight_genes <- c(
  # NFAT5-dependent list
  "Id2", "Lef1", "Bcl6", "Tcf7", "Satb1", "Tbx21", "Hif1a", "Nfatc1", "Nfkb2", 
  "Itgae", "Foxm1", "Serpina3g", "Ccr7", "Myb", "Nt5e", "Bcl2l11", 
  "Pdcd1", "Havcr2", "Tnfsf10", "Il10", "Fasl", "Ifng", "Prf1", "Xcl1", 
  "Ccl4", "Ccl3", "Cd28", "Icos", "Tnfrsf9", "Il2rb",
  
  # NFAT5-independent list
  "Batf", "Runx2", "Tox", "Eomes", "Bach2", "Il2rg", 
  "Fas", "Sell", "Ly6a", "Klf2", "Slamf6", "Igf1r", "Foxp1", 
  "Cd38", "Lag3", "Cd101", "Entpd1", "Tigit", "Gzma", "Tnf", "Cxcl10", 
  "Ccl5", "Cxcr5", "Tnfrsf4", "Tnfsf14", "Klrg1", "Il2"
)

ggplot(sig_in_scr, aes(x = lfc_Scr, y = lfc_KD, color = dependency)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey80") +
  geom_point(alpha = 0.3, size = 2) +
  geom_abline(slope = 1, intercept = 0, color = "#75bd31", linetype = "dashed", linewidth = 0.8) +
  geom_abline(slope = 0.68, intercept = 0, color = "#1c8538", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = c(
    "Strongly NFAT5-dependent"               = "#1c8538",   # dark red
    "Partially dependent"                    = "#75bd31",      # light red
    "NFAT5-independent"                      = "#e8cc3f",   # gray
    "NFAT5-independent (Opposite direction)" = "#e8cc3f"    # gray
  )) +
  geom_point(data = subset(sig_in_scr, gene_name %in% highlight_genes),
             size = 2, shape = 21, color = "gray40", stroke = 1.2) +
  geom_text_repel(data = subset(sig_in_scr, gene_name %in% highlight_genes),
                  aes(label = gene_name),
                  size = 5, fontface = "italic", max.overlaps = Inf, 
                  box.padding = 2, segment.color = "gray40",
                  segment.size = 0.5) +
  labs(title = "SA Effect Dependency on NFAT5\n(Only genes significant in Neu_Scr vs NC_Scr)",
       subtitle = paste("n =", nrow(sig_in_scr), "genes with padj_Scr < 0.05 & |log2FC| > 0.25"),
       x = "log2FC (Neu_Scr vs Un_Scr)",
       y = "log2FC (Neu_KD vs Un_Scr)",
       color = "Dependency on NFAT5") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right",
        panel.grid.minor = element_blank())

# ------------------ 4. Save DI analysis results ------------------

write_tsv(sig_in_scr, "NFAT5_Dependency_of_DEGs_induced_by_Neu5Ac.tsv")


# ------------------ 5. Heatmap of key T cell genes ------------------

df <- read_tsv("gene_count.xls", show_col_types = FALSE)

# Convert to regular data.frame and set gene_name as row names (clean way)
count_df <- as.data.frame(df)
count_df$gene_label <- paste0(df$gene_name, "_", df$gene_id)
rownames(count_df) <- count_df$gene_label

# Extract count matrix (only the sample columns)
sample_cols <- c("NC_Scr_1", "NC_Scr_2", "NC_Scr_3",
                 "NC_KD_1", "NC_KD_2", "NC_KD_3",
                 "SA_Scr_1", "SA_Scr_2", "SA_Scr_3",
                 "SA_KD_1", "SA_KD_2", "SA_KD_3")

count_matrix <- as.matrix(count_df[, sample_cols])

# Create colData for DESeq2
colData <- data.frame(
  Group = factor(c(rep("NC_Scr",3), rep("NC_KD",3), 
                   rep("SA_Scr",3), rep("SA_KD",3)),
                 levels = c("NC_Scr", "NC_KD", "SA_Scr", "SA_KD")),
  row.names = colnames(count_matrix)
)

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),   # must be integers
  colData = colData,
  design = ~ Group
)

# Run DESeq2 normalization (size factor estimation)
dds <- estimateSizeFactors(dds)

# Get normalized counts
norm_counts <- as.data.frame(counts(dds, normalized = TRUE))
group_means <- data.frame(
  NC_Scr = rowMeans(norm_counts[,1:3]),
  NC_KD  = rowMeans(norm_counts[,4:6]),
  SA_Scr = rowMeans(norm_counts[,7:9]),
  SA_KD  = rowMeans(norm_counts[,10:12])
)

norm_counts <- rownames_to_column(norm_counts, 'gene_label')
norm_counts <- norm_counts %>% 
  separate(col = gene_label, into = c("gene_name", "gene_id"), sep = "_")

group_means <- rownames_to_column(group_means, 'gene_label')
group_means <- group_means %>%
  separate(col = gene_label, into = c("gene_name", "gene_id"), sep = "_")

genes <- highlight_genes

# Keep only genes present in the data
present_genes <- genes[genes %in% sig_in_scr$gene_name]
missing <- setdiff(genes, present_genes)
sub_mat_grp_means <- group_means[group_means$gene_name %in% present_genes, , drop = FALSE]
sub_mat_norm_counts <- norm_counts[norm_counts$gene_name %in% present_genes, , drop = FALSE]

# numeric values only
sub_mat_grp_means_val <- sub_mat_grp_means %>% select(-gene_id, -NC_KD)
rownames(sub_mat_grp_means_val) <- sub_mat_grp_means_val$gene_name
sub_mat_grp_means_val$gene_name <- NULL

sub_mat_norm_counts_val <- sub_mat_norm_counts %>% select(-gene_id, -NC_KD_1, -NC_KD_2, -NC_KD_3)
rownames(sub_mat_norm_counts_val) <- sub_mat_norm_counts_val$gene_name
sub_mat_norm_counts_val$gene_name <- NULL


# Log2 + pseudo-count
log_mat_grp <- log2(sub_mat_grp_means_val + 1)
log_mat_sample <- log2(sub_mat_norm_counts_val + 1)

# Row-wise z-score
z_mat_grp <- t(scale(t(log_mat_grp)))
z_mat_sample <- t(scale(t(log_mat_sample)))

# z_mat is row-wise z-score matrix for selected T cell genes
# Columns: NC_Scr, SA_Scr, SA_KD   (or group averages)
# Rows: gene names (use gene_label = GeneName_GeneID)

# Create annotation for dependency (use group means)
row_ann_grp <- data.frame(
  Dependency = sig_in_scr$dependency[match(rownames(z_mat_grp), sig_in_scr$gene_name)],
  row.names = rownames(z_mat_grp)
)
row_ann_sample <- data.frame(
  Dependency = sig_in_scr$dependency[match(rownames(z_mat_sample), sig_in_scr$gene_name)],
  row.names = rownames(z_mat_sample)
)

# Colors for annotation
ann_colors_grp <- list(
  Dependency = c(
    "Strongly NFAT5-dependent" = "#1c8538",   # dark red
    "Partially dependent"      = "#75bd31",      # light red
    "NFAT5-independent"        = "#e8cc3f",   # gray
    "NFAT5-independent (Opposite direction)"       = "#e8cc3f"    # gray
    )
)

# Heatmap
pheatmap(
  z_mat_grp,
  annotation_row = row_ann_grp,
  annotation_colors = ann_colors_grp,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 9,
  fontsize_col = 11,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks = seq(-2, 2, length.out = 101),
  main = "T cell fate genes: Neu5Ac effect and NFAT5 dependency")

# ------------------ 6. Heatmap re-cluster and biological annotation ------------------
# Define gene categories (same as bar plot)
gene_categories <- tribble(
  ~gene_name,          ~Biological_Category,
  "Tcf7",              "Progenitor / Memory-like",
  "Lef1",              "Progenitor / Memory-like",
  "Bcl6",              "Progenitor / Memory-like",
  "Satb1",             "Progenitor / Memory-like",
  "Sell",              "Progenitor / Memory-like",
  "Ccr7",              "Progenitor / Memory-like",
  "Il7r",              "Progenitor / Memory-like",
  "Slamf6",            "Progenitor / Memory-like",
  "Foxp1",             "Progenitor / Memory-like",
  "Klf2",              "Progenitor / Memory-like",
  "Bach2",             "Progenitor / Memory-like",
  "Id3",               "Progenitor / Memory-like",
  "Itgae",             "Progenitor / Memory-like",
  "Bcl2l11",           "Progenitor / Memory-like",
  
  "Tox",               "Terminal Exhaustion",
  "Prdm1",             "Terminal Exhaustion",
  "Eomes",             "Terminal Exhaustion",
  "Id2",               "Terminal Exhaustion",
  "Ccl4",              "Terminal Exhaustion",
  "Ccl3",              "Terminal Exhaustion",
  
  "Ifng",              "Effector / Cytotoxic",
  "Prf1",              "Effector / Cytotoxic",
  "Gzma",              "Effector / Cytotoxic",
  "Tnf",               "Effector / Cytotoxic",
  "Fasl",              "Effector / Cytotoxic",
  "Tnfsf10",           "Effector / Cytotoxic",
  "Batf",              "Effector / Cytotoxic",
  "Xcl1",              "Effector / Cytotoxic",
  "Tnfrsf4",           "Effector / Cytotoxic",
  
  "Pdcd1",             "Inhibitory Receptors",
  "Havcr2",            "Inhibitory Receptors",
  "Lag3",              "Inhibitory Receptors",
  "Tigit",             "Inhibitory Receptors",
  "Entpd1",            "Inhibitory Receptors",
  "Cd38",              "Inhibitory Receptors",
  "Ctla4",             "Inhibitory Receptors",
  "Klrg1",             "Inhibitory Receptors",
  
  "Cd28",              "Co-stimulatory / Others",
  "Icos",              "Co-stimulatory / Others",
  "Tnfrsf9",           "Co-stimulatory / Others",
  "Il2rb",             "Co-stimulatory / Others",
  "Il2rg",             "Co-stimulatory / Others",
  
  
  
)

row_ann_grp_bio <- data.frame(
  Biological_Category = gene_categories$Biological_Category[match(rownames(z_mat_grp), gene_categories$gene_name)],
  Dependency = sig_in_scr$dependency[match(rownames(z_mat_grp), sig_in_scr$gene_name)],
  row.names = rownames(z_mat_grp)
)

# Order rows by Biological Category (this enables clustering within groups)
row_ann_grp_bio <- row_ann_grp_bio %>%
  mutate(Biological_Category = factor(Biological_Category, 
                                      levels = c("Progenitor / Memory-like", 
                                                 "Terminal Exhaustion", 
                                                 "Effector / Cytotoxic", 
                                                 "Inhibitory Receptors", 
                                                 "Co-stimulatory / Others"))) %>%
  arrange(Biological_Category)


# Reorder z_mat to match the annotation order
z_mat_grp_ordered <- z_mat_grp[rownames(row_ann_grp_bio), , drop = FALSE]

ann_colors_bio <- list(
  Biological_Category = c(
    "Progenitor / Memory-like" = "#1f77b4",
    "Terminal Exhaustion"      = "#d62728",
    "Effector / Cytotoxic"     = "#ff7f0e",
    "Inhibitory Receptors"     = "#9467bd",
    "Co-stimulatory / Others"  = "#8c564b"
  ),
  Dependency = c(
    "Strongly NFAT5-dependent" = "#1c8538",   # dark red
    "Partially dependent"      = "#75bd31",      # light red
    "NFAT5-independent"        = "#e8cc3f",   # gray
    "NFAT5-independent (Opposite direction)"       = "#e8cc3f"    # gray
  )
)

pheatmap(z_mat_grp_ordered,
         annotation_row = row_ann_grp_bio,
         annotation_colors = ann_colors_bio,
         cluster_rows = TRUE,      # This will now cluster within categories nicely
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 9.5,
         fontsize_col = 11,
         color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         breaks = seq(-3, 3, length.out = 101),
         main = "T cell fate genes: Neu5Ac effect and NFAT5 dependency",
         border_color = NA,
         cutree_rows = 5)   # optional: cut into 5 clusters if desired

