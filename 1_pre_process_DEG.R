
library(readr)
library(dplyr)
library(DESeq2)
library(tibble)


# ------------------ 1. Load data ------------------
cat("Loading gene_count.xls...\n")
df <- read_tsv("gene_count.xls", show_col_types = FALSE)

count_df <- as.data.frame(df)
count_df$gene_label <- paste0(df$gene_name, "_", df$gene_id)
rownames(count_df) <- count_df$gene_label
df2 <- as.data.frame(df)
df2$gene_label <- paste0(df$gene_name, "_", df$gene_id)

count_matrix <- as.matrix(count_df[, c("NC_Scr_1","NC_Scr_2","NC_Scr_3",
                                       "NC_KD_1","NC_KD_2","NC_KD_3",
                                       "SA_Scr_1","SA_Scr_2","SA_Scr_3",
                                       "SA_KD_1","SA_KD_2","SA_KD_3")])

# Filter very low count genes
count_matrix <- count_matrix[rowSums(count_matrix) >= 10, ]

# ------------------ 2. Metadata ------------------
colData <- data.frame(
  Group = factor(c(rep("NC_Scr",3), rep("NC_KD",3), 
                   rep("SA_Scr",3), rep("SA_KD",3)),
                 levels = c("NC_Scr", "NC_KD", "SA_Scr", "SA_KD")),
  row.names = colnames(count_matrix)
)

# ------------------ 3. Run DESeq2 ------------------
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = colData,
                              design = ~ Group)

dds <- DESeq(dds)

# ------------------ 4. Results SA_KD vs NC_Scr ------------------
cat("\nExtracting results for SA_KD vs NC_Scr...\n")

# Get results without shrinkage first
res <- results(dds, 
               contrast = c("Group", "SA_KD", "NC_Scr"),
               alpha = 0.05)

res_2 <- results(dds,
                 contrast = c("Group", "SA_Scr", "NC_Scr"),
                 alpha = 0.05)

cat("\nAvailable coefficients:\n")
resultsNames(dds)

coef_name <- "Group_SA_KD_vs_NC_Scr"
res_shrunk <- lfcShrink(dds, 
                        coef = coef_name,      
                        type = "apeglm")


coef_name_2 <- "Group_SA_Scr_vs_NC_Scr"
res_shrunk_2 <- lfcShrink(dds,
                          coef = coef_name_2,
                          type = "apeglm")


# ------------------ 5. Add annotations and save ------------------
res_df <- as.data.frame(res_shrunk)
gene_label <- rownames_to_column(res_df, 'gene_label')
res_df2 <- left_join(gene_label, res_df)
annot <- select(df2, c(gene_label, gene_name, gene_id, gene_description))
res_df3 <- left_join(res_df2, annot, by = 'gene_label')

res_2_df <- as.data.frame(res_shrunk_2)
gene_label_2 <- rownames_to_column(res_2_df, 'gene_label')
res_2_df2 <- left_join(gene_label_2, res_2_df)
annot_2 <- select(df2, c(gene_label, gene_name, gene_id, gene_description))
res_2_df3 <- left_join(res_2_df2, annot_2, by = 'gene_label')

# Significant DEGs
sig_df <- res_df3 %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(desc(abs(log2FoldChange)))

cat(sprintf("\nSignificant DEGs (padj < 0.05 and |log2FC| > 1): %d genes\n", nrow(sig_df)))

# Save
write.csv(res_df3,  "DEG_SA_KD_vs_NC_Scr_full.csv", row.names = FALSE)
write.csv(res_2_df3,  "DEG_SA_Scr_vs_NC_Scr_full.csv", row.names = FALSE)


cat("\n✅ Files saved successfully!\n")
cat("   • DEG_SA_KD_vs_NC_Scr_full.csv\n")
cat("   • DEG_SA_KD_vs_NC_Scr_significant.csv\n")

