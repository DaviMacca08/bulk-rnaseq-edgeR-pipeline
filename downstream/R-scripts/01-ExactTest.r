library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)



# ------------------------------------------------
#         Load counts matrix and build DGE list
# ------------------------------------------------ 

# 01-Load counts matrix

count_file <- "downstream/counts.txt"
if(!file.exists(count_file)) {
  stop("counts.txt not found in current directory")
}

count <- read.table(count_file, sep = "\t", header = TRUE, skip = 1)
count <- count[, -c(2:6)]
rownames(count) <- count$Geneid
count <- count[, -1]
sample <- c("CT3", "CT2", "CT1", "CC3", "CC2", "CC1")
colnames(count) <- sample
dim(count)

# 02-Get metadata

meta <- tibble(Geneid = rownames(count))

meta$symbol <- mapIds(org.Hs.eg.db,
                     keys = meta$Geneid,
                     keytype = "ENSEMBL",
                     column = "SYMBOL",
                     multiVals = "first")
meta$entrez <- mapIds(org.Hs.eg.db,
                     keys = meta$Geneid,
                     keytype = "ENSEMBL",
                     column = "ENTREZID",
                     multiVals = "first")
meta <- as.data.frame(meta)

# 03-Create DGEList 

group <- factor(c(rep("NoTreatment", 3), rep("Treatment", 3)))
dgelist <- DGEList(counts = count, group = group, genes = meta)

# 04-Filtering using filterByExpr() 

keep <- filterByExpr(dgelist)
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

background <- data.frame(GeneId = rownames(dgelist$genes))



# ------------------------------------------------
#         Quality Control (QC)
# ------------------------------------------------ 

# 05-Plot library sizes and samples before normalization 

png("01-QC_plots.png", width = 1800, height = 600, res = 150)
par(mfrow = c(1, 3))
barplot(dgelist$samples$lib.size * 1e-6, 
        names.arg = sample, 
        ylab = "Library size before normalization (millions)",
        las = 2,
        main = "Library sizes")
plotMDS(dgelist, col = as.numeric(group), pch = 16, main = "MDS Plot")
legend("topleft", legend = levels(group), col = 1:2, pch = 16)


# 06-Boxplot before normalization

boxplot(cpm(dgelist, log = TRUE, prior.count = 1),
        las = 2,
        main = "Log2 CPM before TMM normalization",
        names = sample,
        ylab= "log2(CPM)")
par(mfrow = c(1, 1))
dev.off()

# 07-Normalization

dgelist <- calcNormFactors(dgelist, method = "TMM")

# 08-Boxplot after normalization

png("02_boxplot_log2CPM_after_TMM.png", width = 1400, height = 800, res = 150)
par(mfrow = c(1, 1))
boxplot(cpm(dgelist, log = TRUE, prior.count = 1),
        las = 2,
        main = "Log2 CPM after TMM normalization",
        names = sample,
        ylab = "log2(CPM)")
par(mfrow = c(1,1))
dev.off()

# ------------------------------------------------
#         Exploratory Data Analysis (EDA)
# ------------------------------------------------ 

# 09-MDS plots per sample 

png("03-MDplots_samples.png", width = 1200, height = 800, res = 150)
par(mfrow = c(2, 3))
for (i in 1:6) {
  plotMD(dgelist, column = i, main = colnames(dgelist)[i])
  abline(h = 0, col = "red", lty = 2, lwd = 2)
}
mtext("Mean-Difference (MD) plots for all samples", 
      outer = TRUE, cex = 1, line = -1.5)
par(mfrow = c(1, 1))
dev.off()

# 10-PCA

logCPM <- cpm(dgelist, log = TRUE, prior.count = 1)
pca <- prcomp(t(logCPM), scale. = TRUE)
pca_var <- pca$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

metadata <- data.frame(
  Sample = colnames(logCPM),
  Group = group
)

pca_df <- data.frame(pca$x[, 1:2], metadata)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 3) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  theme_bw() +
  scale_color_manual(values = c("darkgreen", "darkblue")) +
  ggtitle("PCA of normalized log2 CPM values")
ggsave("04-PCA_log2CPM.png", width = 8, height = 6, dpi = 150)

# 11-Sample distance heatmap

png("05-sample_distance_heatmap.png", width = 1200, height = 1000, res = 150)
distance <- dist(t(logCPM))
pheatmap(as.matrix(distance),
         main = "Heatmap of sample-to-sample distances",
         fontsize_row = 10,
         fontsize_col = 10)
dev.off()


# ------------------------------------------------
#                Model Diagnostics
# ------------------------------------------------ 

# 12-Dispersion estimation

dgelist <- estimateCommonDisp(dgelist)
dgelist <- estimateTagwiseDisp(dgelist)
cat("Common BCV:", sqrt(dgelist$common.dispersion), "\n")

png("06-BCV_and_meanVariance_plots.png", width = 1200, height = 600, res = 150)
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))
plotBCV(dgelist)
plotMeanVar(dgelist)
mtext("Biological Coefficient of Variation and Mean-Variance relationship plots", 
      outer = TRUE, cex = 1, line = -1.5)
par(mfrow = c(1, 1))
dev.off()

# 13-Perform exact test

exact_test <- exactTest(dgelist, pair = c("NoTreatment", "Treatment"))
res <- topTags(exact_test, n = Inf)



# ------------------------------------------------
#             Differential expression
# ------------------------------------------------

# 14-MD plot

png("07-MD_plot_exactTest.png", width = 1000, height = 800, res = 150)
plotMD(exact_test, main = "MD plot of exact-test Data")
abline(h = c(-1, 1), col = "blue", lty = 2)
dev.off()

# 15-Find differential expressed genes (DEG)

deg <- res$table
deg_fil <- deg[deg$FDR < 0.05 & abs(deg$logFC) > 1, ]
deg_up <- deg_fil[deg_fil$logFC > 0, ]
deg_down <- deg_fil[deg_fil$logFC < 0, ]

png("08-PValue_distribution.png", width = 1000, height = 800, res = 150)
hist(deg$PValue, breaks = 25,
     main = "P-value distribution",
     xlab = "P-value",
     col = "lightblue", border = "white"
     )
dev.off()

cat("\n------Differentially expressed genes------\n",
    "Total DEGs: ", nrow(deg_fil), "\n",
    "DEGs up-regulated: ", nrow(deg_up), "\n",
    "DEGs down-regulated: ", nrow(deg_down), "\n\n")

# 16-Volcano Plot

deg$category <- "Not Significant"

deg$category[deg$FDR < 0.05 & abs(deg$logFC) <= 1] <- "FDR < 0.05 only"
deg$category[deg$FDR >= 0.05 & abs(deg$logFC) > 1] <- "|log2FC| > 1 only"
deg$category[deg$FDR < 0.05 & abs(deg$logFC) > 1] <- "Significant (FDR < 0.05 & |logFC| > 1)"


cat("Category distribution:\n")
print(table(deg$category))

VolcanoPlot <- ggplot(deg, aes(x = logFC, y = -log10(FDR), color = category)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c(
      "Significant (FDR < 0.05 & |logFC| > 1)" = "red",
      "FDR < 0.05 only" = "blue",
      "|log2FC| > 1 only" = "green",
      "Not Significant" = "gray70"
    )
  ) + 
  theme_light() +
  labs(
    title = "Volcano plot of DEGs",
    x = "Log2 Fold Change",
    y = "-Log10(FDR)",
    color = "Gene Category"
  )
ggsave("09-volcano_plot_DEGs.png", plot = VolcanoPlot, width = 8, height = 6, dpi = 150)

# 17-HeatMap of DEGs

deg_genes <- rownames(deg_fil)

logCPM_deg <- logCPM[deg_genes, , drop = FALSE]  
logCPM_deg_scaled <- t(scale(t(logCPM_deg)))


if(nrow(logCPM_deg_scaled) > 0) {
  
  png("10-DEGs_heatmap.png", width = 1200, height = 1000, res = 150)
  
  annotation_col <- data.frame(
    Group = group
  )
  rownames(annotation_col) <- colnames(logCPM_deg)
  
  pheatmap(
    logCPM_deg_scaled,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale = "none",
    main = paste0("DEGs Heatmap (", nrow(logCPM_deg_scaled), " genes)"),
    color = colorRampPalette(c("blue", "white", "darkorange"))(50)
  )
  dev.off()
} else {
  warning("No DEGs to plot in heatmap")
}

# 18-Save output for further analysis

output_dir <- "downstream/objects"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}


saveRDS(logCPM, file = file.path(output_dir, "logCPM_counts.rds"))
saveRDS(deg_fil, file = file.path(output_dir, "ExactTestDEG.rds"))
saveRDS(background, file = file.path(output_dir, "Background_genes.rds"))
saveRDS(res$table, file = file.path(output_dir, "ResFromExactTest.rds"))

cat("\nAll files saved successfully to:", output_dir, "\n")

# 19-Save session info

sink("ExactTest_info.txt")
cat("Exact Test Analysis Session Info\n")
cat("===============================\n\n")
sessionInfo()
sink()

cat("\nExact Test analysis completed successfully!\n")
cat("Results saved to current working directory\n")
