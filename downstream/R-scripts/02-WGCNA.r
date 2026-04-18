library(WGCNA)
library(ggrepel)
library(ggplot2)
library(patchwork)
library(tidyverse)



# ------------------------------------------------
#               Quality check (QC)
# ------------------------------------------------

# 01-Load count matrix and organize metadata

input_file <- "downstream/objects/logCPM_counts.rds"
if (!file.exists(input_file)) {
  stop("logCPM_counts.rds not found at: ", input_file, "\n",
       "Please run 01-ExactTest.r first to generate this file.")
}

logCPM <- readRDS(input_file)

traitData <- data.frame(
  Cocultured = c(0, 0, 0, 1, 1, 1),
  Control = c(1, 1, 1, 0, 0, 0)
)
rownames(traitData) <- colnames(logCPM)

# 02-Detecting outliers

count_wgcna <- t(logCPM)


gsg <- goodSamplesGenes(count_wgcna)
if (!gsg$allOK) {
  cat("Removing bad samples/genes...\n")
  count_wgcna <- count_wgcna[gsg$goodSamples, gsg$goodGenes]

  traitData <- traitData[rownames(count_wgcna), , drop = FALSE]
  cat("Removed", sum(!gsg$goodSamples), "samples and", sum(!gsg$goodGenes), "genes\n")
} else{
  print("Samples and Genes are OK")
}

# Sample clustering

png("01-sampleClustering.png", width = 1200, height = 800, res = 150)
htree <- hclust(dist(count_wgcna), method = "average")
plot(htree, 
     main = "Sample clustering dendrogram (WGCNA QC)", 
     xlab = "", 
     sub = "",
     cex = 1.2)
dev.off()

# 03-PCA

pca_wgcna <- prcomp(count_wgcna, scale. = TRUE)
pca_wgcna_data <- as.data.frame(pca_wgcna$x)
pca_wgcna_var <- pca_wgcna$sdev^2
pca_wgcna_perce <- round(pca_wgcna_var / sum(pca_wgcna_var) * 100, 2)

pca_wgcna_data$Treatment <- factor(traitData$Cocultured, 
                                   labels = c("Control", "Cocultured"))
pca_wgcna_data$Sample <- rownames(pca_wgcna_data)

ggplot(pca_wgcna_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 4, alpha = 0.9) +
  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = Inf) +
  stat_ellipse(aes(group = Treatment), linetype = 2, linewidth = 0.8) +
  scale_color_manual(values = c("darkgreen", "darkblue")) +
  labs(
    title = "PCA of WGCNA input data",
    x = paste0("PC1 (", pca_wgcna_perce[1], "% variance)"),
    y = paste0("PC2 (", pca_wgcna_perce[2], "% variance)"),
    color = "Treatment"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("02-PCA_samples.png", width = 8, height = 6, dpi = 150)


# ------------------------------------------------
#                Soft-thresholding
# ------------------------------------------------

# 04-Find the correct power

power <- c(1:10, seq(12, 30, by = 2))

sft <- pickSoftThreshold(count_wgcna,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 0)

sft_data <- sft$fitIndices

# Identify optimal power programmatically

r_sq_cutoff <- 0.8
optimal_power <- sft_data$Power[which(sft_data$SFT.R.sq > r_sq_cutoff)[1]]

if(is.na(optimal_power)) {
  optimal_power <- 28  # default fallback
  warning("No power achieved R² > ", r_sq_cutoff, ". Using default power = ", optimal_power)
} else {
  cat("Optimal power (R² >", r_sq_cutoff, "):", optimal_power, "\n")
}

power_plot <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) + 
  geom_point(color = "red", size = 2) + 
  geom_text(nudge_y = 0.03, size = 3) + 
  geom_hline(yintercept = r_sq_cutoff, color = "blue", linetype = "dashed") + 
  labs(
    x = "Soft Threshold (power)",
    y = "Scale Free Topology Model Fit (R²)"
  ) + 
  theme_bw() +
  ylim(0, 1)

connectivity_plot <- ggplot(sft_data, aes(Power, mean.k., label = Power)) +
  geom_point(color = "red", size = 2) +
  geom_text(nudge_y = 0.5, size = 3) +
  labs(
    x = "Soft Threshold (power)",
    y = "Mean Connectivity"
  ) +
  theme_bw()

png("03-SoftThreshold.png", width = 1200, height = 600, res = 150)
power_plot + connectivity_plot + 
  plot_annotation(title = "Soft-thresholding analysis")
dev.off()


# ------------------------------------------------
#             Clustering genes on modules
# ------------------------------------------------

# 05-Use optimal power

soft_power <- optimal_power

cat("Calculating TOM similarity (this may take a while)...\n")
tom <- TOMsimilarityFromExpr(count_wgcna, 
                             power = soft_power, 
                             TOMType = "signed", 
                             nThreads = 4)
dissTOM <- 1 - tom

# 06-Blockwise modules

bwnet <- blockwiseModules(count_wgcna,
                          maxBlockSize = 5000,
                          minModuleSize = 30,
                          TOMType = "signed",
                          corType = "pearson",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 0)

# 07-Plot dendrogram with colors

geneTree <- hclust(as.dist(dissTOM), method = "average")

png("04-GeneDendro_Modules.png", width = 1400, height = 800, res = 150)
plotDendroAndColors(
  dendro = geneTree,
  colors = cbind(bwnet$unmergedColors, bwnet$colors),
  groupLabels = c("Unmerged", "Merged"),
  dendroLabels = FALSE,
  addGuide = TRUE,
  hang = 0.03,
  main = "Gene dendrogram and module colors"
)
dev.off()

# 08-Relate modules to traits

module_eigengenes <- orderMEs(bwnet$MEs)
nSamples <- nrow(count_wgcna)

module_trait_corr <- cor(module_eigengenes, traitData, use = "p")
module_trait_corr_pvalue <- corPvalueStudent(module_trait_corr, nSamples)

# Ensure dimensions match

if(nrow(module_trait_corr_pvalue) == nrow(module_trait_corr) && 
   ncol(module_trait_corr_pvalue) == ncol(module_trait_corr)) {
  
  # -Create text matrix with significance stars
  
  sigSymbols <- ifelse(module_trait_corr_pvalue < 0.001, "***",
                       ifelse(module_trait_corr_pvalue < 0.01, "**",
                              ifelse(module_trait_corr_pvalue < 0.05, "*", "")))
  
  textMatrix <- paste0(sprintf("%.2f", module_trait_corr), sigSymbols)
  dim(textMatrix) <- dim(module_trait_corr)
  
  # -Heatmap
  
  png("05-Module_trait_heatmap.png", width = 2400, height = 3200, res = 300)
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(
    Matrix = module_trait_corr,
    xLabels = colnames(traitData),
    yLabels = names(module_eigengenes),
    ySymbols = colnames(module_eigengenes),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.8,
    zlim = c(-1, 1),
    main = "Module–trait relationships"
  )
  dev.off()
} else {
  stop("Dimension mismatch in correlation matrices")
}



# -------------------------------------------------------
#   Extract genes from best ME and intersect with DEGs
# -------------------------------------------------------

# 09-Find module with highest correlation to Cocultured

cocultured_cor <- module_trait_corr[, "Cocultured", drop = FALSE]
top_module <- rownames(cocultured_cor)[which.max(abs(cocultured_cor))]
top_module_name <- gsub("^ME", "", top_module)  

cat("\nTop module correlated with Cocultured:", top_module_name, 
    "(r =", round(cocultured_cor[top_module, ], 3), ")\n")

module_of_interest <- top_module_name

# Check if module exists

moduleGeneMapping <- data.frame(
  Gene = colnames(count_wgcna),
  Module = bwnet$colors,
  stringsAsFactors = FALSE
)

if(module_of_interest %in% unique(moduleGeneMapping$Module)) {
  
  gene_module <- moduleGeneMapping |>
    filter(Module == module_of_interest) |>
    pull(Gene)
  
  cat("Genes in", module_of_interest, "module:", length(gene_module), "\n")
  
  if(length(gene_module) > 0) {
    ME_module <- bwnet$MEs[, paste0("ME", module_of_interest)]
    
    MM <- cor(count_wgcna[, gene_module], ME_module, method = "pearson")
    GS <- cor(count_wgcna[, gene_module], traitData$Cocultured, method = "pearson")
    
    module_MM_GS <- data.frame(
      Gene = gene_module,
      MM = as.numeric(MM),
      GS = as.numeric(GS)
    )
    
    # -Plot MM vs GS
    
    p_mm_gs <- ggplot(module_MM_GS, aes(x = MM, y = GS, label = Gene)) +
      geom_point(shape = 1, color = "lightcoral", size = 2, stroke = 1) +
      geom_vline(xintercept = 0.9, linetype = "dashed", color = "gray") +
      geom_hline(yintercept = 0.9, linetype = "dashed", color = "gray") +
      xlim(0, 1) + ylim(0, 1) +
      labs(
        x = paste("Module Membership (", module_of_interest, " module)"),
        y = "Gene Significance for Cocultured",
        title = paste("MM vs GS for", module_of_interest, "module")
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    print(p_mm_gs)
    ggsave(paste0("06-MM_vs_GS_", module_of_interest, ".png"), 
           p_mm_gs, width = 8, height = 6, dpi = 300)
    
    # -Select top genes
    
    topGenesModule <- module_MM_GS |>
      filter(MM > 0.9 & GS > 0.9) |>
      pull(Gene)
    
    cat("\nGenes with MM and GS > 0.9:", length(topGenesModule), "\n")
    
    topGenesModule <- data.frame(GeneId = topGenesModule)
    
    # -Intersect with DEGs
    
    deg_file <- "downstream/objects/ExactTestDEG.rds"
    if (file.exists(deg_file)) {
      exactTestDEG <- readRDS(deg_file)
      intersection <- intersect(topGenesModule$GeneId, exactTestDEG$Geneid)
      cat("Intersection with DEGs:", length(intersection), "genes\n")
    }
    
    # -Save results
    output_dir <- "downstream/objects"
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    saveRDS(topGenesModule, 
            file = file.path(output_dir, paste0("genes_ME_", module_of_interest, ".rds")))
    
    saveRDS(intersection,
            file = file.path(output_dir, "CommonGenes.rds"))
  }
} else {
  warning("Module ", module_of_interest, " not found in results")
}

# 10-Save session info

sink("WGCNA_info.txt")
cat("WGCNA Analysis Session Info\n")
cat("===============================\n\n")
sessionInfo()
sink()

cat("\nWGCNA analysis completed successfully!\n")
cat("Results saved to current working directory\n")
