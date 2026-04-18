library(tidyverse)
library(pheatmap)
library(readxl)
library(ggrepel)
library(hgu133a.db)
library(survival)



# -------------------------------------------------------
#               Load GSE2034 GEO Dataset
# -------------------------------------------------------

# Prepare function

run_cox <- function(expr_mat, meta_df) {
  results <- data.frame(Gene = colnames(expr_mat),
                        HR = NA, lower_CI = NA, upper_CI = NA, pvalue = NA)
  
  for (i in 1:ncol(expr_mat)) {
    gene_expr <- expr_mat[, i]
    if (all(is.na(gene_expr)) || sd(gene_expr) == 0) next
    
    cox_model <- tryCatch(
      coxph(Surv(meta_df$time, meta_df$relapse) ~ gene_expr),
      error = function(e) NULL
    )
    if (is.null(cox_model)) next
    
    cox_summary <- summary(cox_model)
    results$HR[i] <- cox_summary$coef[1, "exp(coef)"]
    results$lower_CI[i] <- cox_summary$conf.int[,"lower .95"][1]
    results$upper_CI[i] <- cox_summary$conf.int[,"upper .95"][1]
    results$pvalue[i] <- cox_summary$coef[1, "Pr(>|z|)"]
  }
  
  return(results)
}


# 01-Load GSE2034 GEO Dataset

base_dir <- "downstream/Dataset_SurvivalAnalysis"

if(!dir.exists(base_dir)) cat("\nThere is not directory!")
count <- read.table(
  file.path(base_dir, "GSE2034_series_matrix.txt"),
  header = TRUE,
  sep = "\t",
  comment.char = "!",
  stringsAsFactors = FALSE
)

rownames(count) <- count[,1]
count <- count[,-1]

count_numb <- as.matrix(count)
mode(count_numb) <- "numeric"

# Check if the matrix is normnalized or not

check <- count[, c("GSM36777", "GSM36778", "GSM36779")]

png("01-CheckNormalization.png", width = 1200, height = 1100, res = 150)
boxplot(check,
        las = 2,       
        main = "Check of Expression Levels",
        ylab = "Expression Level",
        outline = FALSE,      
        col = "lightblue")

paste0("The series matrix is normalized")
dev.off()

logCount <- log2(count + 1)



# -------------------------------------------------------
#          Detecting Outliers
# -------------------------------------------------------

# 02-PCA of samples

pca <- prcomp(t(logCount), scale. = TRUE)
pca_var <- pca$sdev^2
pca_var_percent <- round(100*pca_var / sum(pca_var), 1)

## Use PCs explaining 60% variance

cum_var <- cumsum(pca_var_percent)
n_pc <- which(cum_var >= 60)[1]

pca_scores <- pca$x[,1:n_pc]

## Distance from origin in PC space

dist_center <- sqrt(rowSums(pca_scores^2))
threshold <- mean(dist_center) + 3*sd(dist_center)
outliers <- names(dist_center[dist_center > threshold])

paste0("The outliers detected are: ", paste(outliers, collapse = " - "))

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = colnames(count)
)
pca_df$outlier <- pca_df$sample %in% outliers

pca <- ggplot(pca_df, aes(PC1, PC2, color = outlier)) +
  geom_point(size = 3) +
  geom_text_repel(
    data = subset(pca_df, outlier == TRUE),
    aes(label = sample),
    size = 3
  ) +
  theme_minimal() +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  ggtitle("PCA of Samples Highlighting Outliers") 
ggsave(filename = "02-PCA_outliers.png", plot = pca, width = 8, height = 6, dpi = 150)

# 03-Clustering the samples

sample_dist <- dist(t(logCount), method = "euclidean")
hc <- hclust(sample_dist, method = "average")
plot(hc, main = "Hierarchical Clustering of Samples")

# 04-HeatMap of Samples

ann_row <- data.frame(Outlier = colnames(logCount) %in% outliers)
rownames(ann_row) <- colnames(logCount)
ann_colors <- list(Outlier = c("TRUE" = "red", "FALSE" = "white"))
rownames_to_show <- ifelse(colnames(logCount) %in% outliers, colnames(logCount), "")
ann_row$Outlier <- as.factor(ann_row$Outlier)

pheatmap(as.matrix(sample_dist),
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = FALSE,
         annotation_row = ann_row,
         annotation_col = ann_row,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         labels_row = rownames_to_show,
         show_colnames = FALSE,
         main = "Sample Distance Heatmap with Outliers Highlighted"
)

# 05-Eliminate the outliers from counts matrix and match PROBEID <-> SYMBOL

count_clean <- logCount[,!colnames(logCount) %in% outliers]
dim(count_clean)
class(count_clean)

count_clean_collapse <- count_clean
count_clean_collapse$PROBEID <- rownames(count_clean)

probe2gene <- AnnotationDbi::select(
  hgu133a.db,
  keys = rownames(count_clean),   
  columns = c("SYMBOL"),          
  keytype = "PROBEID"
)

count_clean_collapse <- left_join(count_clean_collapse, probe2gene, by = "PROBEID")

# 06-Load Metadata

meta <- read_excel(file.path(base_dir, "GSE2034_clinical.xlsx"), col_names = TRUE)
meta <- as.data.frame(meta)
meta_clean <- meta[!meta$GEO_asscession_number %in% outliers,]
rownames(meta_clean) <- meta_clean$GEO_asscession_number
meta_clean <- meta_clean[,-2]
meta_clean$time <- meta_clean$time to relapse or last follow-up (months)

status <- meta_clean$relapse
names(status) <- rownames(meta_clean)

# 07-Identified the most hub genes

hub <- read_csv(file.path(base_dir, "MCC_ValueClean.csv"), col_names = TRUE)

hub <- hub |> 
  arrange(desc(Score))

hub_matrix <- count_clean_collapse |> 
  filter(SYMBOL %in% hub$Name)

hub_matrix_ready <- hub_matrix |> 
  group_by(SYMBOL) |> 
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

hub_matrix_ready <- as.data.frame(hub_matrix_ready)
rownames(hub_matrix_ready) <- hub_matrix_ready$SYMBOL
hub_matrix_ready <- hub_matrix_ready[,-1]
hub_matrix_ready <- hub_matrix_ready[, names(status)]

## Check before run statistics

length(status)
ncol(hub_matrix_ready)
all(colnames(hub_matrix_ready) %in% names(status))

## Wilcoxon test

wilcoxon <- apply(hub_matrix_ready, 1, function(x){
  wilcox.test(x ~ status)$p.value
})

hub_wilcoxon <- data.frame(
  SYMBOL = rownames(hub_matrix_ready),
  p_value = wilcoxon,
  FDR = p.adjust(wilcoxon, method = "BH")
)

hub_wilcoxon <- hub_wilcoxon |> 
  arrange(FDR)

significant_hub <- unique(hub_wilcoxon$SYMBOL[1:9])

top_hub_matrix <- hub_matrix_ready[significant_hub,]

## boxplot

top_hub_long <- top_hub_matrix |> 
  rownames_to_column("SYMBOL") |> 
  pivot_longer(
    cols = -SYMBOL,
    names_to = "Sample",
    values_to = "Expression"
  ) |>  
  mutate(
    Relapse = factor(status[Sample]),
    ER = factor(meta_clean$ER_Status[match(Sample, rownames(meta_clean))])
  )

fdr_annot <- top_hub_long |> 
  group_by(SYMBOL, ER) |> 
  summarise(
    p = wilcox.test(Expression ~ Relapse)$p.value,
    FDR = p.adjust(p, method = "BH"),
    .groups = "drop"
  ) |> 
  mutate(label = case_when(
    FDR < 0.01 ~ "**",
    FDR < 0.05 ~ "*",
    TRUE ~ "ns"
  ),
  x_pos = as.numeric(factor(ER)))

y_max <- top_hub_long |> 
  group_by(SYMBOL, ER) |> 
  summarise(
    y_pos = max(Expression, na.rm = TRUE) * 1.05,
    .groups = "drop"
  )

fdr_annot <- left_join(fdr_annot, y_max, by = c("SYMBOL", "ER"))

boxplot_status_hub <- ggplot(top_hub_long, aes(x = ER, y = Expression, fill = Relapse)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ SYMBOL, scales = "free_y") +
  geom_text(
    data = fdr_annot,
    aes(x = x_pos, y = y_pos, label = label),
    inherit.aes = FALSE,
    fontface = "bold", na.rm = TRUE
  ) +
  theme_bw() +
  labs(
    x = "ER Status",
    y = "log2(Expression)",
    fill = "ER Status"
  ) +
  scale_fill_manual(
    values = c("0" = "#69b3a2", "1" = "#5e3c99"),
    labels = c("0" = "ER-", "1" = "ER+")
  ) + 
  theme(
    strip.text = element_text(face = "bold"),
    axis.ticks.x = element_blank()
  )
ggsave(filename = "03-HubGenes_ER_status.png", plot = boxplot_status_hub, width = 12, height = 8,dpi = 300)



# -------------------------------------------------------
#              Survival Analysis
# -------------------------------------------------------

# 08-Prepare data

meta_clean$time <- meta_clean$time to relapse or last follow-up (months)

matrix_surv <- t(hub_matrix_ready)

all(rownames(matrix_surv) %in% rownames(meta_clean))

surv_obj <- Surv(
  time = meta_clean$time to relapse or last follow-up (months),
  event = meta_clean$relapse
)

# 09-Run cox regression function

cox_results <- run_cox(matrix_surv, meta_clean)

cox_results$significant <- cox_results$pvalue < 0.05

cox_plot <- ggplot(cox_results, aes(x = HR, y = reorder(Gene, HR), xmin = lower_CI, xmax = upper_CI, color = significant)) +
  geom_point(size = 3) +
  geom_errorbarh(height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "red"), guide = "none") +
  theme_bw() +
  xlab("Hazard Ratio (HR)") +
  ylab("Gene") +
  theme(
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major.y = element_blank()
  )
ggsave("04-Cox_HR.png", plot = cox_plot, width = 8, height = 6, dpi = 300)

message("NOTE: Survival analysis workflow is complete up to data preparation. ")
message("Access to the Genomic Data Commons (GDC) portal is required for full datasets, which is not available. Workflow stops here.")

