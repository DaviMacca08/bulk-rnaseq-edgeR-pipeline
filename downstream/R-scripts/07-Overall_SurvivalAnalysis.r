library(tidyverse)
library(ggrepel)
library(survival)
library(survminer)
library(org.Hs.eg.db)
library(edgeR)
library(grid)
library(gridExtra)


# ============================================
#        Survival Analysis Pipeline
# ============================================

# -------------------------------------
# 1. Define Helper Functions
# -------------------------------------

# Function: Generate Kaplan-Meier survival curve for a single gene
km_curve <- function(gene, expr_matrix, time, status) {
  
  # Get expression values for this gene
  expr <- expr_matrix[[gene]]
  
  # Create a data frame with survival info
  df <- data.frame(
    time = time,
    status = status,
    expr = expr
  )
  
  # Find the optimal cutpoint to split patients into high/low expression groups
  cut <- survminer::surv_cutpoint(
    df,
    time = "time",
    event = "status",
    variables = "expr"
  )
  
  # Split patients into high and low expression groups
  df_cut <- survminer::surv_categorize(cut)
  
  # Fit the Kaplan-Meier model
  fit <- survival::survfit(Surv(time, status) ~ expr, data = df_cut)
  
  # Create the survival plot
  plot_km <- survminer::ggsurvplot(
    fit,
    data = df_cut,
    pval = TRUE,
    pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    title = gene,
    xlab = "Time (months)",
    ylab = "Survival probability",
    palette = c("#E64B35", "#4DBBD5"),
    ggtheme = theme_minimal(),
    risk.table.height = 0.25,
    risk.table.fontsize = 3,
    break.time.by = 12
  )
  
  # Get summary data
  km_data <- survminer::surv_summary(fit, data = df_cut)
  
  # Return everything as a list
  return(list(
    fit = fit,
    plot = plot_km,
    data = df_cut
  ))
}

# Function: Run univariate Cox regression on all genes
cox_uni <- function(expr_matrix, time, status) {
  
  results <- data.frame()
  
  # Loop through each gene
  for (gene in colnames(expr_matrix)) {
    
    expr <- expr_matrix[[gene]]
    
    df <- data.frame(
      time = time,
      status = status,
      expr = expr
    )
    
    # Fit Cox model
    fit <- survival::coxph(Surv(time, status) ~ expr, data = df)
    
    # Extract results
    summary_fit <- summary(fit)
    
    results <- dplyr::bind_rows(
      results,
      data.frame(
        gene = gene,
        HR = summary_fit$coefficients[, "exp(coef)"],
        pvalue = summary_fit$coefficients[, "Pr(>|z|)"]
      )
    )
  }
  
  return(results)
}

# Function: Univariate Cox regression for clinical variables
cox_uni_clinical <- function(df, vars) {
  
  cox_list <- list()
  
  for (v in vars) {
    
    # Fit model for each variable
    fit <- coxph(as.formula(paste0("Surv(time, status) ~ ", v)), data = df)
    
    s <- summary(fit)
    
    # Extract hazard ratios and confidence intervals
    temp <- data.frame(
      variable = rownames(s$coefficients),
      HR = s$coefficients[, "exp(coef)"],
      lower = s$conf.int[, "lower .95"],
      upper = s$conf.int[, "upper .95"],
      pvalue = s$coefficients[, "Pr(>|z|)"]
    )
    
    cox_list[[v]] <- temp
  }
  
  res <- dplyr::bind_rows(cox_list)
  res$FDR <- p.adjust(res$pvalue, method = "BH")
  
  return(res)
}

# Function: Multivariate Cox regression with multiple clinical variables 
cox_multi_clinical <- function(df,vars) {

  # Build formula with all variables 
  formula <- as.formula(
    paste0("Surv(time, status) ~ ", paste(vars, collapse = " + "))
  )

  # Fit the model 
  multi_fit <- coxph(formula, data = df)
  
  # Check proportional hazards assumption
  ph_test <- cox.zph(multi_fit)
  ph_plot <- survminer::ggcoxzph(ph_test)
  multi_summary <- summary(multi_fit)
  
  # Extract results
  multi_results <- data.frame(
    variable = rownames(multi_summary$coefficients),
    HR = multi_summary$coefficients[, "exp(coef)"],
    lower = multi_summary$conf.int[, "lower .95"],
    upper = multi_summary$conf.int[, "upper .95"],
    pvalue = multi_summary$coefficients[, "Pr(>|z|)"]
  )
  
  multi_results$FDR <- p.adjust(multi_results$pvalue, method = "BH")
  
  return(list(
    model = multi_fit,
    ph_test = ph_test,
    ph_plot = ph_plot,
    results = multi_results
  ))
}

# Function: Create a combined forest plot with results table
forest_plot <- plot_cox_forest_table <- function(df_forest_uni) {
  
  # Format the results table
  table_df <- df_forest_uni |>
    mutate(
      # Format HR with confidence interval 
      `HR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", HR, lower, upper),
      
      # Format p-values nicely
      FDR = ifelse(FDR < 0.001, "<0.001", sprintf("%.3f", FDR)),
      
      # Add significance stars
      Significance = case_when(
        FDR < 0.001 ~ "***",
        FDR < 0.01 ~ "**",
        FDR < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) |>
    dplyr::select(variable, `HR (95% CI)`, FDR, Significance)
  
  # Create the table grob (graphical object)
  table_grid <- tableGrob(
    table_df,
    rows = NULL,
    theme = ttheme_minimal(
      core = list(
        fg_params = list(fontsize = 14, col = "black"),
        padding = unit(c(5, 5), "mm")
      ),
      colhead = list(
        fg_params = list(fontface = "bold", fontsize = 13),
        padding = unit(c(6, 6), "mm")
      )
    )
  )
  
  # Create the forest plot
  forest_plot <- ggplot(df_forest_uni, aes(x = variable, y = HR)) +
    geom_point(size = 4, shape = 18) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
    geom_hline(yintercept = 1, linetype = "dashed", col = "red") +
    coord_flip() +
    labs(x = "", y = "Hazard Ratio (95% CI)") +
    theme_classic(base_size = 13)
  
  #Combine both into one figure
  grid.arrange(
    forest_plot,
    table_grid,
    ncol = 2,
    widths = c(2, 2)
  )
}


# ============================================
#    2. Load and Prepare the Data
# ============================================

# Set up directory paths

base_dir <- "downstream/Dataset_SurvivalAnalysis"

# Make sure all necessary directories exist
if(!dir.exists(base_dir)) {
  cat("\nCreating directory:", base_dir, "\n")
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
}

# Create directory for plots
dir.create("downstream/plots/TCGA-BRCA", recursive = TRUE, showWarnings = FALSE)

# Load theb expression dataset (note: fixed typo in filename)
brca_dataset <- read.table(
  file.path(base_dir, "TGCA-BRCA_dataset.csv"),
  header = TRUE,
  sep = ",",
  comment.char = "!",
  stringsAsFactors = FALSE
)

# Tidy up the expression matrix
rownames(brca_dataset) <- brca_dataset[,1]
brca_dataset <- brca_dataset[,-1]

# Load che clinical metadata
brca_meta <- readRDS("downstream/Dataset_SurvivalAnalysis/TCGA-BRCA_metadata.rds")

# -------------------------------------
# 2.1 Check the Raw Data
# -------------------------------------

# Look at first 10 samples to check normalization
check_norm <- brca_dataset[, 1:10]

# Save a boxplot to check the data distribution
png("downstream/plots/TCGA-BRCA/01-CheckNormalization.png", 
    width = 1200, height = 1500, res = 150)
boxplot(
  check_norm,
  las = 2, 
  cex.axis = 0.8,
  outline = FALSE,
  ylab = "Counts",
  boxwex = 1,
  col = "lightblue",
  main = "Check TCGA-BRCA Dataset"
)
dev.off()

message("The counts matrix is not normalized - proceeding with normalization.")

# -------------------------------------
# 2.2 Normalize the Data (TMM method)
# -------------------------------------

# Create a DGEList object for edgeR
dgelist <- DGEList(brca_dataset)

# Filter out lowly expressed genes
keep <- filterByExpr(dgelist)
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

# Calculate normalization factors using TMM
normFact <- calcNormFactors(dgelist, method = "TMM")

# Get log2 counts per million
logCPM <- cpm(normFact, log = TRUE)

# Check the normalized data
check_norm2 <- logCPM[, 1:10]

png("downstream/plots/TCGA-BRCA/02-Normalization.png", 
    width = 1200, height = 1500, res = 150)
boxplot(
  check_norm2,
  las = 2, 
  cex.axis = 0.8,
  outline = FALSE,
  ylab = "Counts",
  boxwex = 1,
  col = "lightcoral",
  main = "TCGA-BRCA Dataset after normalization"
)
dev.off()


# ============================================
#    3. Quality Control with PCA
# ============================================

# Run PCA on the normalized data
pca <- prcomp(t(logCPM), scale. = TRUE)

# Calculate variance explained by each PC
pca_var <- pca$sdev^2
pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)
cum_var <- cumsum(pca_var_percent)

pca_scores <- pca$x

# Calculate distance from center to find outliers
dist_center <- sqrt(rowSums(pca_scores[, 1:2]^2))

# Identify outlier samples (more than 3 SD from center)
outliers <- names(dist_center[dist_center > mean(dist_center) + 3 * sd(dist_center)])

# Report any outliers found
cat("Potential outliers (visual QC only):\n", paste(outliers, collapse = ", "), "\n")

# Create a data frame for plotting
pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  sample = colnames(logCPM),
  outlier = colnames(logCPM) %in% outliers
)

# Generate PCA plot with outliers highlighted
pca_plot <- ggplot(pca_df, aes(PC1, PC2, color = outlier)) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(
    data = subset(pca_df, outlier == TRUE),
    aes(label = sample),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "#D55E00"),
    labels = c("FALSE" = "Inliers", "TRUE" = "Outliers")
  ) +
  labs(
    title = "PCA of TCGA-BRCA samples (QC only, no filtering)",
    x = paste0("PC1 (", pca_var_percent[1], "% variance)"),
    y = paste0("PC2 (", pca_var_percent[2], "% variance)"),
    color = "Sample status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(filename = "downstream/plots/TCGA-BRCA/03-PCA_outliers.png", 
       plot = pca_plot, width = 8, height = 6, dpi = 150)


# ============================================
#    4. Clean and Match Metadata
# ============================================

# Select only the clinical features we need
meta <- brca_meta |>
  dplyr::select(
    vital_status,
    days_to_death,
    days_to_last_follow_up,
    bcr_patient_barcode,
    age_at_diagnosis,
    ajcc_pathologic_stage,
    ajcc_pathologic_t,
    gender
    )

# Keep only primary tumor samples (code "01)
keep <- substr(colnames(logCPM), 14, 15) == "01"
logCPM <- logCPM[, keep]

# Extract patient IDs from column names
patients <- substr(colnames(logCPM), 1, 12)
colnames(logCPM) <- patients

# Check for any duplicate patient IDs
sum(duplicated(patients))

# Remove duplicates if any
keep_patients <- !duplicated(patients)
logCPM <- logCPM[, keep_patients]
patients <- patients[keep_patients]

# Standardize patient barcode format in metadata
meta$bcr_patient_barcode <- gsub("-", ".", meta$bcr_patient_barcode)
meta$bcr_patient_barcode <- substr(meta$bcr_patient_barcode, 1, 12)

# Clean up metadata
meta <- meta[!is.na(meta$bcr_patient_barcode), ]
meta <- meta[!duplicated(meta$bcr_patient_barcode), ]

# Find patients that are in both expression and metadata
common <- intersect(patients, meta$bcr_patient_barcode)

# Subset both datasets to only common patients 
logCPM <- logCPM[, patients %in% common]
meta <- meta[meta$bcr_patient_barcode %in% common, ]

# Make sure the order matches
patients <- substr(colnames(logCPM), 1, 12)
meta <- meta[match(patients, meta$bcr_patient_barcode), ]

# Final verification
stopifnot(all(colnames(logCPM) == meta$bcr_patient_barcode))
message("Data matching complete - all patient IDs aligned.")


# ============================================
#    5. Prepare Expression Data for Analysis
# ============================================

# Convert ENSEMBL IDs to gene symbols 
rownames(logCPM) <- gsub("\\..*", "", rownames(logCPM))

symbol_tcga <- mapIds(
  org.Hs.eg.db,
  keys = rownames(logCPM),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

logCPM_df <- as.data.frame(logCPM)

# Add gene symbols 
logCPM_df$symbol <- symbol_tcga

# Remove genes that couldn't be mapped
logCPM_df <- logCPM_df[!is.na(logCPM_df$symbol), ]

# If multiple ENSEMBL IDs map to the same symbol, take the median
logCPM_symbol <- logCPM_df |>
  group_by(symbol) |>
  summarise(
    across(where(is.numeric), median, na.rm = TRUE),
    .groups = "drop"
  )

logCPM_symbol <- as.data.frame(logCPM_symbol)
rownames(logCPM_symbol) <- logCPM_symbol$symbol

# Remove the symbol column (now it's rownames)
if ("symbol" %in% colnames(logCPM_symbol)) {
  logCPM_symbol$symbol <- NULL
  message("Gene symbol column removed - symbols are now rownames.")
}

# -------------------------------------
# 5.1 Load and Filter Hub Genes
# -------------------------------------

# Load the hub gene list form network analysis
hub <- read_csv(file.path(base_dir, "MCC_ValueClean.csv"), col_names = TRUE)

# Sort by importance score (highest first)
hub <- hub |> arrange(desc(Score))

# Extract expression data for hub genes only 
hub_matrix <- logCPM_symbol[rownames(logCPM_symbol)%in% hub$Name, ]
hub_matrix <- as.data.frame(hub_matrix)

# Check if any hub genes are missing from our expression data
missing_genes <- setdiff(hub$Name, rownames(logCPM_symbol))

if (length(missing_genes) > 0) {
  cat("Missing hub genes in expression matrix:\n")
  writeLines(missing_genes)
} else {
  cat("All hub genes are present in the expression matrix.\n")
}

# Transpose so genes are columns and patients are rows
hub_matrix_t <- t(hub_matrix)
hub_matrix_t <- as.data.frame(hub_matrix_t)

# Final check that everything aligns
stopifnot(all(rownames(hub_matrix_t) == meta$bcr_patient_barcode))
message("Hub gene expression matrix ready for survival analysis.")


# ============================================
#    6. Survival Analysis
# ============================================

# -------------------------------------
# 6.1 Define Survival Time and Events
# -------------------------------------

# Calculate overall survival time in months
# For deceased patients: use days_to_death
# For living patients: use days_to_last_follow_up
os_time <- ifelse(
  !is.na(meta$days_to_death),
  meta$days_to_death,
  meta$days_to_last_follow_up
) / 30.44

# Define event status (1 = deceased, 0 = censored)
os_event <- ifelse(meta$vital_status == "Dead", 1, 0)

# Print summary statistics
cat("\nSurvival data summary:\n")
cat("  - Total patients:", length(os_time), "\n")
cat("  - Deceased patients:", sum(os_event == 1), "\n")
cat("  - Censored patients:", sum(os_event == 0), "\n")
cat("  - Missing values:", sum(is.na(os_time)), "\n")

# Create survival object
surv_object <- Surv(time = os_time, event = os_event)

# -------------------------------------
# 6.2 Generate Kaplan-Meier Curves
# -------------------------------------

results_km <- list()
plots_km <- list()
km_data <- list()

cat("\nGenerating Kaplan-Meier curves for hub genes...\n")

for (i in seq_along(colnames(hub_matrix_t))) {
  
  gene <- colnames(hub_matrix_t)[i]
  
  # Progress update
  message(sprintf(" Processing gene %d/%d: %s", i, ncol(hub_matrix_t), gene))
  
  res_km <- km_curve(
    gene = gene,
    expr_matrix = hub_matrix_t,
    time = os_time,
    status = os_event
  )
  
  results_km[[gene]] <- res_km$fit
  plots_km[[gene]] <- res_km$plot
}

cat("Kaplan-Meier curves generated for", length(plots_km), "genes.\n")

# -------------------------------------
# 6.3 Univariate Cox Regression (All Genes)
# ------------------------------------- 

cat("\nRunning univariate Cox regression on all hub genes...\n")

cox_results <- cox_uni(
  expr_matrix = hub_matrix_t,
  time = os_time,
  status = os_event
)

# rank genes by p-value and adjust for multiple testing
cox_results <- cox_results |> 
  arrange(pvalue)

cox_results$FDR <- p.adjust(cox_results$pvalue, method = "BH")

# Identify significant genes
significant_genes <- cox_results |> 
  filter(FDR < 0.05)

if (nrow(significant_genes) > 0 ){
  cat("\nSignificant genes (FDR < 0.05):\n")
  print(significant_genes[, c("gene", "HR", "pvalue", "FDR")])
} else {
  cat("\nNo genes reached significance at FDR < 0.05.\n")
}

# -------------------------------------
# 6.4 Plot KM Curve for Top Gene (PGK1)
# ------------------------------------- 

if("PGK1" %in% names(plots_km)){
  png("downstream/plots/TCGA-BRCA/04-KaplanMeier_PGK1.png",
      width = 1500, height = 1500, res = 150)
  print(plots_km[["PGK1"]])
  dev.off()
} else {
  cat("\nWARNING: PGK1 not found in hub genes. Skipping KM plot.\n")
}


# ============================================
#    7. Clinical Cox Regression Analysis
# ============================================

# Only proceed if PGK1 is in our dataset
if ("PGK1" %in% colnames(hub_matrix_t)) {
  
  cat("\n--- Starting clinical Cox regression analysis ---\n")
  
  # -------------------------------------
  # 7.1 Prepare Clinical Variables
  # -------------------------------------
  
  # Combine gene expression with clinical variables
  df_clinical <- data.frame(
    time = os_time,
    status = os_event,
    PGK1 = hub_matrix_t[["PGK1"]],
    age = meta$age_at_diagnosis,
    stage = meta$ajcc_pathologic_stage,
    stageT = meta$ajcc_pathologic_t
  )
  
  # Simplify AJCC stage into major groups
  df_clinical$stage_clean <- forcats::fct_collapse(
    df_clinical$stage,
    I   = c("Stage I", "Stage IA", "Stage IB"),
    II  = c("Stage II", "Stage IIA", "Stage IIB"),
    III = c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"),
    IV  = "Stage IV"
  ) 
  
  df_clinical$stage_clean <- droplevels(df_clinical$stage_clean)
  
  # Simplify T stage into major groups
  df_clinical <- df_clinical[!is.na(df_clinical$stageT), ]
  df_clinical <- df_clinical[df_clinical$stageT != "TX", ]
  
  df_clinical$stageT_clean <- dplyr::case_when(
    df_clinical$stageT %in% c("T1", "T1a", "T1b", "T1c") ~ "T1",
    df_clinical$stageT %in% c("T2", "T2a", "T2b") ~ "T2",
    df_clinical$stageT == "T3" ~ "T3",
    df_clinical$stageT %in% c("T4", "T4b", "T4d") ~ "T4",
    TRUE ~ NA_character_
  )
  
  # Convert to factors and set reference levels
  df_clinical$stage_clean <- factor(df_clinical$stage_clean)
  df_clinical$stageT_clean <- factor(df_clinical$stageT_clean)
  
  df_clinical$stage_clean <- relevel(df_clinical$stage_clean, ref = "I")
  df_clinical$stageT_clean <- relevel(df_clinical$stageT_clean, ref = "T1")
  
  # Remove patients with missing values
  df_clinical <- df_clinical[complete.cases(
    df_clinical[, c("time","status","PGK1","age","stage_clean","stageT_clean")]
    ), ]
  
  cat("\nPatients remaining after filtering:", nrow(df_clinical), "\n")
  
  # -------------------------------------
  # 7.2 Univariate Clinical Cox Regression
  # -------------------------------------
  
  vars <- c("PGK1", "age", "stage_clean", "stageT_clean")
  
  uni_cox_res_clinical <- cox_uni_clinical(df_clinical, vars)
  
  # Clean up variable names for the plot
  uni_cox_res_clinical <- uni_cox_res_clinical |> 
    mutate(
      variable = dplyr::recode(variable,
                               "PGK1" = "PGK1",
                               "age" = "Age",
                               "stage_cleanII" = "Stage II",
                               "stage_cleanIII" = "Stage III",
                               "stage_cleanIV" = "Stage IV",
                               "stageT_cleanT2" = "T2 Stage",
                               "stageT_cleanT3" = "T3 Stage",
                               "stageT_cleanT4" = "T4 Stage"
      )
    )
  
  uni_cox_res_clinical$variable <- factor(
    uni_cox_res_clinical$variable,
    levels = c("PGK1", "Age",
               "Stage II", "Stage III", "Stage IV",
               "T2 Stage", "T3 Stage", "T4 Stage")
  )
  
  # Prepare data for forest plot (including all stages)
  df_forest_uni_clinical <- uni_cox_res_clinical |> 
    dplyr::filter(variable %in% c(
      "PGK1", "Age", 
      "Stage II", "Stage III", "Stage IV",
      "T2 Stage", "T3 Stage", "T4 Stage"
      ))
  
  df_forest_uni_clinical$variable <- forcats::fct_rev(df_forest_uni_clinical$variable)
  
  # Generate univariate forest plot
  png("downstream/plots/TCGA-BRCA/05-ForestPlot_Univariate.png", 
      width = 2000, height = 1500, res = 150)
  forestPlot_uni_clinical <- forest_plot(df_forest_uni_clinical)
  dev.off()
  
  message("Univariate forest plot saved.")
  
  # -------------------------------------
  # 7.3 Multivariate Clinical Cox Regression
  # -------------------------------------
  
  cox_multi <- cox_multi_clinical(df_clinical, vars)
  
  # Display model results 
  cat("\nMultivariate Cox regression results:\n")
  print(cox_multi$results)
  
  # Save proportional hazards diagnostic plot
  png("downstream/plots/TCGA-BRCA/06-phPlot.png",
      width = 1500, height = 1500, res = 150)
  cox_multi$ph_plot
  dev.off()
  
  # Clean up names for forest plot
  multi_cox_res_clinical <- cox_multi$results |>
    mutate(
      variable = dplyr::recode(variable,
        "PGK1" = "PGK1",
        "age" = "Age",
        "stage_cleanII" = "Stage II",
        "stage_cleanIII" = "Stage III",
        "stage_cleanIV" = "Stage IV",
        "stageT_cleanT2" = "T2 Stage",
        "stageT_cleanT3" = "T3 Stage",
        "stageT_cleanT4" = "T4 Stage"
      )
    )
  
  df_forest_multi_clinical <- multi_cox_res_clinical |>
    filter(variable %in% c(
      "PGK1", "Age",
      "Stage II", "Stage III", "Stage IV",
      "T2 Stage", "T3 Stage", "T4 Stage"
    ))
 
  df_forest_multi_clinical$variable <- forcats::fct_rev(df_forest_multi_clinical$variable)
  
  # Generate multivariate forest plot
  png("downstream/plots/TCGA-BRCA/07-ForestPlot_Multivariate.png", 
      width = 2000, height = 1500, res = 150)
  forestPlot_multi_clinical <- forest_plot(df_forest_multi_clinical)
  dev.off()
  
  message("Multivariate forest plot saved.")
} else {
  cat("\nWARNING: PGK1 not found in hub genes. Skipping clinical Cox regression.\n")
}


# ============================================
#    8. Save Session Information
# ============================================

sink("TCGA-BRCA_SurvivalAnalysis_info.txt")
cat("Survival Analysis on TCGA-BRCA dataset - Session Info\n")
cat("=====================================================\n\n")
cat("Analysis completed:", as.character(Sys.time()), "\n\n")
sessionInfo()
sink()








