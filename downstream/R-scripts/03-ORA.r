library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(patchwork)



# ------------------------------------------------
#        Over-representation Analysis (ORA) 
# ------------------------------------------------

# Prepare functions
# 01-Function to add Entrez IDs with error checking
add_entrez <- function(df, id_col = "GeneId") {
  if(is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }
  
  df$EntrezID <- mapIds(org.Hs.eg.db,
                        keys = df[[id_col]],
                        keytype = "ENSEMBL",
                        column = "ENTREZID",
                        multiVals = "first")
  return(df)
}

# 02-DotPlot function with FoldEnrichment calculation
plot_enrichment <- function(df, title, color) {
  if(is.null(df) || nrow(df) == 0) {
    return(NULL)
  }
  
  # -Calculate FoldEnrichment if not present
  if(!"FoldEnrichment" %in% colnames(df)) {
    
    message("Calculate Fold Enrichment")
    gene_ratio_num <- as.numeric(sub("/.*", "", df$GeneRatio))
    gene_ratio_den <- as.numeric(sub(".*/", "", df$GeneRatio))
      
    bg_ratio_num <- as.numeric(sub("/.*", "", df$BgRatio))
    bg_ratio_den <- as.numeric(sub(".*/", "", df$BgRatio))
      
    df$FoldEnrichment <- (gene_ratio_num / gene_ratio_den) / (bg_ratio_num / bg_ratio_den)
    }
  
  df$Set <- factor(df$Set, levels = c("Module", "DEG", "Intersection"))
  
  if(!"ONTOLOGY" %in% colnames(df)){
    Top<- df |>
      group_by(Set) |>
      arrange(p.adjust) |>
      slice_head(n = 20) |>
      ungroup()
  } else{
    Top <- df |> 
      group_by(Set, ONTOLOGY) |> 
      arrange(p.adjust) |> 
      slice_head(n = 10) |> 
      ungroup()
  }
  
  ggplot(Top, aes(x = Set, y = Description)) +
    geom_point(aes(size = FoldEnrichment, color = p.adjust)) +
    scale_color_gradient(low = color, high = "gray60",
                        labels = scales::scientific) +
    scale_size_continuous(range = c(2, 6), breaks = pretty(Top$FoldEnrichment, n = 3)) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5)) +
    labs(title = title,
         x = "",
         y = "",
         size = "Fold Enrichment",
         color = "FDR")
}

# 03-Function to run enrichGO with error handling
run_enrichGO <- function(gene_list, background_list, set_name, min_genes = 10) {
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  result <- tryCatch({
    ego <- enrichGO(gene = gene_list,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   universe = background_list,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   minGSSize = 10,
                   maxGSSize = 500)
    
    if (is.null(ego) || nrow(ego) == 0) {
      message(set_name, ": no significant GO terms found")
      return(NULL)
    }
    
    df <- as.data.frame(ego)
    df$Set <- set_name
    return(df)
  }, error = function(e) {
    message(set_name, ": enrichment failed - ", e$message)
    return(NULL)
  })
  
  return(result)
}

# 04-Function to run enrichKEGG with error handling
run_enrichKEGG <- function(gene_list, background_list, set_name, min_genes = 10) {
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for KEGG enrichment")
    return(NULL)
  }
  
  result <- tryCatch({
    kk <- enrichKEGG(gene = gene_list,
                    organism = "hsa",
                    keyType = "ncbi-geneid",
                    universe = background_list,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05,
                    minGSSize = 10,
                    maxGSSize = 500)
    
    if (is.null(kk) || nrow(kk) == 0) {
      message(set_name, ": no significant KEGG pathways found")
      return(NULL)
    }
    
    df <- as.data.frame(kk)
    df$Set <- set_name
    return(df)
  }, error = function(e) {
    message(set_name, ": KEGG enrichment failed - ", e$message)
    return(NULL)
  })
  
  return(result)
}

# 05-Function to run Reactome enrichment with error handling
run_enrichReactome <- function(gene_list, background_list, set_name, min_genes = 10) {
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for Reactome enrichment")
    return(NULL)
  }
  
  result <- tryCatch({
    rp <- enrichPathway(gene = gene_list,
                        organism = "human",
                        universe = background_list,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        minGSSize = 10,
                        maxGSSize = 500,
                        readable = TRUE)
    
    if (is.null(rp) || nrow(rp) == 0) {
      message(set_name, ": no significant Reactome pathways found")
      return(NULL)
    }
    
    df <- as.data.frame(rp)
    df$Set <- set_name
    return(df)
  }, error = function(e) {
    message(set_name, ": Reactome enrichment failed - ", e$message)
    return(NULL)
  })
  
  return(result)
}





# ------------------------------------------------
#                Start Analysis 
# ------------------------------------------------

# 01-Define path

base_dir <- "downstream/objects"
deg_file <- file.path(base_dir, "ExactTestDEG.rds")
genes_ME_file <- file.path(base_dir, "genes_ME_brown.rds")
background_file <- file.path(base_dir, "Background_genes.rds")

# 02-Load data with error checking

load_data <- function(file, name) {
  if (!file.exists(file)) {
    stop(name, " not found at: ", file, "\n",
         "Please run previous scripts first.")
  }
  readRDS(file)
}

cat("Loading data...\n")
deg <- load_data(deg_file, "ExactTestDEG")
genes_ME <- load_data(genes_ME_file, "genes_ME_brown")
background <- load_data(background_file, "Background_genes")

# 03-Add Entrez IDs

background <- add_entrez(background)
genes_ME <- add_entrez(genes_ME)

# 04-Remove NA and get unique

list_deg <- unique(deg$entrez[!is.na(deg$entrez)])
list_genes_ME <- unique(genes_ME$EntrezID[!is.na(genes_ME$EntrezID)])
list_background <- unique(background$EntrezID[!is.na(background$EntrezID)])

# 05-Common genes

commonGenes <- data.frame(
  GeneId = intersect(genes_ME$GeneId, deg$Geneid)
)
commonGenes <- add_entrez(commonGenes)
list_common <- unique(commonGenes$EntrezID[!is.na(commonGenes$EntrezID)])

cat("\n--- Dataset sizes ---\n",
    "DEGs:", length(list_deg), "\n",
    "Module genes:", length(list_genes_ME), "\n",
    "Intersection:", length(list_common), "\n",
    "Background:", length(list_background), "\n\n")



# ------------------------------------------------
#               GO Terms - enrichGO() 
# ------------------------------------------------

cat("Running GO enrichment analysis...\n")

# 06-Run GO enrichment

GO_deg <- run_enrichGO(list_deg, list_background, "DEG")
GO_ME <- run_enrichGO(list_genes_ME, list_background, "Module")
GO_common <- run_enrichGO(list_common, list_background, "Intersection")

# Combine results 

GO_results <- list(GO_deg, GO_ME, GO_common)
GO_results <- GO_results[!sapply(GO_results, is.null)]

if (length(GO_results) > 0) {
  df_GO <- bind_rows(GO_results)
  
  # -Create plots
  p_bp <- plot_enrichment(filter(df_GO, ONTOLOGY == "BP"), 
                          "GO Biological Process", "firebrick")
  p_mf <- plot_enrichment(filter(df_GO, ONTOLOGY == "MF"), 
                          "GO Molecular Function", "darkblue")
  p_cc <- plot_enrichment(filter(df_GO, ONTOLOGY == "CC"), 
                          "GO Cellular Component", "darkgreen")
  
  # -Arrange and save plots 
  plots <- list(p_bp, p_mf, p_cc)
  plots <- plots[!sapply(plots, is.null)]
  
  if (length(plots) > 0) {
    # -Actually display the plots
    combined_plot <- wrap_plots(plots, ncol = 1) + 
      plot_annotation(title = "GO Enrichment Analysis",
                     theme = theme(plot.title = element_text(hjust = 0.5)))
    
    print(combined_plot)
    ggsave("01-GO_enrichment.png", combined_plot, width = 12, height = 18, dpi = 300)
    cat("GO enrichment plot saved to GO_enrichment.png\n")
  }
} else {
  cat("No significant GO terms found for any gene set\n")
}



# ------------------------------------------------
#               KEGG Pathway Analysis
# ------------------------------------------------

cat("\nRunning KEGG pathway analysis...\n")

# 07-Run KEGG

kegg_deg <- run_enrichKEGG(list_deg, list_background, "DEG")
kegg_ME <- run_enrichKEGG(list_genes_ME, list_background, "Module")
kegg_common <- run_enrichKEGG(list_common, list_background, "Intersection")

# Combine and plot KEGG

kegg_results <- list(kegg_deg, kegg_ME, kegg_common)
kegg_results <- kegg_results[!sapply(kegg_results, is.null)]

if (length(kegg_results) > 0) {
  df_kegg <- bind_rows(kegg_results)
  
  p_kegg <- plot_enrichment(df_kegg, "KEGG Pathway Enrichment", "firebrick")
  if (!is.null(p_kegg)) {
    print(p_kegg)
    ggsave("02-KEGG_enrichment.png", p_kegg, width = 10, height = 8, dpi = 300)
    cat("KEGG enrichment plot saved to KEGG_enrichment.png\n")
  }
} else {
  cat("No significant KEGG pathways found for any gene set\n")
}



# ------------------------------------------------
#               Reactome Analysis
# ------------------------------------------------

cat("\nRunning Reactome pathway analysis...\n")

# 08-Run Reactome

reactome_deg <- run_enrichReactome(list_deg, list_background, "DEG")
reactome_ME <- run_enrichReactome(list_genes_ME, list_background, "Module")
reactome_common <- run_enrichReactome(list_common, list_background, "Intersection")

# Combine and plot Reactome

reactome_results <- list(reactome_deg, reactome_ME, reactome_common)
reactome_results <- reactome_results[!sapply(reactome_results, is.null)]

if (length(reactome_results) > 0) {
  df_reactome <- bind_rows(reactome_results)
  
  p_reactome <- plot_enrichment(df_reactome, "Reactome Pathway Enrichment", "steelblue")
  if (!is.null(p_reactome)) {
    print(p_reactome)
    ggsave("03-Reactome_enrichment.png", p_reactome, width = 14, height = 8, dpi = 300)
    cat("Reactome enrichment plot saved to Reactome_enrichment.png\n")
  }
} else {
  cat("No significant Reactome pathways found for any gene set\n")
}

# 09-Save session info

sink("ORA_info.txt")
cat("ORA Analysis Session Info\n")
cat("===============================\n\n")
sessionInfo()
sink()

cat("\nORA analysis completed successfully!\n")
cat("Results saved to current working directory\n")
