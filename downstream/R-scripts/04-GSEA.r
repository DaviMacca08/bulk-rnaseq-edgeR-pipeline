library(tidyverse)
library(BiocParallel)
library(patchwork)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(pathview)
library(msigdbr)
library(fgsea)

register(SerialParam())
set.seed(123)



# ------------------------------------------------
#      Gene-Set enrichment analysis (GSEA) 
# ------------------------------------------------

# Prepare Functions 
## 01-Load data with error checking
base_dir <- "downstream/objects"
resultsExactTest <- file.path(base_dir, "ResFromExactTest.rds")
degExactTest <- file.path(base_dir, "ExactTestDEG.rds")

load_data <- function(file, name) {
  if (!file.exists(file)) {
    stop(name, " not found at: ", file, "\n",
         "Please run previous scripts first.")
  }
  readRDS(file)
}

## 02-gseGO() and plot for GO terms functions with error handling
run_gseGO_BP <- function(gene_list, set_name, min_genes = 10){
  
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  results <- tryCatch({
    BP <- gseGO(
      geneList = gene_list,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      verbose = TRUE)
    
    df <- as.data.frame(BP)
    
  if(is.null(df) || nrow(df) == 0){
    message(set_name, ": no significant GO-BP terms found")
    return(NULL)
  }
  
 
  df$set <- set_name
  df$Ontology <- "BP"
  df <- df |> 
    extract(
      leading_edge,
      into = c("tags", "list", "signal"),
      regex = "tags=(\\d+)%.*, list=(\\d+)%.*, signal=(\\d+)%"
    ) |> 
    mutate(
      tags = as.numeric(tags),
      list = as.numeric(list),
      signal = as.numeric(signal)
    )
  
  return(df)
  }, error = function(e) {
    message(set_name, ": enrichment failed - ", e$message)
    return(NULL)
})
}
run_gseGO_MF <- function(gene_list, set_name, min_genes = 10){
  
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  results <- tryCatch({
    BP <- gseGO(
      geneList = gene_list,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = "MF",
      pAdjustMethod = "BH",
      verbose = TRUE)
    
    df <- as.data.frame(BP)
    
    if(is.null(df) || nrow(df) == 0){
      message(set_name, ": no significant GO-MF terms found")
      return(NULL)
    }
    
    
    df$set <- set_name
    df$Ontology <- "MF"
    df <- df |> 
      extract(
        leading_edge,
        into = c("tags", "list", "signal"),
        regex = "tags=(\\d+)%.*, list=(\\d+)%.*, signal=(\\d+)%"
      ) |> 
      mutate(
        tags = as.numeric(tags),
        list = as.numeric(list),
        signal = as.numeric(signal)
      )
    
    return(df)
  }, error = function(e) {
    message(set_name, ": enrichment failed - ", e$message)
    return(NULL)
  })
}
run_gseGO_CC <- function(gene_list, set_name, min_genes = 10){
  
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  results <- tryCatch({
    BP <- gseGO(
      geneList = gene_list,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = "CC",
      pAdjustMethod = "BH",
      verbose = TRUE)
    
    df <- as.data.frame(BP)
    
    if(is.null(df) || nrow(df) == 0){
      message(set_name, ": no significant GO-CC terms found")
      return(NULL)
    }
    
    
    df$set <- set_name
    df$Ontology <- "CC"
    df <- df |> 
      extract(
        leading_edge,
        into = c("tags", "list", "signal"),
        regex = "tags=(\\d+)%.*, list=(\\d+)%.*, signal=(\\d+)%"
      ) |> 
      mutate(
        tags = as.numeric(tags),
        list = as.numeric(list),
        signal = as.numeric(signal)
      )
    
    return(df)
  }, error = function(e) {
    message(set_name, ": enrichment failed - ", e$message)
    return(NULL)
  })
}

run_GO_Up <- function(df, ontology, direction = "Up"){
  
  up <- df |> 
    filter(NES > 0 & Ontology == ontology & p.adjust < 0.05 & tags > 30 & signal > 20) |> 
    dplyr::arrange(p.adjust) |> 
    head(n = 12)
  
  if(nrow(up) == 0){
    message("No ",ontology, " ", direction, " terms passed the selection threshold")
    return(NULL)
  }
  
  return(up)
}
run_GO_Down <- function(df, ontology, direction = "Down"){
  
  up <- df |> 
    filter(NES < 0 & Ontology == ontology & p.adjust < 0.05 & tags > 30 & signal > 20) |> 
    dplyr::arrange(p.adjust) |> 
    head(n = 12)
  
  if(nrow(up) == 0){
    message("No ",ontology, " ", direction, " terms passed the selection threshold")
    return(NULL)
  }
  
  return(up)
}
plotGO <- function(up, down, ontology, color){
  
  df <- bind_rows(up, down)
  
  if (is.null(df) || nrow(df) == 0) {
    stop("No data available for plotting")
  }
  
  df <- df |> 
    mutate(
      GeneHit = tags,
      FDR = -log10(p.adjust),
      Pathway = reorder(Description, NES)
    )
  
  plot <- ggplot(df, aes(x = NES, y = Pathway)) +
    geom_segment(aes(x = 0, xend = NES, y = Pathway, yend = Pathway), color = "grey70") + 
    geom_point(aes(size = GeneHit, color = FDR)) + 
    scale_color_gradient(low = "gray60", high = color, name = "-log10(FDR)") + 
    scale_size_continuous(name = "Gene hit") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank()
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "",
      title = paste0("GSEA Plot - ", ontology)
    )
  
  return(plot)
}

## 04-gseKEGG()/gsePathway()/ and plot function with error handling
run_KEGG <- function(gene_list, set_name, min_genes = 10){
  
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  results <- tryCatch({
    kegg <- gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      keyType = "ncbi-geneid",
      nPerm = 10000,
      pAdjustMethod = "BH",
      verbose = TRUE)
    
    df <- as.data.frame(kegg)
    
    if(is.null(df) || nrow(df) == 0){
      message(set_name, ": no significant KEGG pathway found")
      return(NULL)
    }
    
    df$Set <- set_name
    
    df <- df |> 
      extract(
        leading_edge,
        into = c("tags", "list", "signal"),
        regex = "tags=(\\d+)%.*, list=(\\d+)%.*, signal=(\\d+)%"
      ) |> 
      mutate(
        tags = as.numeric(tags),
        list = as.numeric(list),
        signal = as.numeric(signal))
    
    return(list(df = df, kegg = kegg))}, error = function(e) {
      message(set_name, ": enrichment failed - ", e$message)
      return(NULL)
      })
}
run_reactome <- function(gene_list, set_name, min_genes = 10){
  
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  results <- tryCatch({
    reactome <- gsePathway(
      geneList = gene_list,       
      organism = "human",
      pAdjustMethod = "BH",
      nPermSimple = 10000,
      minGSSize = 10,               
      maxGSSize = 500,              
      verbose = TRUE
    )
    
    df <- as.data.frame(reactome)
    
    if(is.null(df) || nrow(df) == 0){
      message(set_name, ": no significant KEGG pathway found")
      return(NULL)
    }
    
    df$Set <- set_name
    
    df <- df |> 
      extract(
        leading_edge,
        into = c("tags", "list", "signal"),
        regex = "tags=(\\d+)%.*, list=(\\d+)%.*, signal=(\\d+)%"
      ) |> 
      mutate(
        tags = as.numeric(tags),
        list = as.numeric(list),
        signal = as.numeric(signal)
      )
    
    return(df)} , error = function(e) {
      message(set_name, ": enrichment failed - ", e$message)
      return(NULL)
  } )
}
run_hallmarks <- function(gene_list, set_name, PathList, min_genes = 10){
  
  # -Check minimum gene set size
  if (length(gene_list) < min_genes) {
    message(set_name, ": too few genes (", length(gene_list), ") for enrichment (min = ", min_genes, ")")
    return(NULL)
  }
  
  results <- tryCatch({
    hallmarks <- fgsea(
      pathways = PathList,
      stats = gene_list,
      nperm = 10000
    )
    
    df <- as.data.frame(hallmarks)
    
    if(is.null(df) || nrow(df) == 0){
      message(set_name, ": no significant HallMarks pathway found")
      return(NULL)
    }
    
    df$Set <- set_name
    
    N <- length(gene_list)
    
    df <- df |> 
      mutate(
        leadingEdge_size = purrr::map_int(leadingEdge, length),
        tags   = 100 * leadingEdge_size / size,
        list   = 100 * leadingEdge_size / N,
        signal = tags * (1 - list/100) * (N/(N - size)),
        p.adjust = padj,
        Description = pathway
      )
      
    return(df)}, error = function(e) {
      message(set_name, ": enrichment failed - ", e$message)
      return(NULL)
  })
}

run_Pathway_Up <- function(df, db, direction = "Up"){
  
  up <- df |> 
    filter(NES > 0 & p.adjust < 0.05 & tags >= 30 & signal > 20) |> 
    dplyr::arrange(p.adjust) |> 
    head(n = 12)
  
  if(nrow(up) == 0){
    message("No ",db, " ", direction, " Pathway passed the selection threshold")
    return(NULL)
  }
  
  return(up)
}
run_Pathway_Down <- function(df, db, direction = "Down"){
  
  down <- df |> 
    filter(NES < 0 & p.adjust < 0.05 & tags > 30 & signal > 20) |> 
    dplyr::arrange(p.adjust) |> 
    head(n = 12)
  
  if(nrow(down) == 0){
    message("No ", db, " ", direction, " Pathway passed the selection threshold")
    return(NULL)
  }
  
  return(down)
}

plotDb <- function(up, down,db, color){
  
  df <- bind_rows(up, down)
  
  if(is.null(df) || nrow(df) == 0){
    stop("No data available for plotting")
  }
  df <- df |> 
    mutate(
      GeneHit = tags,
      FDR = -log10(p.adjust),
      Pathway = reorder(Description, NES)
    )
  
  plot <- ggplot(df, aes(x = NES, y = Pathway)) +
    geom_segment(aes(x = 0, xend = NES, y = Pathway, yend = Pathway), color = "grey70") + 
    geom_point(aes(size = GeneHit, color = FDR)) + 
    scale_color_gradient(low = "gray60", high = color, name = "-log10(FDR)")+ 
    scale_size_continuous(name = "Gene hit") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_blank()
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "",
      title = paste0(db, " Pathway GSEA")
    )
  
  return(plot)
}

# 05-Run enrichment score with error handling

run_enrich_score <- function(ID, title, db){
  
  if(!ID %in% names(db@geneSets)) {
    stop("geneSetID not found in the enrichment result")
  }
  
  path <- gseaplot2(db, 
                    geneSetID = ID,
                    color = "blue",
                    title = title,
                    subplots = 1:3
  )
  
  return(path)
}

# 06-Expression of genes involved in the pathways with error handling

ExpressionToPathway <- function(ID, title, db, genes, color){
  
  if(!ID %in% names(db@geneSets)) {
    stop("geneSetID not found in the enrichment result")
  }
  
  genes_path <- db@result$core_enrichment[db@result$Description == title]
  genes_path <- strsplit(genes_path, "/")[[1]]
  
  intersect <- intersect(genes_path, genes$entrez)
  
  genes_to_plot <- genes |> 
    filter(entrez %in% intersect)
  
  if(nrow(genes_to_plot) == 0){
    paste("No intersected genes identified.\n")
    return(NULL)
  }
  
  if(title == "Oxidative phosphorylation"){
    genes_to_plot <- genes_to_plot |> 
      filter(logFC < 0 ) |> 
      arrange(logFC) |> 
      head(n = 30)
  }
  
  #plot
  
  plot <- ggplot(genes_to_plot, aes(y = reorder(symbol, -logFC), x = logFC, fill = -log10(FDR))) +
    geom_col() +
    scale_fill_gradient(low = "gray60", high = color, name = "-log10(FDR)") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() + 
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),          
      axis.text.x = element_text(angle = 45, hjust = 1),  
      panel.grid.major.y = element_blank()
    ) +
    labs(
      x = "log2 Fold Change",
      y = "Gene",
      title = paste0("Expression of genes involved in ", title))
  
  return(plot)
  
}





# ------------------------------------------------
#                Start Analysis 
# ------------------------------------------------

# 01-Load Data 

cat("Loading data...\n")
res <- load_data(resultsExactTest, "ResFromExactTest.rds")
deg <- load_data(degExactTest, "ExactTestDEG.rds")

rank <- data.frame(
  score = res$logFC,
  Entrez = res$entrez)

rank_df <- rank |> 
  filter(!is.na(Entrez))

rank_df <- rank_df |> 
  group_by(Entrez) |> 
  slice_max(abs(score), n = 1, with_ties = FALSE) |> 
  ungroup() |> 
  arrange(desc(score))

list <- rank_df$score
names(list) <- rank_df$Entrez
list <- list[is.finite(list)]
summary(list)

# ------------------------------------------------
#             Gene Ontology database 
# ------------------------------------------------

# 02-Run gseGO() for all category terms

go_bp <- run_gseGO_BP(list, "All Genes")
go_mf <- run_gseGO_MF(list, "All Genes")
go_cc <- run_gseGO_CC(list, "All Genes")

go_df <- bind_rows(go_bp, go_mf, go_cc)

TopGO_BP_Up <- run_GO_Up(go_df, "BP")
TopGO_BP_Down <- run_GO_Down(go_df, "BP")
TopGO_CC_Up <- run_GO_Up(go_df, "CC")
TopGO_CC_Down <- run_GO_Down(go_df, "CC")
TopGO_MF_Up <- run_GO_Up(go_df, "MF")
TopGO_MF_Down <- run_GO_Down(go_df, "MF")

# Bubble plot - Lollipop style

BioProc_Plot <- plotGO(TopGO_BP_Up, TopGO_BP_Down, "Biological Process", "firebrick")
MolFunc_Plot <- plotGO(TopGO_MF_Up, TopGO_MF_Down, "Molecular Functions", "darkblue")
CellComp_Plot <- plotGO(TopGO_CC_Up, TopGO_CC_Down, "Cellular Component", "darkgreen")

plots_GO <- list(BioProc_Plot, MolFunc_Plot, CellComp_Plot)
plots_GO <- plots_GO[!sapply(plots_GO, is.null)]

if (length(plots_GO) > 0) {
  # -Actually display the plots
  combined_plot <- wrap_plots(plots_GO, ncol = 1) + 
    plot_annotation(title = "GO GSEA Analysis",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  print(combined_plot)
  ggsave("01-GO_GSEA.png", combined_plot, width = 12, height = 18, dpi = 300)
  cat("GO enrichment plots saved to GO_enrichment.png\n")
} else {
  cat("No significant GO terms found for any gene set\n")
}



# ------------------------------------------------
#         KEGG/Reactome/HallMarks database 
# ------------------------------------------------

# 03-KEGG
list_objects <- run_KEGG(list, "All Genes")

kegg <- list_objects$df
kegg_obj <- list_objects$kegg

kegg_up <- run_Pathway_Up(kegg, "KEGG")
kegg_down <- run_Pathway_Down(kegg, "KEGG")

# Bubble plot - Lollipop style

Kegg_Plot <- plotDb(kegg_up, kegg_down, "KEGG", "firebrick")

# 04-Reactome

reactome <- run_reactome(list, "All Genes")

react_up <- run_Pathway_Up(reactome, "Reactome")
react_down <- run_Pathway_Down(reactome, "Reactome")

# Bubble plot - Lollipop style

React_Plot <- plotDb(react_up, react_down, "Reactome", "steelblue")

# 05-HallMarks

m <- msigdbr_collections(db_species = "HS")
print(m)


db_hallmarks <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

hallmark_list <- split(
  x = db_hallmarks$ncbi_gene, 
  f = db_hallmarks$gs_name
)

hallmarks <- run_hallmarks(list, "All Genes", hallmark_list)

hm_up <- run_Pathway_Up(hallmarks, "HallMarks")
hm_down <- run_Pathway_Down(hallmarks, "HallMarks")

# Bubble plot - Lollipop style

HallMarks_plot <- plotDb(hm_up, hm_down, "HallMarks", "darkviolet")

# 06-Combine plots from different database

plots_db <- list(Kegg_Plot, React_Plot, HallMarks_plot)

if (length(plots_db) > 0) {
  # -Actually display the plots
  combined_plot2 <- wrap_plots(plots_db, ncol = 1, nrow = 3) + 
    plot_annotation(title = "Pathway from different Databases",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  print(combined_plot2)
  ggsave("02-Kegg-Reactome-HallMarks.png", combined_plot2, width = 12, height = 18, dpi = 300)
  cat("Pathway plots saved to Kegg-Reactome-HallMarks.png\n")
}  else cat("No significant GO terms found for any gene set\n")



# ------------------------------------------------------------------------
#       HIF-1, Glycolysis/Gluconeogenesis, Oxidative phosphorylation
# ------------------------------------------------------------------------

# 07-Calculate enrichment score

ID_hif1 <- kegg |> 
  filter(Description == "HIF-1 signaling pathway") |> 
  pull(ID)

ID_gly_glu <- kegg |> 
  filter(Description == "Glycolysis / Gluconeogenesis") |> 
  pull(ID)
  
ID_ox <- kegg |> 
  filter(Description == "Oxidative phosphorylation") |> 
  pull(ID)

HIF1 <- run_enrich_score(ID_hif1, "HIF-1 signaling pathway", kegg_obj)
Gly_Glu <- run_enrich_score(ID_gly_glu, "Glycolysis / Gluconeogenesis", kegg_obj)
OxPho <- run_enrich_score(ID_ox, "Oxidative phosphorylation", kegg_obj)

plots_EnrichScore <- list(HIF1, Gly_Glu, OxPho)
plot_names <- c("HIF1", "Gly_Glu", "OxPho")
a <- 3
for (i in seq_along(plots_EnrichScore)) {
  png(filename = paste0(sprintf("%02d", a), "-EnrichmentScore_", plot_names[i], ".png"),
      width = 2400, height = 1800, res = 150)
  grid::grid.newpage()
  grid::grid.draw(plots_EnrichScore[[i]])  # disegna il grob intero
  dev.off()
  a <- a + 1
}

# 08-Expression of genes involved in the three selected pathways

Exp_HIF1 <- ExpressionToPathway(ID_hif1, "HIF-1 signaling pathway", kegg_obj, deg, "darkred")
Exp_Gly_Glu <- ExpressionToPathway(ID_gly_glu, "Glycolysis / Gluconeogenesis", kegg_obj, deg, "darkblue")
Exp_OxPho <- ExpressionToPathway(ID_ox, "Oxidative phosphorylation", kegg_obj, res, "darkgreen")

GeneInPath <- list(Exp_HIF1, Exp_Gly_Glu, Exp_OxPho)

if (length(GeneInPath) > 0) {
  # -Actually display the plots
  combined_plot3 <- wrap_plots(GeneInPath, ncol = 1, nrow = 3) + 
    plot_annotation(title = "Genes involved in the selected pathways",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
  
  print(combined_plot3)
  ggsave("06-GeneInvolved.png", combined_plot3, width = 12, height = 18, dpi = 300)
  cat("Gene Involved plots saved to GeneInvolved.png.png\n")
} else cat("No significant genes found for any pathway set\n")

# 09-UpSet Plot of pathways

kegg_subset <- kegg_obj
kegg_subset@result <- kegg_obj[kegg_obj@result$Description %in% c(
  "HIF-1 signaling pathway",
  "Glycolysis / Gluconeogenesis",
  "Oxidative phosphorylation"
), ]

png("07-UpSet_Plot.png",  width = 1000, height = 900, res = 150)
par(mfrow = c(1,1))
upsetplot(kegg_subset)
dev.off()

# 10-PathView 

list_pathview <- res$logFC
names(list_pathview) <- res$entrez

pathview(
  gene.data = list_pathview,
  pathway.id = "hsa04066",
  species = "hsa",
  res = 600,
  cex = 0.4,
  low  = list(gene = "blue"),
  mid  = list(gene = "white"),
  high = list(gene = "red"),
  new.signature = FALSE,
  out.suffix = "08-HIF-1 signaling pathway"
)

pathview(
  gene.data = list_pathview,
  pathway.id = "hsa00010",
  species = "hsa",
  res = 600,
  cex = 0.4,
  low  = list(gene = "blue"),
  mid  = list(gene = "white"),
  high = list(gene = "red"),
  new.signature = FALSE,
  out.suffix = "09-Glycolysis - Gluconeogenesis"
)

pathview(
  gene.data = list_pathview,
  pathway.id = "hsa00190",
  species = "hsa",
  res = 600,
  cex = 0.4,
  low  = list(gene = "blue"),
  mid  = list(gene = "white"),
  high = list(gene = "red"),
  new.signature = FALSE,
  out.suffix = "10-Oxidative phosphorylation"
)

# 11-Save session info

sink("GSEA_info.txt")
cat("GSEA Analysis Session Info\n")
cat("===============================\n\n")
sessionInfo()
sink()

cat("\nGSEA analysis completed successfully!\n")
cat("Results saved to current working directory\n")

