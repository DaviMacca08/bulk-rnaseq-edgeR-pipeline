library(org.Hs.eg.db)
library(AnnotationDbi)



# -------------------------------------------------------
#            Prepare Input for STRING database
# -------------------------------------------------------

# Function for load data

load_data <- function(file, name) {
  if (!file.exists(file)) {
    stop(name, " not found at: ", file, "\n",
         "Please run previous scripts first.")
  }
  readRDS(file)
}
# 01-Load common genes from DEG and Cluster

base_dir <- "downstream/objects"

DEG_Cluster_common<- file.path(base_dir, "CommonGenes.rds")

genes <- load_data(DEG_Cluster_common, "CommonGenes.rds")

list_SymbolID <- data.frame(
  GeneId = genes
)

list_SymbolID$Symbol <- mapIds(
    org.Hs.eg.db,
    keys = list_SymbolID$GeneId,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
)

list_SymbolID <- list_SymbolID[!is.na(list_SymbolID$Symbol), ]
list_SymbolID <- list_SymbolID[!duplicated(list_SymbolID$Symbol), ]

writeLines(list_SymbolID$Symbol, "downstream/Cytoscape/common_genes_STRING.txt")

# 02-Load MCC table 
if(file.exists("downstream/Cytoscape/MCC_Value.csv")) {
  
MCC <- read.csv("downstream/Cytoscape/MCC_Value.csv", header = FALSE)
} else cat("\nFile not found!")

MCC <- MCC[-1,]
colnames(MCC) <- MCC[1, ]
MCC <- MCC[-1, ]

write.csv(MCC, "downstream/Cytoscape/MCC_ValueClean.csv", row.names = FALSE, col.names = TRUE)
