##Biologial processes

library(ggplot2)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readxl)
library(enrichplot)

# --- 1. Load Data and Process Wong Frequent Genes ---

file_path <- "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al.xlsx"

tryCatch({
  Wong_et_al_Total_Proteome <- read_excel(file_path, sheet = "Wong_et_al_Total_Proteome")
  names(Wong_et_al_Total_Proteome) <- gsub("Vong", "Wong", names(Wong_et_al_Total_Proteome))
  Tarney_data <- read_excel(file_path, sheet = "Tarney_et_al_Total_Proteome")
}, error = function(e) {
  stop(paste("Error loading data:", e))
})

wong_condition_cols <- c(
  "Wong_et_al_TP_LGS125", "Wong_et_al_TP_LGS111", "Wong_et_al_TP_LGS128",
  "Wong_et_al_TP_LGS127", "Wong_et_al_TP_LGS131", "Wong_et_al_TP_LGS107",
  "Wong_et_al_TP_LGS108", "Wong_et_al_TP_LGS101", "Wong_et_al_TP_LGS124",
  "Wong_et_al_TP_LGS102", "Wong_et_al_TP_LGS112", "Wong_et_al_TP_LGS126",
  "Wong_et_al_TP_LGS130", "Wong_et_al_TP_LGS129"
)

logFC_threshold <- 0.5

wong_up_genes <- lapply(wong_condition_cols, function(col) Wong_et_al_Total_Proteome$GS[Wong_et_al_Total_Proteome[[col]] > logFC_threshold])
wong_down_genes <- lapply(wong_condition_cols, function(col) Wong_et_al_Total_Proteome$GS[Wong_et_al_Total_Proteome[[col]] < -logFC_threshold])

tarney_up_genes <- Tarney_data$GS[Tarney_data$TP_Tarney_et_al_logFC_LGSOCvsTube > logFC_threshold]
tarney_down_genes <- Tarney_data$GS[Tarney_data$TP_Tarney_et_al_logFC_LGSOCvsTube < -logFC_threshold]

all_up_genes <- unique(unlist(c(wong_up_genes, list(tarney_up_genes))))
all_down_genes <- unique(unlist(c(wong_down_genes, list(tarney_down_genes))))

# Find frequent Wong genes (upregulated in at least 5 conditions)
wong_all_up_genes <- unlist(wong_up_genes)
wong_up_gene_counts <- table(wong_all_up_genes)
wong_frequent_up_genes <- names(wong_up_gene_counts[wong_up_gene_counts >= 5])

# Find frequent Wong genes (downregulated in at least 5 conditions)
wong_all_down_genes <- unlist(wong_down_genes)
wong_down_gene_counts <- table(wong_all_down_genes)
wong_frequent_down_genes <- names(wong_down_gene_counts[wong_down_gene_counts >= 5])

# --- 2. Find Common Genes and Filter for Frequent Wong Genes ---

common_up_genes <- intersect(tarney_up_genes, unlist(wong_up_genes))
common_down_genes <- intersect(tarney_down_genes, unlist(wong_down_genes))

# Filter for frequent genes
common_up_genes <- intersect(common_up_genes, wong_frequent_up_genes)
common_down_genes <- intersect(common_down_genes, wong_frequent_down_genes)

tarney_up_exclusive <- setdiff(tarney_up_genes, unlist(wong_up_genes))
tarney_down_exclusive <- setdiff(tarney_down_genes, unlist(wong_down_genes))

wong_up_exclusive <- setdiff(unlist(wong_up_genes), tarney_up_genes)
wong_down_exclusive <- setdiff(unlist(wong_down_genes), tarney_down_genes)

# --- 3. Convert Gene Symbols to Entrez IDs ---

tarney_up_ids <- bitr(tarney_up_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
tarney_up_ids <- tarney_up_ids[!is.na(tarney_up_ids$ENTREZID), ]
tarney_down_ids <- bitr(tarney_down_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
tarney_down_ids <- tarney_down_ids[!is.na(tarney_down_ids$ENTREZID), ]
tarney_down_ids <- tarney_down_ids[!duplicated(tarney_down_ids$SYMBOL), ]
wong_up_ids <- bitr(wong_up_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
wong_down_ids <- bitr(wong_down_exclusive, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
common_up_ids <- bitr(common_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
common_down_ids <- bitr(common_down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# --- 4. Create Gene Lists for compareCluster ---

gene.lists <- list(
  Tarney_Up_Exclusive = tarney_up_ids$SYMBOL,
  Tarney_Down_Exclusive = tarney_down_ids$SYMBOL,
  Wong_Up_Exclusive = wong_up_ids$SYMBOL,
  Wong_Down_Exclusive = wong_down_ids$SYMBOL,
  Common_Up = common_up_ids$SYMBOL,
  Common_Down = common_down_ids$SYMBOL
)

# --- 5. Perform Comparative GO Enrichment Analysis ---

ck <- compareCluster(
  geneCluster = gene.lists,
  fun = enrichGO,
  OrgDb = "org.Hs.eg.db",
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  minGSSize = 10,
  pvalueCutoff = 0.005,
  qvalueCutoff = 0.2
)

# --- 6. Emapplot Visualization ---

ck1 <- pairwise_termsim(ck)

# Improved emapplot call
pathway_b <- emapplot(ck1, 
                      showCategory = 10,  # Adjust as needed - reduce if plot is crowded
                      layout = "fr"      # Try different layouts (fr, kk, graphopt)
) +  
  
  ggtitle("GO Biological Process Network - Enriched Terms") + # More specific title
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 8)) # Control label size with theme

print(pathway_b)


# Example using ggrepel (install if needed: install.packages("ggrepel"))
if(require(ggrepel)){
  pathway_b_repel <- emapplot(ck1, 
                              showCategory = 10, 
                              layout = "fr",
                              repel = TRUE) +  # Use ggrepel
    ggtitle("GO Biological Process Network - Enriched Terms (Repel)") +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 8)) # Control label size with theme
  print(pathway_b_repel)
}


# Example adjusting node size by pvalue (requires some data manipulation)
# (This is a more advanced example and might need adaptation)
# ck_df <- as.data.frame(ck) # Convert compareCluster result to data frame
# GO_term_pvalues <- aggregate(p.adjust ~ Description, data = ck_df, FUN = min) # Get min p.value for each GO term

# pathway_b_sized <- emapplot(ck1, 
#                      showCategory = 10, 
#                      layout = "fr") +
#   geom_point(aes(size = -log10(p.adjust))) +  # Size nodes by pvalue
#   ggtitle("GO Biological Process Network - Enriched Terms (Sized by P-value)") +
#   theme(plot.title = element_text(hjust = 0.5))
# print(pathway_b_sized)