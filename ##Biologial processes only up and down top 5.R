##Biologial processes only up and down top 5 

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

# --- 3. Convert Gene Symbols to Entrez IDs ---

common_up_ids <- bitr(common_up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
common_down_ids <- bitr(common_down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# --- 4. Create Gene Lists for compareCluster ---

gene.lists <- list(
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

# --- 6. Emapplot Visualization for Common Genes ---

# Convert compareCluster result to data frame
ck_df <- as.data.frame(ck)

# Filter for Common_Up and Common_Down
ck_common_df <- ck_df[ck_df$Cluster %in% c("Common_Up", "Common_Down"), ]

# Calculate pairwise term similarity
ck1_common <- pairwise_termsim(ck)

# Create emapplot
pathway_common <- emapplot(ck1_common,
                           showCategory = 5,
                           layout = "fr"
                           #node_size = 10 # Removed node_size
) +
  ggtitle("GO Biological Process Network - Top 5 Enriched Terms (Common Genes)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 8)
  )

print(pathway_common)

# Save the plot
ggsave("GO_network_common.png", plot = pathway_common, width = 10, height = 8, units = "in", dpi = 300)