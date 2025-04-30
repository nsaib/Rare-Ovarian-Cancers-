##This dot plot visualizes the results of your Gene Ontology (GO) enrichment analysis for upregulated and downregulated genes from your combined Wong and Tarney datasets.

library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

# --- 1. Load Data and Perform Full Join ---

file_path <- "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al.xlsx"

# Load Wong data
Wong_et_al_Total_Proteome <- read_excel(file_path, sheet = "Wong_et_al_Total_Proteome")
names(Wong_et_al_Total_Proteome) <- gsub("Vong", "Wong", names(Wong_et_al_Total_Proteome))

# Load Tarney data
Tarney_data <- read_excel(file_path, sheet = "Tarney_et_al_Total_Proteome")

# Full join on 'GS' (gene symbol)
combined_data <- full_join(Wong_et_al_Total_Proteome, Tarney_data, by = "GS")

# --- 2. Calculate Average LogFC for Wong Data ---

wong_condition_cols <- c(
  "Wong_et_al_TP_LGS125", "Wong_et_al_TP_LGS111", "Wong_et_al_TP_LGS128",
  "Wong_et_al_TP_LGS127", "Wong_et_al_TP_LGS131", "Wong_et_al_TP_LGS107",
  "Wong_et_al_TP_LGS108", "Wong_et_al_TP_LGS101", "Wong_et_al_TP_LGS124",
  "Wong_et_al_TP_LGS102", "Wong_et_al_TP_LGS112", "Wong_et_al_TP_LGS126",
  "Wong_et_al_TP_LGS130", "Wong_et_al_TP_LGS129"
)

combined_data$Wong_Avg_LogFC <- rowMeans(combined_data[, wong_condition_cols], na.rm = TRUE)

# --- 3. Identify Up/Down Regulated Genes ---

logFC_threshold <- 0.5 # Adjust as needed

combined_data$Wong_Up <- combined_data$Wong_Avg_LogFC > logFC_threshold
combined_data$Wong_Down <- combined_data$Wong_Avg_LogFC < -logFC_threshold
combined_data$Tarney_Up <- combined_data$TP_Tarney_et_al_logFC_LGSOCvsTube > logFC_threshold
combined_data$Tarney_Down <- combined_data$TP_Tarney_et_al_logFC_LGSOCvsTube < -logFC_threshold

# --- 4. Upregulated Gene Analysis ---

upregulated_genes <- combined_data$GS[combined_data$Wong_Up | combined_data$Tarney_Up]

gene.df_up <- bitr(upregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_up <- gene.df_up$ENTREZID

ego_up <- enrichGO(gene = gene_list_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

# --- 5. Downregulated Gene Analysis ---

downregulated_genes <- combined_data$GS[combined_data$Wong_Down | combined_data$Tarney_Down]

gene.df_down <- bitr(downregulated_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
gene_list_down <- gene.df_down$ENTREZID

ego_down <- enrichGO(gene = gene_list_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)

# --- 6. Create Separate Dot Plots ---

plot_up <- if (!is.null(ego_up)) {
  dotplot(ego_up, showCategory = 20, title = "Upregulated Genes") +
    theme(axis.text.y = element_text(size = 7))
} else {
  ggplot() + ggtitle("Upregulated Genes: No significant results")
}

plot_down <- if (!is.null(ego_down)) {
  dotplot(ego_down, showCategory = 20, title = "Downregulated Genes") +
    theme(axis.text.y = element_text(size = 7))
} else {
  ggplot() + ggtitle("Downregulated Genes: No significant results")
}

# --- 7. Combine Plots Using patchwork ---

combined_plot <- plot_up + plot_down
print(combined_plot)