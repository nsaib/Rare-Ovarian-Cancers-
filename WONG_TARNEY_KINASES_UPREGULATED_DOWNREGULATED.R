library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(patchwork)

# --- 1. Read Upregulated Genes from Wong CSV ---

wong_file_path <- "C:/Users/ns01288/Documents/upregulated_genes_by_condition.csv" # Ensure correct path

wong_data_up <- read.csv(wong_file_path, stringsAsFactors = FALSE)

# Initialize an empty list to store all genes
wong_up_genes <- c()

# Loop through the columns containing gene lists (starting from the 3rd column)
for (col_index in 3:ncol(wong_data_up)) {
  for (row_index in 1:nrow(wong_data_up)) {
    gene_list <- wong_data_up[row_index, col_index]
    
    # Check for empty cells (length zero) and NA
    if (length(gene_list) > 0 && !is.na(gene_list)) {
      if (gene_list != "") {
        genes <- trimws(unlist(strsplit(gene_list, ","))) # Split and trim spaces
        wong_up_genes <- c(wong_up_genes, genes)
      }
    }
  }
}

wong_up_genes <- unique(wong_up_genes) # remove duplicates.

# --- 2. Read Downregulated Genes from Wong CSV ---

wong_file_path <- "C:/Users/ns01288/Documents/downregulated_genes_by_condition.csv" # Ensure correct path

wong_data_down <- read.csv(wong_file_path, stringsAsFactors = FALSE)

# Initialize an empty list to store all genes
wong_down_genes <- c()

# Loop through the columns containing gene lists (starting from the 3rd column)
for (col_index in 3:ncol(wong_data_down)) {
  for (row_index in 1:nrow(wong_data_down)) {
    gene_list <- wong_data_down[row_index, col_index]
    
    # Check for empty cells (length zero) and NA
    if (length(gene_list) > 0 && !is.na(gene_list)) {
      if (gene_list != "") {
        genes <- trimws(unlist(strsplit(gene_list, ","))) # Split and trim spaces
        wong_down_genes <- c(wong_down_genes, genes)
      }
    }
  }
}

wong_down_genes <- unique(wong_down_genes) # remove duplicates.


# --- 3. Read Upregulated Genes from Tarney Sheet ---

tarney_file_path <- "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al.xlsx" # Ensure correct path

Tarney_data <- read_excel(tarney_file_path, sheet = "Tarney_et_al_Total_Proteome")

upregulated_tarney_genes <- Tarney_data %>%
  filter(TP_Tarney_et_al_logFC_LGSOCvsTube > 0.5) %>%
  dplyr::select(GS) %>%
  pull(GS)

# --- 4. Read Downregulated Genes from Tarney Sheet ---

tarney_file_path <- "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al.xlsx" # Ensure correct path

Tarney_data <- read_excel(tarney_file_path, sheet = "Tarney_et_al_Total_Proteome")

downregulated_tarney_genes <- Tarney_data %>%
  filter(TP_Tarney_et_al_logFC_LGSOCvsTube < -0.5) %>%
  dplyr::select(GS) %>%
  pull(GS)

# --- 5. Find Common Genes ---

common_up_genes <- intersect(wong_up_genes, upregulated_tarney_genes)
common_down_genes <- intersect(wong_down_genes, downregulated_tarney_genes)


# --- 6. Read Protein Kinase Genes ---

kinase_file_path <- "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al.xlsx" # Ensure correct path

kinase_data <- read_excel(kinase_file_path, sheet = "Protein_Kinase")

kinase_genes <- kinase_data$`Entrez Gene Symbol`

# --- 7. Find Overlap with Protein Kinases ---

common_up_kinase_genes <- intersect(common_up_genes, kinase_genes)
common_down_kinase_genes <- intersect(common_down_genes, kinase_genes)

# --- 8. Pathway Analysis and Dot Plots ---

# Function to perform pathway analysis and generate dot plots
pathway_analysis_dotplot <- function(gene_list, title) {
  if (length(gene_list) == 0) {
    print(paste("No genes for pathway analysis:", title))
    return(NULL)
  }
  
  gene.df <- bitr(gene_list,
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
  
  if (is.null(gene.df) || nrow(gene.df) == 0) {
    print(paste("No Entrez IDs found for pathway analysis:", title))
    return(NULL)
  }
  
  gene_list_entrez <- gene.df$ENTREZID
  
  ego <- enrichGO(gene          = gene_list_entrez,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  if (!is.null(ego) && nrow(ego) > 0) {
    dotplot(ego, showCategory = 10, title = title) +
      theme(axis.text.y = element_text(size = 8))
  } else {
    print(paste("No significant GO terms found for:", title))
    NULL
  }
}

# Perform pathway analysis and generate dot plots
up_plot <- pathway_analysis_dotplot(common_up_kinase_genes, "Upregulated Common Kinase Genes")
down_plot <- pathway_analysis_dotplot(common_down_kinase_genes, "Downregulated Common Kinase Genes")

# --- 9. Combine Plots with Patchwork ---

if (!is.null(up_plot) && !is.null(down_plot)) {
  combined_plot <- up_plot + down_plot + plot_layout(ncol = 2)
  print(combined_plot)
  ggsave("combined_kinase_pathway_analysis.png", combined_plot, width = 12, height = 6)
} else if (!is.null(up_plot)) {
  print(up_plot)
  ggsave("upregulated_kinase_pathway_analysis.png", up_plot, width = 6, height = 6)
} else if (!is.null(down_plot)) {
  print(down_plot)
  ggsave("downregulated_kinase_pathway_analysis.png", down_plot, width = 6, height = 6)
} else {
  print("No pathway analysis plots to combine.")
}