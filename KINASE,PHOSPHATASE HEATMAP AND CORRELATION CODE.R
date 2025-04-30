##Protein Kinase Heatmap 

library(readxl)
library(dplyr)
library(gplots)

# 1. Load the CSV data
file_path <- "C:/Users/ns01288/Documents/full_joined_data.csv"
full_joined_data <- read.csv(file_path)

# 2. Identify "Condition" columns (using your provided names)
condition_cols <- c(
  "Wong_et_al_TP_LGS125",
  "Wong_et_al_TP_LGS111",
  "Wong_et_al_TP_LGS127",
  "Wong_et_al_TP_LGS128",
  "Wong_et_al_TP_LGS131",
  "Wong_et_al_TP_LGS107",
  "Wong_et_al_TP_LGS108",
  "Wong_et_al_TP_LGS101",
  "Wong_et_al_TP_LGS102",
  "Wong_et_al_TP_LGS124",
  "Wong_et_al_TP_LGS112",
  "Wong_et_al_TP_LGS126",
  "Wong_et_al_TP_LGS129",
  "Wong_et_al_TP_LGS130",
  "TP_Tarney_et_al_logFC_LGSOCvsTube"
)

# 3. Filter for logFC > 0.5 and store filtered data frames
filtered_up_dfs <- lapply(condition_cols, function(col) {
  full_joined_data %>%
    dplyr::filter(.data[[col]] > 0.5) %>%
    dplyr::select(GS, .data[[col]])
})

# 4. Count gene occurrences (logFC > 0.5)
gene_counts_up <- table(unlist(lapply(filtered_up_dfs, function(df) df$GS)))

# 5. Identify genes appearing > 3 times (logFC > 0.5)
frequent_up_genes <- names(gene_counts_up[gene_counts_up > 3])

# 6. Repeat for logFC < -0.5
filtered_down_dfs <- lapply(condition_cols, function(col) {
  full_joined_data %>%
    dplyr::filter(.data[[col]] < -0.5) %>%
    dplyr::select(GS, .data[[col]])
})

gene_counts_down <- table(unlist(lapply(filtered_down_dfs, function(df) df$GS)))
frequent_down_genes <- names(gene_counts_down[gene_counts_down > 3])

# 7. Read Protein Kinase Genes from the Excel file
protein_kinase_data <- read_excel(
  "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al2.xlsx",
  sheet = "Protein_Kinase"
)

# Assuming the protein kinase genes are in the first column
protein_kinase_genes <- unique(protein_kinase_data[[1]])
protein_kinase_genes <- protein_kinase_genes[!is.na(protein_kinase_genes) & protein_kinase_genes != ""] # Remove NAs and empty strings

# 8. Filter Your Gene Lists
up_kinase_genes <- intersect(frequent_up_genes, protein_kinase_genes)
down_kinase_genes <- intersect(frequent_down_genes, protein_kinase_genes)

# 9. Combine and Extract LogFC Data
kinase_genes <- unique(c(up_kinase_genes, down_kinase_genes))

if (length(kinase_genes) > 0) {
  kinase_logfc_matrix <- full_joined_data %>%
    dplyr::filter(GS %in% kinase_genes) %>%
    dplyr::select(GS, all_of(condition_cols)) %>%
    tibble::column_to_rownames("GS") %>%
    as.matrix()
  
  # 10. Generate Heatmap with adjusted parameters
  heatmap.2(kinase_logfc_matrix,
            main = "Protein Kinase LogFC Heatmap",
            trace = "none",
            col = colorRampPalette(c("blue", "white", "red"))(256),
            margin = c(10, 20),
            las = 2,
            cexCol = 0.8,
            cex.main = 0.7,
            lhei = c(1.5, 4),
            lwid = c(1.5, 4),
            labCol = gsub("TP_Tarney_et_al_logFC_LGSOCvsTube", "TP_Tarney_et_al_LGSOCvsTube", colnames(kinase_logfc_matrix)) # Change label for display
  )
} else {
  print("No protein kinase genes found.")
}


###Phosphatase heatmap 
library(readxl)
library(dplyr)
library(gplots)

# 1. Load the CSV data
file_path <- "C:/Users/ns01288/Documents/full_joined_data.csv"
full_joined_data <- read.csv(file_path)

# 2. Identify "Condition" columns (using your provided names)
condition_cols <- c(
  "Wong_et_al_TP_LGS125",
  "Wong_et_al_TP_LGS111",
  "Wong_et_al_TP_LGS127",
  "Wong_et_al_TP_LGS128",
  "Wong_et_al_TP_LGS131",
  "Wong_et_al_TP_LGS107",
  "Wong_et_al_TP_LGS108",
  "Wong_et_al_TP_LGS101",
  "Wong_et_al_TP_LGS102",
  "Wong_et_al_TP_LGS124",
  "Wong_et_al_TP_LGS112",
  "Wong_et_al_TP_LGS126",
  "Wong_et_al_TP_LGS129",
  "Wong_et_al_TP_LGS130",
  "TP_Tarney_et_al_logFC_LGSOCvsTube"
)

# 3. Filter for logFC > 0.5 and store filtered data frames
filtered_up_dfs <- lapply(condition_cols, function(col) {
  full_joined_data %>%
    dplyr::filter(.data[[col]] > 0.5) %>%
    dplyr::select(GS, .data[[col]])
})

# 4. Count gene occurrences (logFC > 0.5)
gene_counts_up <- table(unlist(lapply(filtered_up_dfs, function(df) df$GS)))

# 5. Identify genes appearing > 3 times (logFC > 0.5)
frequent_up_genes <- names(gene_counts_up[gene_counts_up > 3])

# 6. Repeat for logFC < -0.5
filtered_down_dfs <- lapply(condition_cols, function(col) {
  full_joined_data %>%
    dplyr::filter(.data[[col]] < -0.5) %>%
    dplyr::select(GS, .data[[col]])
})

gene_counts_down <- table(unlist(lapply(filtered_down_dfs, function(df) df$GS)))
frequent_down_genes <- names(gene_counts_down[gene_counts_down > 3])

# 7. Read Phosphatase Genes from the Excel file
phosphatase_data <- read_excel(
  "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al2.xlsx",
  sheet = "Phosphatase"
)

# Assuming the phosphatase genes are in the first column
phosphatase_genes <- unique(phosphatase_data[[1]])
phosphatase_genes <- phosphatase_genes[!is.na(phosphatase_genes) & phosphatase_genes != ""] # Remove NAs and empty strings

# 8. Filter Your Gene Lists
up_phosphatase_genes <- intersect(frequent_up_genes, phosphatase_genes)
down_phosphatase_genes <- intersect(frequent_down_genes, phosphatase_genes)

# 9. Combine and Extract LogFC Data
phosphatase_genes_final <- unique(c(up_phosphatase_genes, down_phosphatase_genes))

if (length(phosphatase_genes_final) > 0) {
  phosphatase_logfc_matrix <- full_joined_data %>%
    dplyr::filter(GS %in% phosphatase_genes_final) %>%
    dplyr::select(GS, all_of(condition_cols)) %>%
    tibble::column_to_rownames("GS") %>%
    as.matrix()
  
  # 10. Generate Heatmap with adjusted parameters
  heatmap.2(phosphatase_logfc_matrix,
            main = "Phosphatase LogFC Heatmap",
            trace = "none",
            col = colorRampPalette(c("blue", "white", "red"))(256),
            margin = c(10, 20),
            las = 2,
            cexCol = 0.8,
            cex.main = 0.7,
            lhei = c(1.5, 4),
            lwid = c(1.5, 4),
            labCol = gsub("TP_Tarney_et_al_logFC_LGSOCvsTube", "TP_Tarney_et_al_LGSOCvsTube", colnames(phosphatase_logfc_matrix))
  )
} else {
  print("No phosphatase genes found.")
}

### correlation plot between kinase and phosphatases 

library(readxl)
library(dplyr)
library(corrplot)

# 1. Load the CSV data
file_path <- "C:/Users/ns01288/Documents/full_joined_data.csv"
full_joined_data <- read.csv(file_path)

# 2. Identify "Condition" columns
condition_cols <- c(
  "Wong_et_al_TP_LGS125",
  "Wong_et_al_TP_LGS111",
  "Wong_et_al_TP_LGS127",
  "Wong_et_al_TP_LGS128",
  "Wong_et_al_TP_LGS131",
  "Wong_et_al_TP_LGS107",
  "Wong_et_al_TP_LGS108",
  "Wong_et_al_TP_LGS101",
  "Wong_et_al_TP_LGS102",
  "Wong_et_al_TP_LGS124",
  "Wong_et_al_TP_LGS112",
  "Wong_et_al_TP_LGS126",
  "Wong_et_al_TP_LGS129",
  "Wong_et_al_TP_LGS130",
  "TP_Tarney_et_al_logFC_LGSOCvsTube"
)

# 3. Filter for logFC > 0.5 and store filtered data frames
filtered_up_dfs <- lapply(condition_cols, function(col) {
  full_joined_data %>%
    dplyr::filter(.data[[col]] > 0.5) %>%
    dplyr::select(GS, .data[[col]])
})

# 4. Count gene occurrences (logFC > 0.5)
gene_counts_up <- table(unlist(lapply(filtered_up_dfs, function(df) df$GS)))

# 5. Identify genes appearing > 3 times (logFC > 0.5)
frequent_up_genes <- names(gene_counts_up[gene_counts_up > 3])

# 6. Repeat for logFC < -0.5
filtered_down_dfs <- lapply(condition_cols, function(col) {
  full_joined_data %>%
    dplyr::filter(.data[[col]] < -0.5) %>%
    dplyr::select(GS, .data[[col]])
})

# 7. Count gene occurrences (logFC < -0.5)
gene_counts_down <- table(unlist(lapply(filtered_down_dfs, function(df) df$GS)))

# 8. Identify genes appearing > 3 times (logFC < -0.5)
frequent_down_genes <- names(gene_counts_down[gene_counts_down > 3])

# 9. Read Protein Kinase Genes from the Excel file
protein_kinase_data <- read_excel(
  "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al2.xlsx",
  sheet = "Protein_Kinase"
)

# 10. Extract Protein Kinase genes
protein_kinase_genes <- unique(protein_kinase_data[[1]])
protein_kinase_genes <- protein_kinase_genes[!is.na(protein_kinase_genes) & protein_kinase_genes != ""]

# 11. Read Phosphatase Genes from the Excel file
phosphatase_data <- read_excel(
  "C:/Users/ns01288/OneDrive - University of Surrey/data using disertation/Total_Proteome_wong_et_al_Tarney_et_al2.xlsx",
  sheet = "Phosphatase"
)

# 12. Extract Phosphatase genes
phosphatase_genes <- unique(phosphatase_data[[1]])
phosphatase_genes <- phosphatase_genes[!is.na(phosphatase_genes) & phosphatase_genes != ""]

# 13. Filter for relevant genes (up/down and kinase/phosphatase)
up_kinase_genes <- intersect(frequent_up_genes, protein_kinase_genes)
down_kinase_genes <- intersect(frequent_down_genes, protein_kinase_genes)
up_phosphatase_genes <- intersect(frequent_up_genes, phosphatase_genes)
down_phosphatase_genes <- intersect(frequent_down_genes, phosphatase_genes)

# 14. Combine all relevant genes
all_relevant_genes <- unique(c(up_kinase_genes, down_kinase_genes, up_phosphatase_genes, down_phosphatase_genes))

# 15. Extract LogFC data for all relevant genes
if (length(all_relevant_genes) > 0) {
  all_logfc_matrix <- full_joined_data %>%
    dplyr::filter(GS %in% all_relevant_genes) %>%
    dplyr::select(GS, all_of(condition_cols)) %>%
    tibble::column_to_rownames("GS") %>%
    as.matrix()
  
  # 16. Separate Kinase and Phosphatase Data
  kinase_data <- all_logfc_matrix[rownames(all_logfc_matrix) %in% c(up_kinase_genes, down_kinase_genes), ]
  phosphatase_data <- all_logfc_matrix[rownames(all_logfc_matrix) %in% c(up_phosphatase_genes, down_phosphatase_genes), ]
  
  # 17. Calculate Spearman Correlation
  correlation_matrix <- cor(t(kinase_data), t(phosphatase_data), method = "spearman")
  
  # 18. Filter out genes with completely zero correlations and NA values
  row_zeros <- apply(correlation_matrix, 1, function(row) all(row == 0 | is.na(row)))
  col_zeros <- apply(correlation_matrix, 2, function(col) all(col == 0 | is.na(col)))
  
  correlation_matrix_filtered <- correlation_matrix[!row_zeros, !col_zeros]
  
  # 19. Save the filtered correlation matrix to a CSV file
  output_file_path <- "C:/Users/ns01288/Documents/kinase_phosphatase_correlation_filtered.csv"
  write.csv(correlation_matrix_filtered, file = output_file_path)
  
  print(paste("Filtered correlation matrix saved to:", output_file_path))
  
  # Function to create corrplot with adjustable parameters (PDF output, LANDSCAPE, WITH DENDROGRAMS and CLUSTERS)
  create_corrplot_pdf_landscape_dendro_clusters <- function(corr_matrix, 
                                                            filename = "Kinase_Phosphatase_Correlation_landscape_dendro_clusters.pdf",
                                                            width = 18, height = 9, 
                                                            tl_cex = 1.2, cl_cex = 1.5,
                                                            mar = c(5, 5, 5, 5),
                                                            tl_srt = 45, 
                                                            hclust_method = "ward.D2",
                                                            num_clusters = 4,
                                                            title_cex = 3.5,
                                                            kinase_genes, 
                                                            phosphatase_genes) {
    
    pdf(filename, width = width, height = height)
    
    # Create custom labels with "K" for kinases and "P" for phosphatases
    custom_labels <- c()
    for (gene in rownames(corr_matrix)) {
      if (gene %in% kinase_genes) {
        custom_labels <- c(custom_labels, paste0(gene, " (K)"))
      } else if (gene %in% phosphatase_genes) {
        custom_labels <- c(custom_labels, paste0(gene, " (P)"))
      } else {
        custom_labels <- c(custom_labels, gene)
      }
    }
    
    suppressWarnings({
      corrplot(corr_matrix,
               method = "color",
               type = "full",
               diag = TRUE,
               tl.col = "black",
               tl.cex = tl_cex,
               cl.lim = c(-1, 1),
               title = "Kinase-Phosphatase Correlation",
               col = colorRampPalette(c("blue", "cyan", "white", "yellow", "red"))(200),
               cl.cex = cl_cex,
               mar = mar,
               tl.srt = tl_srt,
               order = "hclust",
               hclust.method = hclust_method,
               addrect = num_clusters,
               title.cex = title_cex,
               tl.labels = custom_labels,
               cl.pos = "b"
      )
      
      # **NO MTEXT HERE** - We'll rely on the gene labels themselves.
      
    })
    dev.off()
  }
  
  # 20. Corrplot the filtered correlation matrix using the function (PDF, LANDSCAPE, WITH DENDROGRAMS and CLUSTERS)
  create_corrplot_pdf_landscape_dendro_clusters(correlation_matrix_filtered,
                                                filename = "Kinase_Phosphatase_Correlation_landscape_dendro_clusters.pdf",
                                                width = 23, height = 10,
                                                tl_cex = 1.3, cl_cex = 1.6,
                                                mar = c(0.5, 0.5, 0.5, 0.0), tl_srt = 45, # Reduced right margin
                                                hclust_method = "ward.D2",
                                                num_clusters = 5,
                                                title_cex = 4.0,
                                                protein_kinase_genes,
                                                phosphatase_genes)
  
} else {
  print("No kinase or phosphatase genes found.")
}