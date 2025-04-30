library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

# --- 1. Genes of Interest (Tumor Suppressor and Oncogenes) ---

tumor_suppressor_genes <- c("MYH11", "SAMHD1")
oncogenes <- c("ALDH1L2", "CD276", "WWTR1")

# --- 2. Read Gene Occurrence Counts from CSVs ---

oncogene_counts_file_path <- "C:/Users/ns01288/Documents/gene_occurrences_wong_plus_tarney_counts.csv" # Upregulated
tumor_suppressor_counts_file_path <- "C:/Users/ns01288/Documents/gene_occurrences_wong_tarney_downregulated.csv" # Downregulated

oncogene_counts_df <- read_csv(oncogene_counts_file_path)
oncogene_counts_df <- as.data.frame(oncogene_counts_df)

tumor_suppressor_counts_df <- read_csv(tumor_suppressor_counts_file_path)
tumor_suppressor_counts_df <- as.data.frame(tumor_suppressor_counts_df)

# --- 3. Filter for Tumor Suppressor and Oncogene Counts ---

tumor_suppressor_counts <- tumor_suppressor_counts_df %>%
  filter(Gene %in% tumor_suppressor_genes)
oncogene_counts <- oncogene_counts_df %>%
  filter(Gene %in% oncogenes)

# --- 4. Create Bar Graphs ---

tumor_suppressor_plot <- NULL
oncogene_plot <- NULL

if (nrow(tumor_suppressor_counts) > 0) {
  tumor_suppressor_plot <- ggplot(tumor_suppressor_counts, aes(x = Gene, y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", width = 0.4) +
    labs(title = "Tumor Suppressor Gene Counts (Downregulated)", x = "Gene", y = "Count") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme(plot.title = element_text(hjust = 0.5)) # Center title
}

if (nrow(oncogene_counts) > 0) {
  oncogene_plot <- ggplot(oncogene_counts, aes(x = Gene, y = Count)) +
    geom_bar(stat = "identity", fill = "salmon", width = 0.4) +
    labs(title = "Oncogene Counts (Upregulated)", x = "Gene", y = "Count") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme(plot.title = element_text(hjust = 0.5)) # Center title
}

# --- 5. Patchwork of Bar Graphs ---

if (!is.null(tumor_suppressor_plot) || !is.null(oncogene_plot)) {
  if (!is.null(tumor_suppressor_plot) && !is.null(oncogene_plot)) {
    combined_plot <- tumor_suppressor_plot + oncogene_plot
    print(combined_plot)
  } else if (!is.null(tumor_suppressor_plot)) {
    print(tumor_suppressor_plot)
  } else {
    print(oncogene_plot)
  }
} else {
  cat("No graphs to display (no genes present).\n")
}