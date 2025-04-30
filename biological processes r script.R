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

wong_up_genes <- lapply(wong_condition_cols, function(col) {
  indices <- which(Wong_et_al_Total_Proteome[[col]] > logFC_threshold)
  if (length(indices) > 0) {
    return(Wong_et_al_Total_Proteome$GS[indices])
  } else {
    return(character(0)) # Return an empty character vector if no genes meet the criteria
  }
})

wong_down_genes <- lapply(wong_condition_cols, function(col) {
  indices <- which(Wong_et_al_Total_Proteome[[col]] < -logFC_threshold)
  if (length(indices) > 0) {
    return(Wong_et_al_Total_Proteome$GS[indices])
  } else {
    return(character(0)) # Return an empty character vector if no genes meet the criteria
  }
})

# ... (rest of your code) ...

# --- 6. Emapplot Visualization (Reduced Terms) ---

ck1 <- pairwise_termsim(ck)

pathway_b <- emapplot(ck1, showCategory = 5) # Show 5 terms (adjust as needed)

pathway_b <- pathway_b +
  ggtitle("Biological process networks (Top 5)") +
  theme(plot.title = element_text(hjust = 0.5))

print(pathway_b)