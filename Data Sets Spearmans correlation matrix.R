library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)

# (Data loading code - same as before)
# ...

# Calculate the correlation matrix
cor_matrix <- cor(numeric_data, method = "pearson", use = "complete.obs")

# Prepare data for ggplot2
melted_cor_matrix <- melt(cor_matrix)
melted_cor_matrix <- melted_cor_matrix %>%
  filter(as.numeric(Var1) < as.numeric(Var2))

# Create the ggplot2 correlation plot
ggplot(data = melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limit = c(-1, 1), space = "Lab",
                       name = "Pearson\nCorrelation") +
  coord_fixed() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"), # Bold x-axis labels
        axis.text.y = element_text(face = "bold"), # Bold y-axis labels
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold")) + # Bold and center title
  geom_text(aes(label = round(value, 2)), color = "black", size = 2, fontface = "bold") + # Bold numbers
  ggtitle("Correlation matrix of Wong et al. and Tarney et al. Datasets")

# Save the ggplot2 plot
ggsave("Correlation_matrix_ggplot2_bold.png", width = 8, height = 8, dpi = 300)