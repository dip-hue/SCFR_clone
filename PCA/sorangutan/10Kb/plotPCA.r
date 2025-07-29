# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(RColorBrewer)

# Read RSCU data
rscu <- read.table("rscu.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Pivot to wide format
rscu_wide <- rscu %>%
  select(clade, species, codon, rscu) %>%
  pivot_wider(names_from = codon, values_from = rscu, values_fill = 0)

# Extract metadata and numeric matrix
metadata <- rscu_wide %>% select(species, clade)
rscu_matrix <- rscu_wide %>% select(-species, -clade)

# Remove columns with zero variance
rscu_matrix_filtered <- rscu_matrix[, apply(rscu_matrix, 2, var) != 0]

# PCA
pca <- prcomp(rscu_matrix_filtered, scale. = TRUE)
expl <- summary(pca)$importance[2, ]
pca_df <- as.data.frame(pca$x)
pca_df$species <- metadata$species
pca_df$clade <- metadata$clade

# Loadings
loadings <- as.data.frame(pca$rotation)
loadings$codon <- rownames(loadings)

# Define color palette
clades <- unique(pca_df$clade)
n_clades <- length(clades)
palette_colors <- setNames(brewer.pal(max(n_clades, 15), "Set3")[1:n_clades], clades)

# Function to generate individual PCA plots and save as high-res PDF
plot_and_save_pca <- function(df, loadings, x, y, outname, top_n = 5) {
  pcx <- as.numeric(sub("PC", "", x))
  pcy <- as.numeric(sub("PC", "", y))

  # Codon loadings
  loading_df <- loadings %>%
    select(codon, !!x, !!y) %>%
    rename(xload = !!x, yload = !!y) %>%
    mutate(magnitude = xload^2 + yload^2) %>%
    arrange(desc(magnitude)) %>%
    slice_head(n = top_n) %>%
    mutate(xend = xload * 5, yend = yload * 5)

  # Calculate axis expansion
  xlim <- range(df[[x]]) * 1.1
  ylim <- range(df[[y]]) * 1.1

  # PCA plot
  p <- ggplot(df, aes_string(x = x, y = y, color = "clade", label = "species")) +
    geom_point(size = 3, alpha = 0.9) +
    geom_text_repel(size = 3, max.overlaps = 100) +
    scale_color_manual(values = palette_colors) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(x, "vs", y),
      x = paste0(x, " (", round(expl[pcx]*100, 1), "%)"),
      y = paste0(y, " (", round(expl[pcy]*100, 1), "%)"),
      color = "Clade"
    ) +
    geom_segment(
      data = loading_df,
      aes(x = 0, y = 0, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.2, "cm")),
      inherit.aes = FALSE,
      color = "gray30"
    ) +
    geom_text_repel(
      data = loading_df,
      aes(x = xend, y = yend, label = codon),
      size = 3,
      inherit.aes = FALSE,
      max.overlaps = 100
    ) +
    xlim(xlim) + ylim(ylim)

  # Save high-resolution PDF
  ggsave(outname, plot = p, width = 12, height = 10, dpi = 600)
}

# List of all PC pair combinations up to PC5
pcs <- paste0("PC", 1:5)
combinations <- combn(pcs, 2, simplify = FALSE)

# Generate and save plots
for (pair in combinations) {
  x <- pair[1]
  y <- pair[2]
  filename <- paste0("PCA_", x, "_", y, ".pdf")
  plot_and_save_pca(pca_df, loadings, x, y, filename, top_n = 5)
}
