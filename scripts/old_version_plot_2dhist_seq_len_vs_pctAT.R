#!/usr/bin/env Rscript

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(tools)
  library(ggrepel)
})

# Command line args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1 || length(args) > 2) {
  stop("Usage: Rscript plot_2dhist_seq_len_vs_pctAT.R <main_file.txt> [<genes_file.txt>]")
}

main_file <- args[1]
genes_file <- ifelse(length(args) == 2, args[2], NA)
output_prefix <- file_path_sans_ext(basename(main_file))
output_png <- paste0(output_prefix, "_2Dhist_seqLen_vs_pctAT.png")

# Read main data
df <- read.table(main_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(df) <- c("col1", "col2", "col3", "col4", "pct_at", "pct_gc",
                  "num_A", "num_C", "num_G", "num_T", "num_N", "num_oth", "seq_len")
df$pct_at <- as.numeric(df$pct_at)
df$seq_len <- as.numeric(df$seq_len)

# Optional gene label data
gene_df <- NULL
if (!is.na(genes_file)) {
  gene_df <- read.table(genes_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_df) <- c("pct_at", "seq_len", "gene")
  gene_df$pct_at <- as.numeric(gene_df$pct_at)
  gene_df$seq_len <- as.numeric(gene_df$seq_len)
}

# Define number of bins
num_bins <- 100

# Define bin edges
at_breaks <- seq(min(df$pct_at, na.rm = TRUE), max(df$pct_at, na.rm = TRUE), length.out = num_bins + 1)
len_breaks <- seq(min(df$seq_len, na.rm = TRUE), max(df$seq_len, na.rm = TRUE), length.out = num_bins + 1)

# Bin the data
df$at_bin <- cut(df$pct_at, breaks = at_breaks, include.lowest = TRUE, right = FALSE)
df$len_bin <- cut(df$seq_len, breaks = len_breaks, include.lowest = TRUE, right = FALSE)

# Count entries per bin
bin_counts <- aggregate(seq_len ~ at_bin + len_bin, data = df, FUN = length)
colnames(bin_counts)[3] <- "count"

# If gene_df is available, bin those too
if (!is.null(gene_df)) {
  gene_df$at_bin <- cut(gene_df$pct_at, breaks = at_breaks, include.lowest = TRUE, right = FALSE)
  gene_df$len_bin <- cut(gene_df$seq_len, breaks = len_breaks, include.lowest = TRUE, right = FALSE)

  # Collapse gene names per bin
  gene_bins <- aggregate(gene ~ at_bin + len_bin, data = gene_df, FUN = function(x) paste(unique(x), collapse = ","))
  
  # Merge with bin_counts
  bin_counts <- merge(bin_counts, gene_bins, by = c("at_bin", "len_bin"), all.x = TRUE)
} else {
  bin_counts$gene <- NA
}

# Write to output file
output_txt <- paste0(output_prefix, "_bin_counts_with_genes.txt")
write.table(bin_counts, file = output_txt, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Bin count summary written to:", output_txt, "\n")

# Create plot
p <- ggplot(df, aes(x = pct_at, y = seq_len)) +
  geom_bin2d(bins = 100) +
  scale_fill_gradient(name = "Count", low = "red", high = "blue") +
  labs(
    title = "2D Histogram: Sequence Length vs %AT Content",
    x = "%AT content",
    y = "Sequence Length"
  ) +
  theme_bw(base_size = 14)

# Add gene labels if present
if (!is.null(gene_df)) {
  p <- p + 
    geom_text_repel(
      data = gene_df,
      aes(x = pct_at, y = seq_len, label = gene),
      size = 3.5,
      color = "black",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "gray50",
      min.segment.length = 0,
      inherit.aes = FALSE
    )
}

# Save plot
ggsave(output_png, plot = p, width = 10, height = 8, dpi = 300)

cat("Plot saved to:", output_png, "\n")

