#!/usr/bin/env Rscript

# Load libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(optparse)
})

# ---------------------------
# Command-line arguments
# ---------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character", help="CSV file with strand asymmetry per window"),
  make_option(c("-p", "--pdf"), type="character", default="strand_asymmetry_plots.pdf", help="Output PDF file"),
  make_option(c("-b", "--bed"), type="character", default="high_low_asymmetry_regions.bed", help="Output BED file"),
  make_option(c("--high_thresh"), type="double", default=0.8, help="Threshold for high asymmetry"),
  make_option(c("--low_thresh"), type="double", default=-0.8, help="Threshold for low asymmetry"),
  make_option(c("--min_region_size"), type="integer", default=100000, help="Minimum region size in bp to report in BED and highlight in plots")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---------------------------
# Read data
# ---------------------------
df <- read_csv(opt$input, show_col_types = FALSE)

# ---------------------------
# Classify windows
# ---------------------------
df <- df %>%
  mutate(region_type = case_when(
    strand_count_asym >= opt$high_thresh ~ "high",
    strand_count_asym <= opt$low_thresh ~ "low",
    TRUE ~ NA_character_
  ))

# ---------------------------
# Identify contiguous high/low regions
# ---------------------------
regions <- df %>%
  filter(!is.na(region_type)) %>%
  arrange(chrom, window_start) %>%
  group_by(chrom, region_type, grp = cumsum(lag(region_type, default = first(region_type)) != region_type | lag(chrom, default = first(chrom)) != chrom)) %>%
  summarise(
    start = min(window_start),
    end = max(window_end),
    .groups = "drop"
  ) %>%
  mutate(length = end - start) %>%
  filter(length >= opt$min_region_size)

# Save BED file
write_tsv(regions %>% select(chrom, start, end, region_type), opt$bed, col_names = FALSE)
cat("BED file with filtered regions written to:", opt$bed, "\n")

# ---------------------------
# Generate plots
# ---------------------------
pdf(opt$pdf, width = 10, height = 6)

for (chr in unique(df$chrom)) {
  chr_df <- df %>% filter(chrom == chr)
  chr_regions <- regions %>% filter(chrom == chr)

  p <- ggplot(chr_df, aes(x = window_start, y = strand_count_asym)) +
    geom_line(color = "grey40") +
    labs(
      title = paste("Strand Count Asymmetry -", chr),
      x = "Window Start Position (bp)",
      y = "Strand Count Asymmetry"
    ) +
    theme_minimal(base_size = 14) +
    geom_hline(yintercept = c(opt$high_thresh, opt$low_thresh), linetype = "dashed", color = "red")

  if (nrow(chr_regions) > 0) {
    p <- p + 
      geom_rect(
        data = chr_regions,
        aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = region_type),
        alpha = 0.3,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c(high = "tomato", low = "steelblue"), name = "Region Type")
  }

  print(p)
}

dev.off()
cat("PDF plots written to:", opt$pdf, "\n")
