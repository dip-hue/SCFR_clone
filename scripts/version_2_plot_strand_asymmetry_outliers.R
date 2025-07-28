#!/usr/bin/env Rscript

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
  make_option(c("-b", "--bed"), type="character", default="asymmetry_outlier_regions.bed", help="Output BED file"),
  make_option(c("--min_region_size"), type="integer", default=1000000, help="Minimum region size in bp to report")
)

opt <- parse_args(OptionParser(option_list=option_list))

# ---------------------------
# Read data
# ---------------------------
df <- read_csv(opt$input, show_col_types = FALSE)

# ---------------------------
# Identify high/low outliers per chromosome
# ---------------------------
outlier_regions <- list()
outlier_marks <- list()

chromosomes <- unique(df$chrom)

for (chr in chromosomes) {
  chr_df <- df %>% filter(chrom == chr)
  q1 <- quantile(chr_df$strand_count_asym, 0.25, na.rm = TRUE)
  q3 <- quantile(chr_df$strand_count_asym, 0.75, na.rm = TRUE)
  iqr_val <- q3 - q1

  high_cutoff <- q3 + 1.5 * iqr_val
  low_cutoff  <- q1 - 1.5 * iqr_val

  chr_df <- chr_df %>%
    mutate(region_type = case_when(
      strand_count_asym >= high_cutoff ~ "high",
      strand_count_asym <= low_cutoff  ~ "low",
      TRUE ~ NA_character_
    ))

  # Store marked windows
  outlier_marks[[chr]] <- chr_df

  # Identify contiguous regions
  chr_regions <- chr_df %>%
    filter(!is.na(region_type)) %>%
    arrange(window_start) %>%
    mutate(grp = cumsum(lag(region_type, default = first(region_type)) != region_type)) %>%
    group_by(region_type, grp) %>%
    summarise(
      chrom = chr,
      start = min(window_start),
      end = max(window_end),
      length = end - start,
      .groups = "drop"
    ) %>%
    filter(length >= opt$min_region_size) %>%
    select(chrom, start, end, region_type)

  outlier_regions[[chr]] <- chr_regions
}

# ---------------------------
# Write BED file
# ---------------------------
bed_df <- bind_rows(outlier_regions)
write_tsv(bed_df, opt$bed, col_names = FALSE)
cat("Outlier BED file written to:", opt$bed, "\n")

# ---------------------------
# Plot per chromosome
# ---------------------------
pdf(opt$pdf, width = 10, height = 6)

for (chr in chromosomes) {
  chr_df <- outlier_marks[[chr]]
  chr_regions <- outlier_regions[[chr]]

  p <- ggplot(chr_df, aes(x = window_start, y = strand_count_asym)) +
    geom_line(color = "grey30") +
    labs(
      title = paste("Strand Count Asymmetry (Outliers) -", chr),
      x = "Window Start Position (bp)",
      y = "Strand Count Asymmetry"
    ) +
    theme_minimal(base_size = 14)

  if (nrow(chr_regions) > 0) {
    p <- p +
      geom_rect(
        data = chr_regions,
        aes(xmin = start, xmax = end, ymin = -1, ymax = 1, fill = region_type),
        alpha = 0.3,
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c(high = "tomato", low = "steelblue"), name = "Outlier Type")
  }

  print(p)
}

dev.off()
cat("PDF plots written to:", opt$pdf, "\n")
