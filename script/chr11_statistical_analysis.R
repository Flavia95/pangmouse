#!/usr/bin/env Rscript
# Statistical analysis for chr11 pangenome nodes
# Companion script to analyze_chr11_nodes.py
# Performs advanced statistical tests and creates publication-quality figures

library(tidyverse)
library(ggplot2)
library(scales)
library(viridis)
library(gridExtra)

# Mouse chr11 reference data (mm10)
CHR11_LENGTH <- 121843856

# HOX clusters on chr11
HOX_CLUSTERS <- data.frame(
  name = c("Hoxa"),
  start = c(96161468),
  end = c(96273464),
  stringsAsFactors = FALSE
)

# Centromeric region
CENTROMERE <- list(start = 0, end = 3000000)

# Known duplicated regions
DUPLICATED_REGIONS <- data.frame(
  name = c("Ig_lambda", "Olfactory_1", "Olfactory_2"),
  start = c(113500000, 49000000, 58000000),
  end = c(115000000, 52000000, 60000000),
  stringsAsFactors = FALSE
)

#' Load nodes data
#' @param filename Path to the nodes file (can be gzipped)
load_nodes <- function(filename) {
  cat("Loading data from", filename, "...\n")

  if (grepl("\\.gz$", filename)) {
    df <- read.csv(gzfile(filename), stringsAsFactors = FALSE)
  } else {
    df <- read.csv(filename, stringsAsFactors = FALSE)
  }

  cat("Total nodes loaded:", nrow(df), "\n")
  return(df)
}

#' Filter for chr11 nodes
filter_chr11 <- function(df) {
  chr11_df <- df %>% filter(chr == "chr11")
  cat("Chr11 nodes:", nrow(chr11_df), "\n")
  return(chr11_df)
}

#' Statistical comparison of node types
analyze_node_statistics <- function(chr11_df) {
  cat("\n", rep("=", 60), "\n", sep="")
  cat("STATISTICAL ANALYSIS OF NODE TYPES\n")
  cat(rep("=", 60), "\n", sep="")

  # Summary statistics by tag
  stats <- chr11_df %>%
    group_by(tag) %>%
    summarise(
      count = n(),
      mean_length = mean(length, na.rm = TRUE),
      median_length = median(length, na.rm = TRUE),
      sd_length = sd(length, na.rm = TRUE),
      min_length = min(length, na.rm = TRUE),
      max_length = max(length, na.rm = TRUE),
      total_length_mb = sum(length, na.rm = TRUE) / 1e6,
      .groups = "drop"
    )

  print(stats)

  # Kruskal-Wallis test for length differences between types
  if (length(unique(chr11_df$tag)) > 1) {
    kw_test <- kruskal.test(length ~ tag, data = chr11_df)
    cat("\nKruskal-Wallis test for length differences between node types:\n")
    cat("  Chi-squared =", kw_test$statistic, "\n")
    cat("  df =", kw_test$parameter, "\n")
    cat("  p-value =", format.pval(kw_test$p.value), "\n")

    if (kw_test$p.value < 0.05) {
      cat("  Result: Significant differences in node lengths between types\n")
    }
  }

  return(stats)
}

#' Analyze clustering of duplicates along chromosome
analyze_duplicate_clustering <- function(chr11_df, window_size = 1e6) {
  cat("\n", rep("=", 60), "\n", sep="")
  cat("DUPLICATE CLUSTERING ANALYSIS\n")
  cat(rep("=", 60), "\n", sep="")

  # Filter for nodes with position and duplicates
  pos_df <- chr11_df %>% filter(!is.na(pos))

  if (nrow(pos_df) == 0) {
    cat("No positional information available\n")
    return(NULL)
  }

  # Create bins
  pos_df <- pos_df %>%
    mutate(bin = cut(pos,
                     breaks = seq(0, CHR11_LENGTH, by = window_size),
                     labels = seq(0, CHR11_LENGTH - window_size, by = window_size) / 1e6,
                     include.lowest = TRUE))

  # Count by bin and tag
  bin_counts <- pos_df %>%
    group_by(bin, tag) %>%
    summarise(count = n(), total_length = sum(length) / 1e6, .groups = "drop") %>%
    pivot_wider(names_from = tag, values_from = c(count, total_length), values_fill = 0)

  # Calculate duplicate enrichment ratio
  if ("count_DUPLICATED" %in% colnames(bin_counts)) {
    bin_counts <- bin_counts %>%
      mutate(
        total_count = count_DUPLICATED +
                     ifelse("count_UNIQUE" %in% colnames(bin_counts), count_UNIQUE, 0) +
                     ifelse("count_UNPLACED" %in% colnames(bin_counts), count_UNPLACED, 0),
        dup_ratio = count_DUPLICATED / total_count
      )

    # Find hotspots (bins with >20% duplicates)
    hotspots <- bin_counts %>%
      filter(dup_ratio > 0.2) %>%
      arrange(desc(dup_ratio))

    if (nrow(hotspots) > 0) {
      cat("\nDuplication hotspots (>20% duplicated nodes):\n")
      for (i in 1:min(10, nrow(hotspots))) {
        cat(sprintf("  %s-%s Mb: %.1f%% duplicated (%d nodes)\n",
                   hotspots$bin[i],
                   as.numeric(hotspots$bin[i]) + (window_size/1e6),
                   hotspots$dup_ratio[i] * 100,
                   hotspots$count_DUPLICATED[i]))
      }
    }
  }

  return(bin_counts)
}

#' Create comprehensive visualization
create_chr11_plot <- function(chr11_df, bin_counts, output_file = "chr11_detailed_analysis.png") {
  cat("\nCreating detailed visualization...\n")

  # Set theme
  theme_set(theme_bw(base_size = 12))

  # Prepare data with positions
  pos_df <- chr11_df %>% filter(!is.na(pos))

  if (nrow(pos_df) == 0) {
    cat("Cannot create positional plots without position data\n")
    return(NULL)
  }

  # Plot 1: Node density along chromosome
  p1 <- ggplot(pos_df, aes(x = pos / 1e6, fill = tag)) +
    geom_histogram(bins = 100, alpha = 0.7, position = "stack") +
    labs(x = "Position on Chr11 (Mb)",
         y = "Node Count",
         title = "Node Density Along Chr11") +
    scale_fill_manual(values = c("DUPLICATED" = "#E15759",
                                  "UNIQUE" = "#4E79A7",
                                  "UNPLACED" = "#BAB0AC")) +
    theme(legend.position = "bottom") +
    # Add vertical lines for HOX
    geom_vline(data = HOX_CLUSTERS,
               aes(xintercept = start / 1e6),
               color = "red", linetype = "dashed", alpha = 0.7) +
    geom_vline(data = HOX_CLUSTERS,
               aes(xintercept = end / 1e6),
               color = "red", linetype = "dashed", alpha = 0.7) +
    # Add shading for centromere
    annotate("rect", xmin = CENTROMERE$start / 1e6, xmax = CENTROMERE$end / 1e6,
             ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "gray")

  # Plot 2: Duplication ratio along chromosome
  if (!is.null(bin_counts) && "dup_ratio" %in% colnames(bin_counts)) {
    p2 <- ggplot(bin_counts, aes(x = as.numeric(as.character(bin)), y = dup_ratio * 100)) +
      geom_line(color = "#E15759", size = 1) +
      geom_point(color = "#E15759", size = 2) +
      labs(x = "Position on Chr11 (Mb)",
           y = "% Duplicated Nodes",
           title = "Duplication Rate Along Chr11") +
      geom_hline(yintercept = 20, linetype = "dashed", color = "gray50") +
      # Add regions of interest
      geom_rect(data = DUPLICATED_REGIONS,
                aes(xmin = start / 1e6, xmax = end / 1e6, ymin = -Inf, ymax = Inf),
                alpha = 0.1, fill = "orange", inherit.aes = FALSE)
  } else {
    p2 <- ggplot() +
      annotate("text", x = 60, y = 50, label = "No duplication data available") +
      theme_void()
  }

  # Plot 3: Node length distribution by type
  p3 <- ggplot(chr11_df, aes(x = length, fill = tag)) +
    geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
    labs(x = "Node Length (bp)",
         y = "Count",
         title = "Node Length Distribution by Type") +
    scale_x_log10(labels = comma) +
    scale_y_log10() +
    scale_fill_manual(values = c("DUPLICATED" = "#E15759",
                                  "UNIQUE" = "#4E79A7",
                                  "UNPLACED" = "#BAB0AC")) +
    theme(legend.position = "bottom")

  # Plot 4: Cumulative sequence by type
  seq_summary <- chr11_df %>%
    group_by(tag) %>%
    summarise(total_mb = sum(length) / 1e6, .groups = "drop") %>%
    mutate(tag = factor(tag, levels = c("UNIQUE", "DUPLICATED", "UNPLACED")))

  p4 <- ggplot(seq_summary, aes(x = tag, y = total_mb, fill = tag)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%.1f Mb", total_mb)),
              vjust = -0.5, size = 4) +
    labs(x = "Node Type",
         y = "Total Sequence (Mb)",
         title = "Cumulative Sequence by Type") +
    scale_fill_manual(values = c("DUPLICATED" = "#E15759",
                                  "UNIQUE" = "#4E79A7",
                                  "UNPLACED" = "#BAB0AC")) +
    theme(legend.position = "none")

  # Combine plots
  combined <- grid.arrange(p1, p2, p3, p4, ncol = 2)

  # Save
  ggsave(output_file, combined, width = 16, height = 12, dpi = 300)
  cat("Saved:", output_file, "\n")

  return(combined)
}

#' Main analysis function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    cat("Usage: Rscript chr11_statistical_analysis.R <DBmm10.nodes.tag.txt.gz>\n")
    quit(status = 1)
  }

  input_file <- args[1]

  cat(rep("=", 80), "\n", sep="")
  cat("CHR11 STATISTICAL ANALYSIS\n")
  cat(rep("=", 80), "\n", sep="")

  # Load and filter data
  df <- load_nodes(input_file)
  chr11_df <- filter_chr11(df)

  if (nrow(chr11_df) == 0) {
    cat("ERROR: No chr11 nodes found!\n")
    quit(status = 1)
  }

  # Statistical analyses
  stats <- analyze_node_statistics(chr11_df)
  bin_counts <- analyze_duplicate_clustering(chr11_df)

  # Create visualizations
  create_chr11_plot(chr11_df, bin_counts)

  cat("\n", rep("=", 80), "\n", sep="")
  cat("ANALYSIS COMPLETE\n")
  cat(rep("=", 80), "\n", sep="")
}

# Run main function
if (!interactive()) {
  main()
}
