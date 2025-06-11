#!/usr/bin/env Rscript

# Advanced visualization of variable importance for viral metagenomic assembly strategy
# Creates comprehensive plots to guide co-assembly decisions

library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(corrplot)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

# Function to load PERMANOVA results
load_permanova_results <- function(results_file) {
  cat("Loading PERMANOVA results from:", results_file, "\n")
  
  if (!file.exists(results_file)) {
    stop("PERMANOVA results file not found: ", results_file)
  }
  
  results <- read.csv(results_file, stringsAsFactors = FALSE)
  
  # Add derived variables
  results$neg_log_p <- -log10(results$p_value)
  results$effect_size_category <- cut(
    results$R_squared,
    breaks = c(0, 0.05, 0.1, 0.2, 1),
    labels = c("Small (<0.05)", "Medium (0.05-0.1)", "Large (0.1-0.2)", "Very Large (>0.2)"),
    include.lowest = TRUE
  )
  
  results$significance_category <- ifelse(
    results$p_value <= 0.001, "p ≤ 0.001",
    ifelse(results$p_value <= 0.01, "p ≤ 0.01",
           ifelse(results$p_value <= 0.05, "p ≤ 0.05",
                  ifelse(results$p_value <= 0.1, "p ≤ 0.1", "ns")))
  )
  
  cat("Loaded", nrow(results), "variables\n")
  return(results)
}

# Function to load similarity matrix and metadata for additional plots
load_additional_data <- function(similarity_file, metadata_file) {
  similarity_matrix <- NULL
  metadata <- NULL
  
  if (file.exists(similarity_file)) {
    similarity_matrix <- read.csv(similarity_file, row.names = 1, check.names = FALSE)
    cat("Loaded similarity matrix with", nrow(similarity_matrix), "samples\n")
  }
  
  if (file.exists(metadata_file)) {
    metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
    cat("Loaded metadata with", nrow(metadata), "samples\n")
  }
  
  return(list(similarity_matrix = similarity_matrix, metadata = metadata))
}

# Function to create variable importance ranking plot
create_importance_ranking <- function(results, output_dir) {
  cat("Creating variable importance ranking plot...\n")
  
  # Prepare data
  plot_data <- results %>%
    arrange(desc(R_squared)) %>%
    mutate(Variable = factor(Variable, levels = Variable))
  
  # Create main ranking plot
  ranking_plot <- ggplot(plot_data, aes(x = Variable, y = R_squared)) +
    geom_col(aes(fill = effect_size_category), alpha = 0.8, width = 0.7) +
    geom_text(aes(label = sprintf("%.3f%s", R_squared, Significance)), 
              hjust = -0.1, size = 3, fontface = "bold") +
    coord_flip() +
    labs(
      title = "Variable Importance Ranking",
      subtitle = "Proportion of variance in viral community structure explained by each variable",
      x = "Variables",
      y = "R² (Proportion of variance explained)",
      fill = "Effect Size",
      caption = "Numbers show R² values with significance: *** p≤0.001, ** p≤0.01, * p≤0.05, . p≤0.1"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_fill_viridis_d(option = "plasma", direction = -1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  
  ranking_file <- file.path(output_dir, "variable_importance_ranking.png")
  ggsave(ranking_file, ranking_plot, width = 12, height = 8, dpi = 300)
  cat("Ranking plot saved to:", ranking_file, "\n")
}

# Function to create effect size vs significance bubble plot
create_bubble_plot <- function(results, output_dir) {
  cat("Creating effect size vs significance bubble plot...\n")
  
  bubble_plot <- ggplot(results, aes(x = R_squared, y = neg_log_p)) +
    geom_point(aes(size = F_statistic, color = Variable_Type), alpha = 0.7) +
    geom_text_repel(aes(label = Variable), size = 3, max.overlaps = 15) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.01), color = "darkred", linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = 0.1, color = "blue", linetype = "dashed", alpha = 0.7) +
    geom_vline(xintercept = 0.05, color = "lightblue", linetype = "dashed", alpha = 0.7) +
    labs(
      title = "Variable Effect Size vs Statistical Significance",
      subtitle = "Bubble size represents F-statistic magnitude",
      x = "R² (Effect Size)",
      y = "-log₁₀(p-value)",
      color = "Variable Type",
      size = "F-statistic",
      caption = "Dashed lines: red p=0.05, dark red p=0.01, blue R²=0.1, light blue R²=0.05"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    ) +
    scale_color_brewer(type = "qual", palette = "Set2") +
    scale_size_continuous(range = c(3, 12)) +
    guides(
      color = guide_legend(override.aes = list(size = 5)),
      size = guide_legend(title = "F-statistic")
    )
  
  bubble_file <- file.path(output_dir, "effect_size_significance_bubble.png")
  ggsave(bubble_file, bubble_plot, width = 12, height = 8, dpi = 300)
  cat("Bubble plot saved to:", bubble_file, "\n")
}

# Function to create variable type comparison
create_variable_type_comparison <- function(results, output_dir) {
  cat("Creating variable type comparison plot...\n")
  
  # Summary by variable type
  type_summary <- results %>%
    group_by(Variable_Type) %>%
    summarise(
      n_variables = n(),
      mean_r_squared = mean(R_squared),
      median_r_squared = median(R_squared),
      n_significant = sum(p_value <= 0.05),
      prop_significant = n_significant / n_variables,
      .groups = 'drop'
    )
  
  # Box plot comparing variable types
  type_plot1 <- ggplot(results, aes(x = Variable_Type, y = R_squared)) +
    geom_boxplot(aes(fill = Variable_Type), alpha = 0.7, width = 0.6) +
    geom_jitter(aes(color = significance_category), width = 0.2, size = 2, alpha = 0.8) +
    labs(
      title = "R² Distribution by Variable Type",
      x = "Variable Type",
      y = "R² (Effect Size)",
      fill = "Variable Type",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "right"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    scale_color_manual(values = c("p ≤ 0.001" = "darkred", "p ≤ 0.01" = "red", 
                                  "p ≤ 0.05" = "orange", "p ≤ 0.1" = "yellow", "ns" = "grey"))
  
  # Bar plot of significance proportions
  type_plot2 <- ggplot(type_summary, aes(x = Variable_Type, y = prop_significant)) +
    geom_col(aes(fill = Variable_Type), alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", prop_significant*100, n_significant, n_variables)),
              vjust = -0.1, size = 3, fontface = "bold") +
    labs(
      title = "Proportion of Significant Variables by Type",
      x = "Variable Type",
      y = "Proportion Significant (p ≤ 0.05)",
      fill = "Variable Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "none"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.15)))
  
  # Combine plots
  combined_type_plot <- ggarrange(type_plot1, type_plot2, ncol = 2, nrow = 1)
  
  type_file <- file.path(output_dir, "variable_type_comparison.png")
  ggsave(type_file, combined_type_plot, width = 14, height = 6, dpi = 300)
  cat("Variable type comparison saved to:", type_file, "\n")
}

# Function to create assembly strategy recommendations plot
create_strategy_recommendations <- function(results, output_dir) {
  cat("Creating assembly strategy recommendations plot...\n")
  
  # Identify different variable categories for recommendations
  results$recommendation <- case_when(
    results$R_squared >= 0.1 & results$p_value <= 0.05 ~ "Primary grouping variable",
    results$R_squared >= 0.05 & results$p_value <= 0.1 ~ "Secondary grouping variable",
    results$R_squared >= 0.02 & results$p_value <= 0.05 ~ "Consider for stratification",
    TRUE ~ "Use k-mer similarity only"
  )
  
  # Count recommendations
  rec_summary <- results %>%
    count(recommendation) %>%
    mutate(recommendation = factor(recommendation, levels = c(
      "Primary grouping variable",
      "Secondary grouping variable", 
      "Consider for stratification",
      "Use k-mer similarity only"
    )))
  
  # Create recommendation plot
  rec_plot <- ggplot(results, aes(x = R_squared, y = neg_log_p)) +
    geom_point(aes(color = recommendation, size = F_statistic), alpha = 0.8) +
    geom_text_repel(
      data = filter(results, R_squared >= 0.05 | neg_log_p >= -log10(0.1)),
      aes(label = Variable), size = 3, max.overlaps = 10
    ) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.1), color = "orange", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0.1, color = "blue", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0.05, color = "lightblue", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0.02, color = "grey", linetype = "dashed", alpha = 0.5) +
    labs(
      title = "Assembly Strategy Recommendations Based on Variable Importance",
      subtitle = "Variables are categorized by their effect size and statistical significance",
      x = "R² (Proportion of variance explained)",
      y = "-log₁₀(p-value)",
      color = "Recommendation",
      size = "F-statistic",
      caption = "Dashed lines indicate thresholds for different recommendation categories"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    ) +
    scale_color_manual(values = c(
      "Primary grouping variable" = "darkgreen",
      "Secondary grouping variable" = "forestgreen",
      "Consider for stratification" = "orange",
      "Use k-mer similarity only" = "grey60"
    )) +
    scale_size_continuous(range = c(2, 8)) +
    guides(
      color = guide_legend(override.aes = list(size = 4), ncol = 2),
      size = guide_legend(title = "F-statistic")
    )
  
  strategy_file <- file.path(output_dir, "assembly_strategy_recommendations.png")
  ggsave(strategy_file, rec_plot, width = 14, height = 10, dpi = 300)
  cat("Strategy recommendations plot saved to:", strategy_file, "\n")
  
  # Create summary table plot
  rec_table <- results %>%
    filter(recommendation != "Use k-mer similarity only") %>%
    select(Variable, Variable_Type, R_squared, p_value, recommendation) %>%
    arrange(desc(R_squared))
  
  if (nrow(rec_table) > 0) {
    # Create a text-based table plot
    table_plot <- ggplot() +
      annotation_custom(
        tableGrob(rec_table, rows = NULL, 
                 theme = ttheme_minimal(base_size = 10)),
        xmin = 0, xmax = 1, ymin = 0, ymax = 1
      ) +
      labs(title = "Important Variables for Assembly Strategy") +
      theme_void() +
      theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    
    table_file <- file.path(output_dir, "important_variables_table.png")
    ggsave(table_file, table_plot, width = 12, height = max(6, nrow(rec_table) * 0.3 + 2), dpi = 300)
    cat("Important variables table saved to:", table_file, "\n")
  }
}

# Function to create comprehensive summary dashboard
create_summary_dashboard <- function(results, output_dir) {
  cat("Creating comprehensive summary dashboard...\n")
  
  # Key metrics
  n_total <- nrow(results)
  n_significant <- sum(results$p_value <= 0.05)
  n_large_effect <- sum(results$R_squared >= 0.1)
  n_actionable <- sum(results$R_squared >= 0.1 & results$p_value <= 0.05)
  
  # Create summary text plot
  summary_text <- sprintf(
    "VARIABLE IMPORTANCE ANALYSIS SUMMARY
    
Total Variables Analyzed: %d
Statistically Significant (p ≤ 0.05): %d (%.1f%%)
Large Effect Size (R² ≥ 0.1): %d (%.1f%%)
Actionable Variables (R² ≥ 0.1 & p ≤ 0.05): %d (%.1f%%)

ASSEMBLY STRATEGY RECOMMENDATIONS:
%s",
    n_total,
    n_significant, n_significant/n_total*100,
    n_large_effect, n_large_effect/n_total*100,
    n_actionable, n_actionable/n_total*100,
    ifelse(n_actionable > 0,
           paste("Use", n_actionable, "important variable(s) to guide co-assembly grouping"),
           "Use k-mer similarity as primary grouping criterion")
  )
  
  # Simple text plot
  text_plot <- ggplot() +
    annotate("text", x = 0.05, y = 0.5, label = summary_text, 
             hjust = 0, vjust = 0.5, size = 4, family = "mono") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = "black", size = 1))
  
  summary_file <- file.path(output_dir, "analysis_summary.png")
  ggsave(summary_file, text_plot, width = 10, height = 6, dpi = 300)
  cat("Analysis summary saved to:", summary_file, "\n")
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    results_file <- "results/variable_analysis/permanova_results.csv"
  } else {
    results_file <- args[1]
  }
  
  if (length(args) < 2) {
    output_dir <- dirname(results_file)
  } else {
    output_dir <- args[2]
  }
  
  cat("PERMANOVA results file:", results_file, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load results
  results <- load_permanova_results(results_file)
  
  # Create all visualizations
  create_importance_ranking(results, output_dir)
  create_bubble_plot(results, output_dir)
  create_variable_type_comparison(results, output_dir)
  create_strategy_recommendations(results, output_dir)
  create_summary_dashboard(results, output_dir)
  
  cat("\nAll visualization plots created successfully!\n")
  cat("Plots saved to:", output_dir, "\n")
}

# Run main function
if (!interactive()) {
  main()
}