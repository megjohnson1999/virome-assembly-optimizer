#!/usr/bin/env Rscript

# PERMANOVA analysis to identify important variables affecting viral community structure
# Helps determine which metadata variables should guide co-assembly strategy

library(vegan)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# Function to load and validate data
load_data <- function(similarity_file, metadata_file) {
  cat("Loading similarity matrix from:", similarity_file, "\n")
  
  # Load similarity matrix
  if (!file.exists(similarity_file)) {
    stop("Similarity matrix file not found: ", similarity_file)
  }
  
  sim_matrix <- read.csv(similarity_file, row.names = 1, check.names = FALSE)
  
  # Convert to distance matrix
  dist_matrix <- as.dist(1 - sim_matrix)
  
  cat("Loaded similarity matrix with", nrow(sim_matrix), "samples\n")
  
  # Load metadata
  if (!file.exists(metadata_file)) {
    stop("Metadata file not found: ", metadata_file)
  }
  
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
  
  cat("Loaded metadata with", nrow(metadata), "samples and", ncol(metadata), "variables\n")
  
  # Ensure sample names match
  sample_names <- rownames(sim_matrix)
  metadata <- metadata[match(sample_names, metadata$Sample), ]
  
  if (any(is.na(metadata$Sample))) {
    missing_samples <- sample_names[is.na(metadata$Sample)]
    warning("Missing metadata for samples: ", paste(missing_samples, collapse = ", "))
    
    # Remove samples with missing metadata
    valid_samples <- !is.na(metadata$Sample)
    dist_matrix <- as.dist(as.matrix(dist_matrix)[valid_samples, valid_samples])
    metadata <- metadata[valid_samples, ]
    sample_names <- sample_names[valid_samples]
  }
  
  cat("Final dataset:", length(sample_names), "samples with complete data\n")
  
  return(list(distance_matrix = dist_matrix, metadata = metadata, sample_names = sample_names))
}

# Function to prepare variables for analysis
prepare_variables <- function(metadata) {
  cat("Preparing variables for analysis...\n")
  
  # Identify variable types
  factor_vars <- c()
  numeric_vars <- c()
  
  for (col in names(metadata)) {
    if (col == "Sample") next
    
    # Check if variable should be treated as factor
    if (is.character(metadata[[col]]) || is.factor(metadata[[col]]) || 
        length(unique(metadata[[col]])) <= 10) {
      factor_vars <- c(factor_vars, col)
      metadata[[col]] <- as.factor(metadata[[col]])
    } else {
      numeric_vars <- c(numeric_vars, col)
      metadata[[col]] <- as.numeric(metadata[[col]])
    }
  }
  
  cat("Factor variables:", paste(factor_vars, collapse = ", "), "\n")
  cat("Numeric variables:", paste(numeric_vars, collapse = ", "), "\n")
  
  return(list(
    metadata = metadata,
    factor_vars = factor_vars,
    numeric_vars = numeric_vars
  ))
}

# Function to run PERMANOVA for individual variables
run_individual_permanova <- function(dist_matrix, metadata, variables, n_perm = 999) {
  cat("Running individual PERMANOVA tests...\n")
  
  results <- data.frame(
    Variable = character(),
    R_squared = numeric(),
    F_statistic = numeric(),
    p_value = numeric(),
    Significance = character(),
    Variable_Type = character(),
    stringsAsFactors = FALSE
  )
  
  for (var in variables) {
    cat("  Testing variable:", var, "\n")
    
    # Skip if variable has too many missing values
    if (sum(is.na(metadata[[var]])) > nrow(metadata) * 0.5) {
      cat("    Skipping due to too many missing values\n")
      next
    }
    
    # Remove samples with missing values for this variable
    valid_samples <- !is.na(metadata[[var]])
    if (sum(valid_samples) < 10) {
      cat("    Skipping due to insufficient valid samples\n")
      next
    }
    
    temp_dist <- as.dist(as.matrix(dist_matrix)[valid_samples, valid_samples])
    temp_metadata <- metadata[valid_samples, ]
    
    # Run PERMANOVA
    tryCatch({
      if (is.factor(temp_metadata[[var]]) && length(levels(temp_metadata[[var]])) == 1) {
        cat("    Skipping: only one level in factor\n")
        next
      }
      
      adonis_result <- adonis2(temp_dist ~ temp_metadata[[var]], permutations = n_perm)
      
      # Extract results
      r_squared <- adonis_result$R2[1]
      f_stat <- adonis_result$F[1]
      p_val <- adonis_result$Pr[1]
      
      # Determine significance
      sig_level <- ""
      if (p_val <= 0.001) sig_level <- "***"
      else if (p_val <= 0.01) sig_level <- "**"
      else if (p_val <= 0.05) sig_level <- "*"
      else if (p_val <= 0.1) sig_level <- "."
      
      # Determine variable type
      var_type <- ifelse(is.factor(metadata[[var]]), "Factor", "Numeric")
      
      results <- rbind(results, data.frame(
        Variable = var,
        R_squared = r_squared,
        F_statistic = f_stat,
        p_value = p_val,
        Significance = sig_level,
        Variable_Type = var_type,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat("    Error testing variable", var, ":", e$message, "\n")
    })
  }
  
  # Sort by R-squared value
  results <- results[order(results$R_squared, decreasing = TRUE), ]
  
  return(results)
}

# Function to run combined PERMANOVA with multiple variables
run_combined_permanova <- function(dist_matrix, metadata, significant_vars, n_perm = 999) {
  cat("Running combined PERMANOVA with significant variables...\n")
  
  if (length(significant_vars) < 2) {
    cat("Not enough significant variables for combined analysis\n")
    return(NULL)
  }
  
  # Create formula
  formula_str <- paste("dist_matrix ~", paste(significant_vars, collapse = " + "))
  cat("Formula:", formula_str, "\n")
  
  # Run combined PERMANOVA
  tryCatch({
    combined_result <- adonis2(as.formula(formula_str), data = metadata, permutations = n_perm)
    return(combined_result)
  }, error = function(e) {
    cat("Error in combined PERMANOVA:", e$message, "\n")
    return(NULL)
  })
}

# Function to create visualization of results
create_permanova_plots <- function(results, output_dir) {
  cat("Creating PERMANOVA visualization plots...\n")
  
  # R-squared plot
  r2_plot <- ggplot(results, aes(x = reorder(Variable, R_squared), y = R_squared)) +
    geom_col(aes(fill = Variable_Type), alpha = 0.8) +
    geom_text(aes(label = Significance), hjust = -0.1, size = 4) +
    coord_flip() +
    labs(
      title = "Variable Importance in Explaining Viral Community Structure",
      subtitle = "Based on PERMANOVA R² values",
      x = "Variables",
      y = "R² (Proportion of variance explained)",
      fill = "Variable Type",
      caption = "Significance: *** p≤0.001, ** p≤0.01, * p≤0.05, . p≤0.1"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  r2_file <- file.path(output_dir, "permanova_r_squared.png")
  ggsave(r2_file, r2_plot, width = 10, height = 8, dpi = 300)
  cat("R-squared plot saved to:", r2_file, "\n")
  
  # p-value plot
  pval_plot <- ggplot(results, aes(x = reorder(Variable, -log10(p_value)), y = -log10(p_value))) +
    geom_col(aes(fill = Variable_Type), alpha = 0.8) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = 1) +
    geom_hline(yintercept = -log10(0.01), color = "darkred", linetype = "dashed", size = 1) +
    coord_flip() +
    labs(
      title = "Statistical Significance of Variables",
      subtitle = "Based on PERMANOVA p-values",
      x = "Variables",
      y = "-log10(p-value)",
      fill = "Variable Type",
      caption = "Dashed lines: red = p=0.05, dark red = p=0.01"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom"
    ) +
    scale_fill_brewer(type = "qual", palette = "Set2")
  
  pval_file <- file.path(output_dir, "permanova_p_values.png")
  ggsave(pval_file, pval_plot, width = 10, height = 8, dpi = 300)
  cat("P-value plot saved to:", pval_file, "\n")
  
  # Combined plot
  combined_plot <- ggplot(results, aes(x = R_squared, y = -log10(p_value))) +
    geom_point(aes(color = Variable_Type, size = F_statistic), alpha = 0.7) +
    geom_text(aes(label = Variable), hjust = 0, vjust = 0, nudge_x = 0.001, size = 3) +
    geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
    geom_vline(xintercept = 0.1, color = "blue", linetype = "dashed") +
    labs(
      title = "Variable Importance: Effect Size vs Statistical Significance",
      x = "R² (Effect Size)",
      y = "-log10(p-value)",
      color = "Variable Type",
      size = "F-statistic",
      caption = "Dashed lines: red = p=0.05, blue = R²=0.1 (recommended threshold)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "right"
    ) +
    scale_color_brewer(type = "qual", palette = "Set2")
  
  combined_file <- file.path(output_dir, "permanova_combined.png")
  ggsave(combined_file, combined_plot, width = 12, height = 8, dpi = 300)
  cat("Combined plot saved to:", combined_file, "\n")
}

# Main function
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("Usage: Rscript 01_permanova_analysis.R <similarity_matrix.csv> <metadata.csv> [output_dir]\n")
    cat("Using default paths...\n")
    similarity_file <- "results/similarity_analysis/similarity_matrix.csv"
    metadata_file <- "examples/metadata_template.csv"
  } else {
    similarity_file <- args[1]
    metadata_file <- args[2]
  }
  
  if (length(args) < 3) {
    output_dir <- "results/variable_analysis"
  } else {
    output_dir <- args[3]
  }
  
  cat("Similarity matrix:", similarity_file, "\n")
  cat("Metadata file:", metadata_file, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  data <- load_data(similarity_file, metadata_file)
  
  # Prepare variables
  var_info <- prepare_variables(data$metadata)
  all_vars <- c(var_info$factor_vars, var_info$numeric_vars)
  
  # Run individual PERMANOVA tests
  results <- run_individual_permanova(
    data$distance_matrix, 
    var_info$metadata, 
    all_vars
  )
  
  # Save results
  results_file <- file.path(output_dir, "permanova_results.csv")
  write.csv(results, results_file, row.names = FALSE)
  cat("PERMANOVA results saved to:", results_file, "\n")
  
  # Print summary
  cat("\nPERMANOVA Results Summary:\n")
  cat("=========================\n")
  print(results)
  
  # Identify significant variables
  significant_vars <- results$Variable[results$p_value <= 0.05 & results$R_squared >= 0.1]
  
  cat("\nSignificant variables (p ≤ 0.05, R² ≥ 0.1):\n")
  if (length(significant_vars) > 0) {
    cat(paste(significant_vars, collapse = ", "), "\n")
    
    # Save significant variables
    sig_vars_file <- file.path(output_dir, "significant_variables.txt")
    writeLines(significant_vars, sig_vars_file)
    cat("Significant variables saved to:", sig_vars_file, "\n")
    
    # Run combined PERMANOVA
    combined_result <- run_combined_permanova(
      data$distance_matrix,
      var_info$metadata,
      significant_vars
    )
    
    if (!is.null(combined_result)) {
      combined_file <- file.path(output_dir, "combined_permanova.txt")
      sink(combined_file)
      cat("Combined PERMANOVA Results\n")
      cat("==========================\n\n")
      print(combined_result)
      sink()
      cat("Combined PERMANOVA results saved to:", combined_file, "\n")
    }
  } else {
    cat("No variables meet the significance threshold\n")
  }
  
  # Create visualizations
  create_permanova_plots(results, output_dir)
  
  # Generate recommendations
  recommendations_file <- file.path(output_dir, "assembly_recommendations.txt")
  sink(recommendations_file)
  
  cat("Assembly Strategy Recommendations Based on PERMANOVA Analysis\n")
  cat("=============================================================\n\n")
  
  cat("Important Variables (R² ≥ 0.1, p ≤ 0.05):\n")
  if (length(significant_vars) > 0) {
    for (var in significant_vars) {
      var_result <- results[results$Variable == var, ]
      cat(sprintf("- %s: R² = %.3f, p = %.3f %s\n", 
                  var, var_result$R_squared, var_result$p_value, var_result$Significance))
    }
    
    cat("\nRecommendations:\n")
    cat("1. Use these variables to guide co-assembly group formation\n")
    cat("2. Samples with similar values for these variables are good candidates for co-assembly\n")
    cat("3. Consider stratifying analyses by these variables\n")
    cat("4. Prioritize biological variables over technical ones when possible\n")
    
  } else {
    cat("No variables show strong association with community structure.\n")
    cat("\nRecommendations:\n")
    cat("1. Consider individual sample assemblies\n")
    cat("2. Use k-mer similarity as the primary grouping criterion\n")
    cat("3. Investigate potential batch effects or data quality issues\n")
  }
  
  sink()
  
  cat("Assembly recommendations saved to:", recommendations_file, "\n")
  cat("\nPERMANOVA analysis completed successfully!\n")
}

# Run main function
if (!interactive()) {
  main()
}