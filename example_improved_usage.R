# ==============================================================================
# EXAMPLE: How to use the improved DPY mixture model
# ==============================================================================

source("dpy_mixture_improved.R")

set.seed(42)

# ==============================================================================
# Generate example data (same structure as your plot)
# ==============================================================================

cat("Generating example data...\n")

# Group 1: Multimodal distribution (uniform-like)
Y1 <- c(
  rnorm(30, -10, 1.5),
  rnorm(30, -5, 1.5),
  rnorm(30, 0, 1.5),
  rnorm(30, 5, 1.5),
  rnorm(30, 10, 1.5)
)

# Group 2: Bimodal distribution (peak at 0 and peak at 10)
Y2 <- c(
  rnorm(100, 0, 1.5),   # Main peak at 0
  rnorm(60, 10, 1.5),   # Secondary peak at 10
  rnorm(20, -10, 2)     # Small peak at -10
)

data <- list(Y1 = Y1, Y2 = Y2)

cat(sprintf("Group 1: N=%d, mean=%.2f, sd=%.2f\n", length(Y1), mean(Y1), sd(Y1)))
cat(sprintf("Group 2: N=%d, mean=%.2f, sd=%.2f\n", length(Y2), mean(Y2), sd(Y2)))

# ==============================================================================
# OPTION 1: Run with SIMPLIFIED model (RECOMMENDED)
# ==============================================================================

cat("\n========================================\n")
cat("OPTION 1: SIMPLIFIED MODEL (Recommended)\n")
cat("========================================\n\n")

# This uses a simpler parametrization:
# Y_{it} ~ N(μ_k, σ²_k)  instead of  N(μ_0k + μ_{ik}, σ²_0k * σ²_{ik})

cat("Running MCMC sampler (simplified model)...\n")

res_simplified <- run_DPY_sampler_improved(
  data = data,
  hyperparams = initialize_hyperparameters_improved(),
  use_simplified_model = TRUE,
  n_iter = 5000,    # Increase to 10000 for better results
  burnin = 2500,
  thin = 5,
  stampa = TRUE
)

cat("\nDiagnosing results...\n")
diagnose_results(res_simplified)

cat("\nComputing best partition (Dahl's method)...\n")
best_simplified <- dahl(res_simplified)

cat("\nPlotting density estimates...\n")
plot_density_improved(best_simplified, res_simplified, threshold = 0.005,
                     title = "Simplified Model - Best Partition")

# ==============================================================================
# OPTION 2: Run with ORIGINAL complex model + improved hyperparameters
# ==============================================================================

cat("\n\n========================================\n")
cat("OPTION 2: COMPLEX MODEL (Original)\n")
cat("========================================\n\n")

cat("Running MCMC sampler (complex model with improved hyperparameters)...\n")

res_complex <- run_DPY_sampler_improved(
  data = data,
  hyperparams = initialize_hyperparameters_improved(),
  use_simplified_model = FALSE,
  n_iter = 5000,
  burnin = 2500,
  thin = 5,
  stampa = TRUE
)

cat("\nDiagnosing results...\n")
diagnose_results(res_complex)

cat("\nComputing best partition (Dahl's method)...\n")
best_complex <- dahl(res_complex)

cat("\nPlotting density estimates...\n")
plot_density_improved(best_complex, res_complex, threshold = 0.005,
                     title = "Complex Model - Best Partition")

# ==============================================================================
# OPTION 3: Quick test with default original sampler
# ==============================================================================

cat("\n\n========================================\n")
cat("OPTION 3: ORIGINAL SAMPLER (for comparison)\n")
cat("========================================\n\n")

cat("Running original sampler (for comparison)...\n")

res_original <- run_DPY_sampler(
  data = data,
  hyperparams = initialize_hyperparameters(),  # Old hyperparameters
  n_iter = 5000,
  burnin = 2500,
  thin = 5,
  stampa = TRUE
)

cat("\nDiagnosing results...\n")
diagnose_results(res_original)

cat("\nComputing best partition (Dahl's method)...\n")
best_original <- dahl(res_original)

cat("\nPlotting density estimates...\n")
plot_density_improved(best_original, res_original, threshold = 0.005,
                     title = "Original Model - Best Partition")

# ==============================================================================
# COMPARISON SUMMARY
# ==============================================================================

cat("\n\n============================================\n")
cat("COMPARISON SUMMARY\n")
cat("============================================\n\n")

cat("Active Clusters:\n")
cat(sprintf("Simplified - Group 1: %.1f, Group 2: %.1f\n",
            mean(res_simplified$n_active[,1]), mean(res_simplified$n_active[,2])))
cat(sprintf("Complex    - Group 1: %.1f, Group 2: %.1f\n",
            mean(res_complex$n_active[,1]), mean(res_complex$n_active[,2])))
cat(sprintf("Original   - Group 1: %.1f, Group 2: %.1f\n",
            mean(res_original$n_active[,1]), mean(res_original$n_active[,2])))

cat("\nAcceptance Rates:\n")
cat(sprintf("Simplified: %.3f\n", mean(res_simplified$acc_rate, na.rm=TRUE)))
cat(sprintf("Complex:    %.3f\n", mean(res_complex$acc_rate, na.rm=TRUE)))
cat(sprintf("Original:   %.3f\n", mean(res_original$acc_rate, na.rm=TRUE)))

cat("\n\n============================================\n")
cat("RECOMMENDATIONS:\n")
cat("============================================\n\n")
cat("1. Try SIMPLIFIED MODEL first - it's more interpretable\n")
cat("2. If fit is still poor, increase n_iter to 10000-20000\n")
cat("3. Lower threshold in plot_density (try 0.001 or 0.002)\n")
cat("4. Check diagnostic plots for convergence\n")
cat("5. For multimodal data, may need more clusters (increase K_start)\n\n")

# ==============================================================================
# ADDITIONAL DIAGNOSTIC PLOTS
# ==============================================================================

if(requireNamespace("ggplot2", quietly = TRUE)) {
  cat("Creating additional diagnostic plots...\n\n")

  # Plot 1: Trace of concentration parameters
  library(ggplot2)

  n_iter <- nrow(res_simplified$psi)
  df_psi <- data.frame(
    Iteration = rep(1:n_iter, 3),
    Value = c(res_simplified$psi[,1], res_simplified$psi[,2], res_simplified$psi[,3]),
    Parameter = rep(c("alpha1", "alpha2", "l"), each = n_iter)
  )

  p_psi <- ggplot(df_psi, aes(x = Iteration, y = Value, color = Parameter)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(title = "Trace of Concentration Parameters") +
    theme(legend.position = "none")

  print(p_psi)

  # Plot 2: Number of active clusters over iterations
  df_K <- data.frame(
    Iteration = rep(1:n_iter, 2),
    K = c(res_simplified$n_active[,1], res_simplified$n_active[,2]),
    Group = rep(c("Group 1", "Group 2"), each = n_iter)
  )

  p_K <- ggplot(df_K, aes(x = Iteration, y = K, color = Group)) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Number of Active Components Over Iterations") +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e"))

  print(p_K)

  # Plot 3: Variance trace (if simplified model)
  if(!is.null(res_simplified$trace_variances)) {
    df_var <- data.frame(
      Iteration = rep(1:nrow(res_simplified$trace_variances), 5),
      Variance = c(res_simplified$trace_variances[,1],
                   res_simplified$trace_variances[,2],
                   res_simplified$trace_variances[,3],
                   res_simplified$trace_variances[,4],
                   res_simplified$trace_variances[,5]),
      Component = rep(paste0("K=", 1:5), each = nrow(res_simplified$trace_variances))
    )

    p_var <- ggplot(df_var, aes(x = Iteration, y = Variance, color = Component)) +
      geom_line(alpha = 0.7) +
      theme_minimal() +
      labs(title = "Variance Trace (First 5 Components)", y = "sigma^2")

    print(p_var)
  }
}

cat("\n\nExample completed!\n")
cat("========================================\n\n")

# Save results (optional)
# saveRDS(res_simplified, "results_simplified.rds")
# saveRDS(res_complex, "results_complex.rds")
