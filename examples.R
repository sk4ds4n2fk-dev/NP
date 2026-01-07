# Examples and Applications for Beta-Product Dependent Pitman-Yor Processes
# Based on Bassetti, Casarin, and Leisen (2014)

source("mcmc_sampler.R")

# Set seed for reproducibility
set.seed(123)

################################################################################
# Example 1: Simulated Data from DPY Mixture
################################################################################

example1_simulated_data <- function() {
  cat("=== Example 1: Simulated Data from DPY Mixture ===\n\n")

  # Simulation parameters
  d <- 3  # Number of dependent processes
  K_true <- 5  # True number of clusters
  n <- c(200, 250, 180)  # Sample sizes for each process

  # True parameters
  alpha_true <- c(0.5, 0.5, 0.5)
  theta_true <- c(1, 1, 1)
  alpha0_true <- 0.5
  theta0_true <- 1
  sigma2_true <- 0.5

  # Generate true DPY process
  cat("Generating true DPY process...\n")
  dpy_true <- generate_dpy(d, K_true, alpha_true, theta_true, alpha0_true, theta0_true)

  # Generate data
  cat("Generating data from mixture model...\n")
  y <- vector("list", d)
  true_allocations <- vector("list", d)

  for (j in 1:d) {
    true_allocations[[j]] <- sample(1:K_true, n[j], replace = TRUE,
                                     prob = dpy_true$weights[, j])
    y[[j]] <- numeric(n[j])
    for (i in 1:n[j]) {
      k <- true_allocations[[j]][i]
      y[[j]][i] <- rnorm(1, mean = dpy_true$atoms[k, j], sd = sqrt(sigma2_true))
    }
  }

  # Print data summary
  cat("\nData summary:\n")
  for (j in 1:d) {
    cat(sprintf("Process %d: n=%d, mean=%.2f, sd=%.2f\n",
                j, n[j], mean(y[[j]]), sd(y[[j]])))
  }

  # Run MCMC sampler
  cat("\nRunning MCMC sampler...\n")
  results <- gibbs_sampler_dpy(
    y = y,
    K = 20,
    n_iter = 2000,
    burn_in = 1000,
    thin = 2,
    alpha_vec = c(0.5, 0.5, 0.5),
    theta_vec = c(1, 1, 1),
    alpha0 = 0.5,
    theta0 = 1,
    sigma2 = 1,
    estimate_sigma = TRUE,
    verbose = TRUE
  )

  # Posterior summaries
  cat("\n=== Posterior Summaries ===\n")
  cat(sprintf("Posterior mean of sigma2: %.3f (True: %.3f)\n",
              mean(results$sigma2), sigma2_true))
  cat(sprintf("95%% CI for sigma2: [%.3f, %.3f]\n",
              quantile(results$sigma2, 0.025), quantile(results$sigma2, 0.975)))

  # Number of clusters
  cat("\nPosterior number of clusters:\n")
  for (j in 1:d) {
    n_clust <- sapply(results$n_clusters, function(x) x[j])
    cat(sprintf("Process %d: mean=%.1f, median=%d (True: %d)\n",
                j, mean(n_clust), median(n_clust),
                length(unique(true_allocations[[j]]))))
  }

  # Model comparison
  cat("\n=== Model Fit ===\n")
  waic <- compute_waic(y, results)
  cat(sprintf("WAIC: %.2f\n", waic$waic))
  cat(sprintf("Effective number of parameters: %.2f\n", waic$p_waic))

  # Plotting (if available)
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_results_example1(y, results, dpy_true)
  }

  return(list(y = y, results = results, truth = dpy_true))
}


################################################################################
# Example 2: Business Cycle Application
################################################################################

example2_business_cycles <- function() {
  cat("\n\n=== Example 2: Business Cycle Application ===\n\n")

  # Simulate business cycle data (GDP growth rates)
  # This mimics the US and EU business cycle analysis in the paper

  n_years <- 100
  quarters_per_year <- 4
  n <- n_years * quarters_per_year

  # Simulate two correlated business cycles (US and EU)
  cat("Simulating business cycle data (US and EU GDP growth)...\n")

  # Generate regime switches (expansion vs. recession)
  regime_switches <- c(1, 80, 150, 280, 350)  # Recession start points
  regimes <- rep(1, n)  # 1 = expansion, 2 = recession
  for (i in seq(1, length(regime_switches), by = 2)) {
    if (i < length(regime_switches)) {
      regimes[regime_switches[i]:regime_switches[i+1]] <- 2
    } else {
      regimes[regime_switches[i]:n] <- 2
    }
  }

  # Generate growth rates with different means in different regimes
  us_growth <- numeric(n)
  eu_growth_independent <- numeric(n)

  for (i in 1:n) {
    if (regimes[i] == 1) {  # Expansion
      us_growth[i] <- rnorm(1, mean = 3.0, sd = 1.5)
      eu_growth_independent[i] <- rnorm(1, mean = 2.5, sd = 1.2)
    } else {  # Recession
      us_growth[i] <- rnorm(1, mean = -1.5, sd = 2.0)
      eu_growth_independent[i] <- rnorm(1, mean = -1.0, sd = 1.5)
    }
  }

  # Add correlation between US and EU
  correlation <- 0.7
  eu_growth <- correlation * us_growth + sqrt(1 - correlation^2) * eu_growth_independent

  y <- list(us_growth, eu_growth)

  # Print data summary
  cat("\nData summary:\n")
  cat(sprintf("US: mean=%.2f, sd=%.2f\n", mean(us_growth), sd(us_growth)))
  cat(sprintf("EU: mean=%.2f, sd=%.2f\n", mean(eu_growth), sd(eu_growth)))
  cat(sprintf("Correlation: %.2f\n", cor(us_growth, eu_growth)))

  # Run MCMC sampler
  cat("\nRunning MCMC sampler for business cycle model...\n")
  results <- gibbs_sampler_dpy(
    y = y,
    K = 15,
    n_iter = 3000,
    burn_in = 1500,
    thin = 3,
    alpha_vec = c(0.3, 0.3),
    theta_vec = c(2, 2),
    alpha0 = 0.3,
    theta0 = 2,
    sigma2 = 2,
    estimate_sigma = TRUE,
    verbose = TRUE
  )

  # Posterior summaries
  cat("\n=== Posterior Summaries ===\n")
  cat(sprintf("Posterior mean of sigma2: %.3f\n", mean(results$sigma2)))

  # Identify business cycle regimes
  cat("\nIdentified business cycle regimes:\n")
  for (j in 1:2) {
    region <- ifelse(j == 1, "US", "EU")
    n_clust <- sapply(results$n_clusters, function(x) x[j])
    cat(sprintf("%s: mean clusters=%.1f, median=%d\n",
                region, mean(n_clust), median(n_clust)))
  }

  # Compute cluster means (regime characteristics)
  cat("\nCluster characteristics (regime means):\n")
  final_atoms <- results$atoms[dim(results$atoms)[1], , ]
  for (j in 1:2) {
    region <- ifelse(j == 1, "US", "EU")
    cat(sprintf("%s regimes: %s\n", region,
                paste(sprintf("%.2f", final_atoms[, j]), collapse=", ")))
  }

  return(list(y = y, results = results, regimes = regimes))
}


################################################################################
# Example 3: Vector Autoregressive (VAR) Model with DPY Prior
################################################################################

example3_var_model <- function() {
  cat("\n\n=== Example 3: VAR Model with DPY Prior ===\n\n")

  # Simulate VAR(1) data with regime changes
  d <- 3  # Number of time series
  n <- 200  # Length of time series

  cat("Simulating VAR(1) data with regime changes...\n")

  # Initialize
  Y <- matrix(0, nrow = n, ncol = d)
  Y[1, ] <- rnorm(d)

  # Regime-dependent VAR coefficients
  # Regime 1: Positive persistence
  A1 <- matrix(c(0.7, 0.1, 0.0,
                 0.1, 0.6, 0.1,
                 0.0, 0.1, 0.5), nrow = d, byrow = TRUE)

  # Regime 2: Negative persistence
  A2 <- matrix(c(-0.5, 0.2, 0.1,
                 0.2, -0.4, 0.2,
                 0.1, 0.2, -0.3), nrow = d, byrow = TRUE)

  # Generate data with regime switches
  regimes <- rep(1, n)
  regimes[100:150] <- 2  # Regime switch

  for (t in 2:n) {
    if (regimes[t] == 1) {
      Y[t, ] <- A1 %*% Y[t-1, ] + rnorm(d, sd = 0.5)
    } else {
      Y[t, ] <- A2 %*% Y[t-1, ] + rnorm(d, sd = 0.5)
    }
  }

  # Prepare data for DPY mixture (using differences or residuals)
  y <- lapply(1:d, function(j) Y[, j])

  cat("\nData summary:\n")
  for (j in 1:d) {
    cat(sprintf("Series %d: mean=%.2f, sd=%.2f\n", j, mean(y[[j]]), sd(y[[j]])))
  }

  # Run MCMC sampler
  cat("\nRunning MCMC sampler for VAR model...\n")
  results <- gibbs_sampler_dpy(
    y = y,
    K = 12,
    n_iter = 2000,
    burn_in = 1000,
    thin = 2,
    alpha_vec = c(0.4, 0.4, 0.4),
    theta_vec = c(1.5, 1.5, 1.5),
    alpha0 = 0.4,
    theta0 = 1.5,
    sigma2 = 1,
    estimate_sigma = TRUE,
    verbose = TRUE
  )

  # Posterior summaries
  cat("\n=== Posterior Summaries ===\n")
  for (j in 1:d) {
    n_clust <- sapply(results$n_clusters, function(x) x[j])
    cat(sprintf("Series %d: mean clusters=%.1f\n", j, mean(n_clust)))
  }

  # Compute weight correlations
  cat("\nWeight correlations across series:\n")
  mean_weights <- apply(results$weights, c(2, 3), mean)
  cor_matrix <- cor(mean_weights)
  print(round(cor_matrix, 3))

  return(list(Y = Y, y = y, results = results, regimes = regimes))
}


################################################################################
# Plotting Functions
################################################################################

plot_results_example1 <- function(y, results, truth) {
  library(ggplot2)

  # Plot trace of sigma2
  df_sigma <- data.frame(
    iteration = 1:length(results$sigma2),
    sigma2 = results$sigma2
  )

  p1 <- ggplot(df_sigma, aes(x = iteration, y = sigma2)) +
    geom_line() +
    labs(title = "Trace Plot of sigma2", x = "Iteration", y = "sigma2") +
    theme_minimal()

  print(p1)

  # Plot histogram of sigma2
  p2 <- ggplot(df_sigma, aes(x = sigma2)) +
    geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
    labs(title = "Posterior Distribution of sigma2", x = "sigma2", y = "Frequency") +
    theme_minimal()

  print(p2)

  # Plot data with estimated clusters for first process
  if (length(y) > 0) {
    final_alloc <- results$allocations[[length(results$allocations)]][[1]]
    df_data <- data.frame(
      index = 1:length(y[[1]]),
      value = y[[1]],
      cluster = as.factor(final_alloc)
    )

    p3 <- ggplot(df_data, aes(x = index, y = value, color = cluster)) +
      geom_point() +
      labs(title = "Data with Estimated Clusters (Process 1)",
           x = "Observation", y = "Value") +
      theme_minimal()

    print(p3)
  }
}


################################################################################
# Run all examples
################################################################################

run_all_examples <- function() {
  cat("Running all examples for Beta-Product Dependent Pitman-Yor Processes\n")
  cat("====================================================================\n\n")

  # Example 1
  ex1 <- example1_simulated_data()

  # Example 2
  ex2 <- example2_business_cycles()

  # Example 3
  ex3 <- example3_var_model()

  cat("\n\nAll examples completed!\n")

  return(list(
    example1 = ex1,
    example2 = ex2,
    example3 = ex3
  ))
}


################################################################################
# Quick test function
################################################################################

quick_test <- function() {
  cat("Quick test of Beta-Product DPY implementation\n")
  cat("=============================================\n\n")

  # Generate simple test data
  set.seed(42)
  d <- 2
  n <- 100

  # Two processes with some dependence
  y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
  y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
  y <- list(y1, y2)

  cat("Running quick MCMC test...\n")
  results <- gibbs_sampler_dpy(
    y = y,
    K = 10,
    n_iter = 500,
    burn_in = 250,
    thin = 1,
    verbose = FALSE
  )

  cat("\nResults:\n")
  cat(sprintf("Mean sigma2: %.3f\n", mean(results$sigma2)))
  cat(sprintf("Mean clusters (Process 1): %.1f\n",
              mean(sapply(results$n_clusters, function(x) x[1]))))
  cat(sprintf("Mean clusters (Process 2): %.1f\n",
              mean(sapply(results$n_clusters, function(x) x[2]))))

  cat("\nQuick test completed successfully!\n")

  return(results)
}


# Uncomment to run examples:
# results <- run_all_examples()
# Or run quick test:
# test_results <- quick_test()
