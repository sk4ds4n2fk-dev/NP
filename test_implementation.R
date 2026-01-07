#!/usr/bin/env Rscript
# Test script for Beta-Product Dependent Pitman-Yor implementation
# This script runs basic tests to verify the implementation is working correctly

cat("Testing Beta-Product Dependent Pitman-Yor Implementation\n")
cat("========================================================\n\n")

# Source the implementation
source("beta_product_dpy.R")
source("mcmc_sampler.R")

# Set seed for reproducibility
set.seed(42)

################################################################################
# Test 1: Stick-Breaking Process
################################################################################

cat("Test 1: Stick-Breaking Process\n")
cat("-------------------------------\n")

alpha <- 0.5
theta <- 1.0
K <- 10

weights <- stick_breaking_py(alpha, theta, K)

cat(sprintf("Generated %d stick-breaking weights\n", length(weights)))
cat(sprintf("Sum of weights: %.4f (should be close to 1)\n", sum(weights)))
cat(sprintf("First weight: %.4f\n", weights[1]))
cat(sprintf("Last weight: %.4f\n", weights[K]))

if (abs(sum(weights) - 1) < 0.1 && all(weights >= 0) && all(weights <= 1)) {
  cat("✓ Test 1 PASSED\n\n")
} else {
  cat("✗ Test 1 FAILED\n\n")
  stop("Stick-breaking test failed")
}


################################################################################
# Test 2: Beta-Product Construction
################################################################################

cat("Test 2: Beta-Product Construction\n")
cat("----------------------------------\n")

d <- 3
K <- 10
alpha <- 0.5
theta <- 1

# Generate V with proper stick-breaking structure for each process
V <- matrix(0, nrow = K, ncol = d)
for (j in 1:d) {
  for (k in 1:K) {
    V[k, j] <- rbeta(1, 1 - alpha, theta + k * alpha)
  }
}
W <- stick_breaking_py(0.5, 1, K)

dependent_weights <- beta_product_construction(V, W)

cat(sprintf("Generated dependent weights for %d processes\n", d))
cat(sprintf("Dimensions: %d x %d\n", nrow(dependent_weights), ncol(dependent_weights)))

# Check that weights are positive and matrix has correct dimensions
sums <- colSums(dependent_weights)
cat(sprintf("Column sums: %s\n", paste(sprintf("%.3f", sums), collapse=", ")))

if (all(dim(dependent_weights) == c(K, d)) && all(dependent_weights >= 0) && all(sums > 0)) {
  cat("✓ Test 2 PASSED\n\n")
} else {
  cat("✗ Test 2 FAILED\n\n")
  stop("Beta-product construction test failed")
}


################################################################################
# Test 3: DPY Generation
################################################################################

cat("Test 3: DPY Generation\n")
cat("----------------------\n")

d <- 2
K <- 15
alpha_vec <- c(0.5, 0.5)
theta_vec <- c(1, 1)
alpha0 <- 0.5
theta0 <- 1

dpy <- generate_dpy(d, K, alpha_vec, theta_vec, alpha0, theta0)

cat(sprintf("Generated DPY with %d processes and %d components\n", d, K))
cat(sprintf("Weight matrix dimensions: %d x %d\n",
            nrow(dpy$weights), ncol(dpy$weights)))
cat(sprintf("Atom matrix dimensions: %d x %d\n",
            nrow(dpy$atoms), ncol(dpy$atoms)))

# Check structure
if (all(dim(dpy$weights) == c(K, d)) &&
    all(dim(dpy$atoms) == c(K, d)) &&
    length(dpy$base_weights) == K) {
  cat("✓ Test 3 PASSED\n\n")
} else {
  cat("✗ Test 3 FAILED\n\n")
  stop("DPY generation test failed")
}


################################################################################
# Test 4: Sample Allocations
################################################################################

cat("Test 4: Sample Allocations\n")
cat("--------------------------\n")

n <- 50
y <- rnorm(n)
atoms <- rnorm(K)
weights <- stick_breaking_py(0.5, 1, K)
sigma2 <- 1

allocations <- sample_allocations(y, atoms, weights, sigma2)

cat(sprintf("Allocated %d observations to %d clusters\n", n, length(unique(allocations))))
cat(sprintf("Allocation range: [%d, %d]\n", min(allocations), max(allocations)))
cat(sprintf("Number of unique clusters used: %d\n", length(unique(allocations))))

if (length(allocations) == n &&
    all(allocations >= 1) &&
    all(allocations <= K)) {
  cat("✓ Test 4 PASSED\n\n")
} else {
  cat("✗ Test 4 FAILED\n\n")
  stop("Sample allocations test failed")
}


################################################################################
# Test 5: Update Atoms
################################################################################

cat("Test 5: Update Atoms\n")
cat("--------------------\n")

updated_atoms <- update_atoms(y, allocations, K, sigma2)

cat(sprintf("Updated %d atom locations\n", length(updated_atoms)))
cat(sprintf("Atom range: [%.2f, %.2f]\n", min(updated_atoms), max(updated_atoms)))

if (length(updated_atoms) == K && all(is.finite(updated_atoms))) {
  cat("✓ Test 5 PASSED\n\n")
} else {
  cat("✗ Test 5 FAILED\n\n")
  stop("Update atoms test failed")
}


################################################################################
# Test 6: Small MCMC Run
################################################################################

cat("Test 6: Small MCMC Run\n")
cat("----------------------\n")

# Generate simple test data
d <- 2
n <- 100

# Two clusters per process
y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
y <- list(y1, y2)

cat("Running MCMC with 200 iterations...\n")

results <- gibbs_sampler_dpy(
  y = y,
  K = 10,
  n_iter = 200,
  burn_in = 100,
  thin = 1,
  alpha_vec = c(0.5, 0.5),
  theta_vec = c(1, 1),
  alpha0 = 0.5,
  theta0 = 1,
  sigma2 = 1,
  estimate_sigma = TRUE,
  verbose = FALSE
)

cat(sprintf("Number of posterior samples: %d\n", length(results$sigma2)))
cat(sprintf("Posterior mean sigma2: %.3f\n", mean(results$sigma2)))
cat(sprintf("Mean clusters (Process 1): %.1f\n",
            mean(sapply(results$n_clusters, function(x) x[1]))))
cat(sprintf("Mean clusters (Process 2): %.1f\n",
            mean(sapply(results$n_clusters, function(x) x[2]))))

if (length(results$sigma2) == 100 &&
    all(is.finite(results$sigma2)) &&
    all(results$sigma2 > 0)) {
  cat("✓ Test 6 PASSED\n\n")
} else {
  cat("✗ Test 6 FAILED\n\n")
  stop("MCMC run test failed")
}


################################################################################
# Test 7: Model Diagnostics
################################################################################

cat("Test 7: Model Diagnostics\n")
cat("-------------------------\n")

# Compute WAIC
waic <- compute_waic(y, results)
cat(sprintf("WAIC: %.2f\n", waic$waic))
cat(sprintf("lppd: %.2f\n", waic$lppd))
cat(sprintf("p_waic: %.2f\n", waic$p_waic))

# Compute DIC
dic <- compute_dic(results)
cat(sprintf("DIC: %.2f\n", dic$DIC))

if (is.finite(waic$waic) && is.finite(dic$DIC)) {
  cat("✓ Test 7 PASSED\n\n")
} else {
  cat("✗ Test 7 FAILED\n\n")
  stop("Model diagnostics test failed")
}


################################################################################
# Test 8: Posterior Predictive
################################################################################

cat("Test 8: Posterior Predictive\n")
cat("-----------------------------\n")

pred_samples <- posterior_predictive(results, n_pred = 50, process_idx = 1)

cat(sprintf("Generated predictive samples: %d x %d\n",
            nrow(pred_samples), ncol(pred_samples)))
cat(sprintf("Mean prediction: %.2f\n", mean(pred_samples)))
cat(sprintf("SD of predictions: %.2f\n", sd(as.vector(pred_samples))))

if (all(dim(pred_samples) == c(100, 50)) && all(is.finite(pred_samples))) {
  cat("✓ Test 8 PASSED\n\n")
} else {
  cat("✗ Test 8 FAILED\n\n")
  stop("Posterior predictive test failed")
}


################################################################################
# Test 9: Weight Correlation
################################################################################

cat("Test 9: Weight Correlation\n")
cat("--------------------------\n")

mean_weights <- apply(results$weights, c(2, 3), mean)
cor_matrix <- weight_correlation(mean_weights)

cat("Weight correlation matrix:\n")
print(round(cor_matrix, 3))

if (all(dim(cor_matrix) == c(d, d)) && all(diag(cor_matrix) == 1)) {
  cat("✓ Test 9 PASSED\n\n")
} else {
  cat("✗ Test 9 FAILED\n\n")
  stop("Weight correlation test failed")
}


################################################################################
# Test 10: ESS Computation
################################################################################

cat("Test 10: ESS Computation\n")
cat("------------------------\n")

# Compute ESS for posterior mean weights (effective number of mixture components)
ess <- compute_ess(mean_weights)

cat(sprintf("Effective number of components for each process: %s\n",
            paste(sprintf("%.2f", ess), collapse=", ")))

if (length(ess) == d && all(ess > 0) && all(ess <= K)) {
  cat("✓ Test 10 PASSED\n\n")
} else {
  cat("✗ Test 10 FAILED\n\n")
  stop("ESS computation test failed")
}


################################################################################
# Summary
################################################################################

cat("\n")
cat("========================================================\n")
cat("ALL TESTS PASSED!\n")
cat("========================================================\n")
cat("\nThe Beta-Product Dependent Pitman-Yor implementation\n")
cat("is working correctly.\n\n")
cat("You can now run the examples:\n")
cat("  source('examples.R')\n")
cat("  quick_test()\n")
cat("  run_all_examples()\n\n")
