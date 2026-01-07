# Script to plot data colored by cluster assignments
# Run this locally in R after running the MCMC sampler

# Function to plot cluster assignments
plot_cluster_results <- function(y, results, process_idx = 1, iteration = "last") {

  # Get allocation for specified iteration
  if (iteration == "last") {
    alloc_idx <- length(results$allocations)
  } else {
    alloc_idx <- iteration
  }

  allocations <- results$allocations[[alloc_idx]][[process_idx]]
  atoms <- results$atoms[alloc_idx, , process_idx]

  # Get unique clusters
  unique_clusters <- sort(unique(allocations))
  n_clusters <- length(unique_clusters)

  # Create color palette
  colors <- rainbow(n_clusters)
  cluster_colors <- colors[match(allocations, unique_clusters)]

  # Plot
  plot(1:length(y[[process_idx]]), y[[process_idx]],
       col = cluster_colors,
       pch = 19,
       cex = 1.2,
       xlab = "Observation Index",
       ylab = "Value",
       main = sprintf("Process %d - Cluster Assignments (Iteration %d)",
                      process_idx, alloc_idx))

  # Add horizontal lines for cluster means
  for (i in 1:n_clusters) {
    cluster <- unique_clusters[i]
    cluster_mean <- atoms[cluster]
    abline(h = cluster_mean, col = colors[i], lwd = 2, lty = 2)
  }

  # Add legend
  legend("topright",
         legend = paste("Cluster", unique_clusters,
                       sprintf("(n=%d, μ=%.2f)",
                               sapply(unique_clusters, function(k) sum(allocations == k)),
                               atoms[unique_clusters])),
         col = colors,
         pch = 19,
         cex = 0.8)

  # Print summary
  cat(sprintf("\nProcess %d Summary:\n", process_idx))
  cat(sprintf("Number of clusters used: %d\n", n_clusters))
  for (i in 1:n_clusters) {
    cluster <- unique_clusters[i]
    cat(sprintf("  Cluster %d: n=%d, mean=%.2f\n",
                cluster,
                sum(allocations == cluster),
                atoms[cluster]))
  }
}


# Function to plot both processes side-by-side
plot_both_processes <- function(y, results, iteration = "last") {
  par(mfrow = c(1, 2))
  plot_cluster_results(y, results, process_idx = 1, iteration = iteration)
  plot_cluster_results(y, results, process_idx = 2, iteration = iteration)
  par(mfrow = c(1, 1))
}


# Function to plot data distribution with cluster overlay
plot_histogram_clusters <- function(y, results, process_idx = 1, iteration = "last") {

  # Get allocation for specified iteration
  if (iteration == "last") {
    alloc_idx <- length(results$allocations)
  } else {
    alloc_idx <- iteration
  }

  allocations <- results$allocations[[alloc_idx]][[process_idx]]
  atoms <- results$atoms[alloc_idx, , process_idx]

  # Get unique clusters
  unique_clusters <- sort(unique(allocations))
  n_clusters <- length(unique_clusters)

  # Create histogram
  hist(y[[process_idx]],
       breaks = 30,
       col = "lightgray",
       border = "white",
       main = sprintf("Process %d - Data Distribution with Clusters", process_idx),
       xlab = "Value",
       ylab = "Frequency",
       prob = FALSE)

  # Add vertical lines for cluster means
  colors <- rainbow(n_clusters)
  for (i in 1:n_clusters) {
    cluster <- unique_clusters[i]
    cluster_mean <- atoms[cluster]
    abline(v = cluster_mean, col = colors[i], lwd = 3, lty = 2)
  }

  # Add legend
  legend("topright",
         legend = paste("Cluster", unique_clusters,
                       sprintf("(μ=%.2f)", atoms[unique_clusters])),
         col = colors,
         lwd = 3,
         lty = 2,
         cex = 0.8)
}


# Function to create a comprehensive plot
comprehensive_plot <- function(y, results, process_idx = 1) {
  par(mfrow = c(2, 2))

  # Plot 1: Cluster assignments over observations
  plot_cluster_results(y, results, process_idx, iteration = "last")

  # Plot 2: Histogram with cluster means
  plot_histogram_clusters(y, results, process_idx, iteration = "last")

  # Plot 3: Number of clusters over iterations
  n_clusters_vec <- sapply(results$n_clusters, function(x) x[process_idx])
  plot(1:length(n_clusters_vec), n_clusters_vec,
       type = "l",
       col = "blue",
       lwd = 2,
       xlab = "MCMC Iteration",
       ylab = "Number of Clusters",
       main = sprintf("Process %d - Clusters Over Time", process_idx))
  abline(h = mean(n_clusters_vec), col = "red", lty = 2, lwd = 2)
  legend("topright",
         legend = c("N clusters", sprintf("Mean = %.1f", mean(n_clusters_vec))),
         col = c("blue", "red"),
         lty = c(1, 2),
         lwd = 2)

  # Plot 4: Sigma^2 trace
  plot(1:length(results$sigma2), results$sigma2,
       type = "l",
       col = "darkgreen",
       lwd = 2,
       xlab = "MCMC Iteration",
       ylab = expression(sigma^2),
       main = "Variance Parameter Trace")
  abline(h = mean(results$sigma2), col = "red", lty = 2, lwd = 2)

  par(mfrow = c(1, 1))
}


# Example usage:
# ----------------
# After running MCMC:
#
# source("beta_product_dpy.R")
# source("mcmc_sampler.R")
#
# # Generate data
# set.seed(42)
# y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
# y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
# y <- list(y1, y2)
#
# # Run sampler
# results <- gibbs_sampler_dpy(
#   y = y,
#   K = 25,
#   n_iter = 5000,
#   burn_in = 2500,
#   alpha_vec = c(0.25, 0.25),
#   theta_vec = c(5, 5),
#   alpha0 = 0.25,
#   theta0 = 5,
#   verbose = TRUE
# )
#
# # Load plotting functions
# source("plot_clusters.R")
#
# # Simple plots
# plot_cluster_results(y, results, process_idx = 1)
# plot_cluster_results(y, results, process_idx = 2)
#
# # Side-by-side
# plot_both_processes(y, results)
#
# # Comprehensive diagnostics
# comprehensive_plot(y, results, process_idx = 1)
# comprehensive_plot(y, results, process_idx = 2)
#
# # Save to PDF
# pdf("cluster_results.pdf", width = 10, height = 8)
# comprehensive_plot(y, results, process_idx = 1)
# comprehensive_plot(y, results, process_idx = 2)
# dev.off()

cat("Plotting functions loaded!\n")
cat("\nAvailable functions:\n")
cat("  plot_cluster_results(y, results, process_idx = 1)\n")
cat("  plot_both_processes(y, results)\n")
cat("  plot_histogram_clusters(y, results, process_idx = 1)\n")
cat("  comprehensive_plot(y, results, process_idx = 1)\n")
cat("\nSee comments at the end of this file for usage examples.\n")
