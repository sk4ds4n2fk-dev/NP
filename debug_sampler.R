# Debug script to understand the cluster collapse issue

source("beta_product_dpy.R")
source("mcmc_sampler.R")

set.seed(42)

# Generate simple test data with 2 clear clusters
d <- 2
n <- 100

y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
y <- list(y1, y2)

cat("Data summary:\n")
cat(sprintf("y1: mean=%.2f, sd=%.2f\n", mean(y1), sd(y1)))
cat(sprintf("y2: mean=%.2f, sd=%.2f\n", mean(y2), sd(y2)))
cat("\n")

# Run a few iterations manually to see what's happening
K <- 10
alpha_vec <- c(0.5, 0.5)
theta_vec <- c(1, 1)
alpha0 <- 0.5
theta0 <- 1
sigma2 <- 1

# Initialize
dpy <- generate_dpy(d, K, alpha_vec, theta_vec, alpha0, theta0)

cat("Initial weights (first 5 components):\n")
print(dpy$weights[1:5, ])
cat("\n")

# Initialize allocations
allocations <- vector("list", d)
for (j in 1:d) {
  allocations[[j]] <- sample(1:K, 100, replace = TRUE, prob = dpy$weights[, j])
}

cat("Initial allocation counts:\n")
cat("Process 1:", table(factor(allocations[[1]], levels=1:K)), "\n")
cat("Process 2:", table(factor(allocations[[2]], levels=1:K)), "\n")
cat("\n")

# Run 5 iterations
for (iter in 1:5) {
  cat(sprintf("===  Iteration %d ===\n", iter))

  # Update allocations
  for (j in 1:d) {
    allocations[[j]] <- sample_allocations(
      y[[j]], dpy$atoms[, j], dpy$weights[, j], sigma2
    )
  }

  cat("Allocation counts after sampling:\n")
  cat("Process 1:", table(factor(allocations[[1]], levels=1:K)), "\n")
  cat("Process 2:", table(factor(allocations[[2]], levels=1:K)), "\n")

  # Update atoms
  for (j in 1:d) {
    dpy$atoms[, j] <- update_atoms(y[[j]], allocations[[j]], K, sigma2)
  }

  cat("Atom means (first 3):\n")
  print(dpy$atoms[1:3, ])

  # Update weights
  dpy <- update_weights_dpy(y, allocations, dpy, alpha_vec, theta_vec, alpha0, theta0)

  cat("Weights (first 5 components):\n")
  print(dpy$weights[1:5, ])
  cat("\n")
}
