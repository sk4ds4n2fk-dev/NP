# Troubleshooting Guide for Beta-Product DPY Implementation

## Issue: Cluster Collapse (Sampler assigns all data to 1 cluster)

### Description
The MCMC sampler may collapse to a single cluster, assigning all observations to cluster 1. This is a known challenge with truncated stick-breaking and Pitman-Yor priors.

### Why This Happens
1. **Rich-get-richer property**: Pitman-Yor processes have a self-reinforcing mechanism where clusters with more observations become more likely to receive new observations
2. **Truncation effects**: With finite K, once a cluster dominates, the remaining probability mass is redistributed in a way that reinforces the dominant cluster
3. **Initialization sensitivity**: Poor initialization can lead the sampler into a local mode

### Solutions to Try

#### 1. Increase the Truncation Level K
```r
results <- gibbs_sampler_dpy(
  y = y,
  K = 30,  # Try 20-30 instead of 10
  n_iter = 2000,
  ...
)
```

#### 2. Adjust Hyperparameters
**Smaller alpha** (discount parameter) and **larger theta** (strength parameter) encourage more clusters:

```r
results <- gibbs_sampler_dpy(
  y = y,
  K = 20,
  alpha_vec = c(0.25, 0.25),  # Smaller alpha
  theta_vec = c(5, 5),          # Larger theta
  alpha0 = 0.25,
  theta0 = 5,
  ...
)
```

**Guideline:**
- `alpha` in [0, 0.5): Smaller values → more clusters
- `theta` > 0: Larger values → more clusters
- Try alpha=0.1-0.3 and theta=2-10

#### 3. Better Initialization
The sampler now initializes according to DPY weights, but you can also try initializing based on k-means:

```r
# In mcmc_sampler.R, modify initialization:
library(stats)
# For each process, use k-means for initialization
for (j in 1:d) {
  if (n[j] >= K) {
    km <- kmeans(matrix(y[[j]], ncol=1), centers=min(K, 10))
    allocations[[j]] <- km$cluster
  } else {
    allocations[[j]] <- sample(1:K, n[j], replace=TRUE, prob=dpy$weights[,j])
  }
}
```

#### 4. Run Multiple Chains
Run multiple chains with different random seeds and check for convergence:

```r
set.seed(123)
results1 <- gibbs_sampler_dpy(...)

set.seed(456)
results2 <- gibbs_sampler_dpy(...)

set.seed(789)
results3 <- gibbs_sampler_dpy(...)

# Compare results
mean(sapply(results1$n_clusters, mean))
mean(sapply(results2$n_clusters, mean))
mean(sapply(results3$n_clusters, mean))
```

#### 5. Use the Debug Script
Run the debug script to see what's happening:

```r
source("debug_sampler.R")
```

This will show you:
- Initial weight distribution
- How allocations change each iteration
- Whether atoms are being updated correctly

#### 6. Increase Burn-in and Iterations
Allow more time for the sampler to explore:

```r
results <- gibbs_sampler_dpy(
  y = y,
  K = 20,
  n_iter = 5000,   # More iterations
  burn_in = 2500,  # Longer burn-in
  thin = 5,
  ...
)
```

### Example: Complete Setup for Better Mixing

```r
# Generate test data
set.seed(42)
y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
y <- list(y1, y2)

# Run with improved parameters
results <- gibbs_sampler_dpy(
  y = y,
  K = 25,                      # More clusters
  n_iter = 5000,              # More iterations
  burn_in = 2500,             # Longer burn-in
  thin = 5,
  alpha_vec = c(0.25, 0.25),  # Smaller alpha
  theta_vec = c(5, 5),         # Larger theta
  alpha0 = 0.25,
  theta0 = 5,
  sigma2 = 1,
  estimate_sigma = TRUE,
  verbose = TRUE
)

# Check results
cat(sprintf("Mean clusters (Process 1): %.1f\n",
            mean(sapply(results$n_clusters, function(x) x[1]))))
cat(sprintf("Mean clusters (Process 2): %.1f\n",
            mean(sapply(results$n_clusters, function(x) x[2]))))

# Should see 2-3 clusters instead of 1
```

### Advanced: Implementing Slice Sampling
For better handling of truncation, consider implementing the slice sampler variant (see Kalli et al. 2011, "Slice sampling mixture models"). This adaptively adjusts K during sampling.

### References
- Bassetti, Casarin, and Leisen (2014). "Beta-Product dependent Pitman-Yor processes for Bayesian inference"
- Kalli, Griffin, and Walker (2011). "Slice sampling mixture models"
- Pitman and Yor (1997). "The two-parameter Poisson-Dirichlet distribution"
