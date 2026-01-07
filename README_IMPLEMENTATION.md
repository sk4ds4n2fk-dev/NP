# Beta-Product Dependent Pitman-Yor Processes for Bayesian Inference

## Implementation in R

This repository contains a complete R implementation of the Beta-Product Dependent Pitman-Yor (DPY) processes as described in:

**Bassetti, F., Casarin, R., & Leisen, F. (2014)**
*Beta-product dependent Pitman–Yor processes for Bayesian inference*
Journal of Econometrics, 180(1), 49-72

## Overview

The Beta-Product Dependent Pitman-Yor process is a Bayesian non-parametric model for analyzing multiple dependent time series. It extends the Pitman-Yor process to allow for dependence across multiple processes while maintaining clustering properties.

### Key Features

- **Stick-breaking representation** of Pitman-Yor processes
- **Beta-product construction** for inducing dependence across processes
- **MCMC algorithms** for posterior inference (Gibbs sampling)
- **Model selection** tools (WAIC, DIC)
- **Posterior predictive sampling**
- Applications to:
  - Mixture models
  - Business cycle analysis
  - Vector autoregressive models

## Files

- `beta_product_dpy.R` - Core implementation of Beta-Product DPY processes
- `mcmc_sampler.R` - MCMC algorithms for posterior inference
- `examples.R` - Example applications and demonstrations
- `README_IMPLEMENTATION.md` - This documentation file

## Installation

### Required R Packages

```r
install.packages(c("MASS", "MCMCpack", "ggplot2"))
```

### Loading the Code

```r
source("beta_product_dpy.R")
source("mcmc_sampler.R")
source("examples.R")
```

## Quick Start

### Example 1: Simple Two-Process Mixture

```r
# Load the implementation
source("mcmc_sampler.R")

# Generate synthetic data
set.seed(123)
y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
y <- list(y1, y2)

# Run MCMC sampler
results <- gibbs_sampler_dpy(
  y = y,
  K = 20,                    # Truncation level
  n_iter = 2000,             # Number of iterations
  burn_in = 1000,            # Burn-in period
  thin = 2,                  # Thinning
  alpha_vec = c(0.5, 0.5),   # Discount parameters
  theta_vec = c(1, 1),       # Strength parameters
  alpha0 = 0.5,              # Base discount
  theta0 = 1,                # Base strength
  sigma2 = 1,                # Initial variance
  estimate_sigma = TRUE,     # Estimate variance
  verbose = TRUE
)

# Examine results
summary(results$sigma2)
sapply(results$n_clusters, function(x) x[1])  # Number of clusters in process 1
```

### Example 2: Business Cycle Application

```r
# Run the business cycle example
source("examples.R")
bc_results <- example2_business_cycles()

# The example simulates US and EU GDP growth rates
# with regime switches (expansion/recession)
```

### Example 3: Quick Test

```r
# Run a quick test
source("examples.R")
test_results <- quick_test()
```

## Core Functions

### Stick-Breaking Process

```r
# Generate stick-breaking weights for Pitman-Yor process
weights <- stick_breaking_py(alpha = 0.5, theta = 1, K = 20)
```

### Beta-Product Construction

```r
# Create dependent random measures
V <- matrix(rbeta(20 * 3, 1, 2), nrow = 20, ncol = 3)
W <- stick_breaking_py(0.5, 1, 20)
dependent_weights <- beta_product_construction(V, W)
```

### Generate Dependent Pitman-Yor Process

```r
# Generate DPY for d=3 processes
dpy <- generate_dpy(
  d = 3,                    # Number of processes
  K = 20,                   # Truncation level
  alpha_vec = c(0.5, 0.5, 0.5),
  theta_vec = c(1, 1, 1),
  alpha0 = 0.5,
  theta0 = 1
)

# Access components
dpy$weights       # Dependent weights (K x d matrix)
dpy$atoms         # Atom locations (K x d matrix)
dpy$base_weights  # Base measure weights (length K)
```

### MCMC Sampling

```r
# Full Gibbs sampler
results <- gibbs_sampler_dpy(
  y,                        # List of d data vectors
  K = 20,                   # Truncation level
  n_iter = 2000,            # MCMC iterations
  burn_in = 1000,           # Burn-in
  thin = 2,                 # Thinning
  alpha_vec = rep(0.5, d),  # Discount parameters
  theta_vec = rep(1, d),    # Strength parameters
  estimate_sigma = TRUE     # Estimate variance
)

# Access posterior samples
results$weights           # Posterior weights
results$atoms            # Posterior atom locations
results$sigma2           # Posterior variance
results$log_likelihood   # Log-likelihood trace
results$n_clusters       # Number of clusters
```

### Model Comparison

```r
# Compute WAIC
waic <- compute_waic(y, results)
cat("WAIC:", waic$waic, "\n")

# Compute DIC
dic <- compute_dic(results)
cat("DIC:", dic$DIC, "\n")
```

### Posterior Predictive Sampling

```r
# Generate predictions for process 1
predictions <- posterior_predictive(
  mcmc_results = results,
  n_pred = 100,
  process_idx = 1
)

# Posterior mean prediction
mean_pred <- colMeans(predictions)

# 95% credible intervals
lower <- apply(predictions, 2, quantile, 0.025)
upper <- apply(predictions, 2, quantile, 0.975)
```

## Mathematical Background

### Pitman-Yor Process

A Pitman-Yor process PY(α, θ, G₀) is characterized by:
- Discount parameter: 0 ≤ α < 1
- Strength parameter: θ > -α
- Base measure: G₀

### Stick-Breaking Representation

The weights are constructed as:
- vₖ ~ Beta(1 - α, θ + kα)
- πₖ = vₖ ∏ⱼ₌₁^(k-1) (1 - vⱼ)

### Beta-Product Construction

For dependent processes, the weights are:
- πₖʲ = wₖ × vₖʲ × ∏ₗ₌₁^(k-1) (1 - vₗʲ)

where wₖ are common base weights and vₖʲ are process-specific beta variables.

## Model Parameters

### Hyperparameters

- **α (alpha)**: Discount parameter, controls the tail behavior
  - Smaller α → more clusters
  - Larger α → fewer, larger clusters
  - Range: [0, 1)

- **θ (theta)**: Strength parameter, controls concentration
  - Smaller θ → fewer clusters
  - Larger θ → more clusters
  - Range: (−α, ∞)

- **K**: Truncation level for stick-breaking
  - Practical values: 10-50
  - Should be large enough that unused clusters have negligible weight

- **σ² (sigma2)**: Variance of mixture components
  - Can be fixed or estimated from data
  - Use inverse-gamma prior if estimating

### Choosing Hyperparameters

1. **For clustering applications**:
   - α ∈ [0.3, 0.7]
   - θ ∈ [0.5, 2]
   - K = 20 (adjust based on data complexity)

2. **For business cycle applications**:
   - α ∈ [0.2, 0.5] (fewer, persistent regimes)
   - θ ∈ [1, 3]
   - K = 10-15

3. **General advice**:
   - Start with α = 0.5, θ = 1
   - Increase K if many clusters are used
   - Use model selection (WAIC, DIC) to compare

## MCMC Diagnostics

### Convergence Checks

```r
# Trace plots
plot(results$sigma2, type = "l", main = "Trace of sigma2")

# Effective sample size for weights
ess <- compute_ess(results$weights[, , 1])  # For process 1
mean(ess)

# Number of clusters over iterations
n_clust <- sapply(results$n_clusters, function(x) x[1])
plot(n_clust, type = "l", main = "Number of Clusters")
```

### Recommendations

- **Burn-in**: At least 50% of iterations
- **Thinning**: Every 2-5 iterations to reduce autocorrelation
- **Total iterations**: 2000-5000 for most applications
- **Multiple chains**: Run 2-3 chains with different initializations

## Applications

### 1. Mixture Modeling

Use DPY for flexible density estimation with dependent processes:

```r
# Example: Clustering multiple related datasets
results <- gibbs_sampler_dpy(y, K = 20, ...)
final_clusters <- results$allocations[[length(results$allocations)]]
```

### 2. Business Cycle Analysis

Identify regimes (expansion/recession) in economic time series:

```r
# US and EU GDP growth rates
bc_results <- example2_business_cycles()

# Estimated regimes are given by cluster allocations
# Clusters with positive mean = expansion
# Clusters with negative mean = recession
```

### 3. Time Series with Regime Switching

Model multivariate time series with structural breaks:

```r
# VAR model with DPY prior
var_results <- example3_var_model()
```

## Computational Considerations

### Memory Usage

- Storing all MCMC samples can be memory-intensive
- For large K and d, consider:
  - Increasing thinning
  - Reducing number of saved iterations
  - Saving only summary statistics

### Speed Optimization

- The bottleneck is typically the allocation update step
- For large datasets (n > 1000), consider:
  - Reducing K if fewer clusters are needed
  - Parallelizing across processes
  - Using compiled code (Rcpp) for inner loops

### Numerical Stability

- Log-space computations for probability calculations
- Normalization of weights to prevent underflow
- Check for degenerate allocations (all data in one cluster)

## Extensions

### Possible Enhancements

1. **Hierarchical hyperparameters**: Estimate α and θ
2. **Time-varying processes**: Extend to non-stationary data
3. **Covariates**: Include predictors in the mixture
4. **Slice sampling**: Adaptive truncation without fixing K
5. **Alternative base measures**: Non-Gaussian components

## Troubleshooting

### Common Issues

**Issue**: Too many/few clusters identified

**Solution**: Adjust α and θ hyperparameters
- Too many → decrease θ or increase α
- Too few → increase θ or decrease α

---

**Issue**: Slow mixing of MCMC chains

**Solution**:
- Increase number of iterations
- Check for label switching
- Use better initialization

---

**Issue**: Memory errors

**Solution**:
- Reduce K (truncation level)
- Increase thinning
- Process data in batches

## References

1. **Bassetti, F., Casarin, R., & Leisen, F. (2014)**. Beta-product dependent Pitman–Yor processes for Bayesian inference. *Journal of Econometrics*, 180(1), 49-72.

2. **Pitman, J., & Yor, M. (1997)**. The two-parameter Poisson-Dirichlet distribution derived from a stable subordinator. *The Annals of Probability*, 855-900.

3. **Teh, Y. W. (2006)**. A hierarchical Bayesian language model based on Pitman-Yor processes. *Proceedings of the 21st International Conference on Computational Linguistics*.

4. **Griffin, J. E., & Walker, S. G. (2011)**. Posterior simulation of normalized random measure mixtures. *Journal of Computational and Graphical Statistics*, 20(1), 241-259.

## Citation

If you use this implementation, please cite:

```
Bassetti, F., Casarin, R., & Leisen, F. (2014).
Beta-product dependent Pitman–Yor processes for Bayesian inference.
Journal of Econometrics, 180(1), 49-72.
```

## License

This implementation is provided for educational and research purposes.

## Contact

For questions or issues with this implementation, please open an issue on the repository.

---

**Note**: This is an educational implementation. For production use, consider additional optimizations and validation.
