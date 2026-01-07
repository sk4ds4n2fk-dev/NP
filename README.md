# Beta-Product Dependent Pitman-Yor Processes

R implementation of the Beta-Product Dependent Pitman-Yor (DPY) processes from:

**Bassetti, F., Casarin, R., & Leisen, F. (2014)**
*Beta-product dependent Pitman–Yor processes for Bayesian inference*
Journal of Econometrics, 180(1), 49-72
[SSRN Link](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2290959)

## Overview

This repository provides a complete implementation of Beta-Product Dependent Pitman-Yor processes for Bayesian non-parametric analysis of multiple dependent time series.

## Files

- `beta_product_dpy.R` - Core DPY implementation
- `mcmc_sampler.R` - MCMC algorithms for posterior inference
- `examples.R` - Example applications (mixture models, business cycles, VAR)
- `test_implementation.R` - Test suite
- `README_IMPLEMENTATION.md` - Detailed documentation

## Quick Start

```r
# Run tests
source("test_implementation.R")

# Run quick example
source("examples.R")
quick_test()

# Run all examples
run_all_examples()
```

## Requirements

```r
install.packages(c("MASS", "MCMCpack", "ggplot2"))
```

## Documentation

See `README_IMPLEMENTATION.md` for comprehensive documentation including:
- Mathematical background
- Function reference
- Usage examples
- Parameter tuning guide
- Troubleshooting

## Citation

```
Bassetti, F., Casarin, R., & Leisen, F. (2014).
Beta-product dependent Pitman–Yor processes for Bayesian inference.
Journal of Econometrics, 180(1), 49-72.
```