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

### 🚀 Test Minimo (30 secondi)
```r
# Apri R e scrivi:
source("START_HERE.R")
```

### 📊 Versione Semplificata (2 minuti)
```r
# Carica versione semplificata
source("beta_dpy_simple.R")

# Test rapido
test_quick()

# Esempio completo con grafici
example_simple()
```

### 🔬 Versione Completa
```r
# Prima installa i pacchetti
install.packages(c("MASS", "MCMCpack", "ggplot2"))

# Esegui test completi
source("test_implementation.R")

# Esegui esempi
source("examples.R")
quick_test()
run_all_examples()
```

## Requirements

Solo R base per la versione semplificata. Per la versione completa:

```r
install.packages(c("MASS", "MCMCpack", "ggplot2"))
```

## Guida in Italiano

Vedi **GUIDA_RAPIDA.md** per istruzioni dettagliate in italiano.

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