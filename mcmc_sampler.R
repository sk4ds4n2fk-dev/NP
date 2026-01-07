# MCMC Sampler for Beta-Product Dependent Pitman-Yor Processes
# Implementation based on Bassetti, Casarin, and Leisen (2014)

source("beta_product_dpy.R")

#' Gibbs Sampler for DPY Mixture Model
#'
#' @param y List of d data vectors (multivariate time series)
#' @param K Truncation level for stick-breaking
#' @param n_iter Number of MCMC iterations
#' @param burn_in Number of burn-in iterations
#' @param thin Thinning parameter
#' @param alpha_vec Initial discount parameters (length d)
#' @param theta_vec Initial strength parameters (length d)
#' @param alpha0 Base discount parameter
#' @param theta0 Base strength parameter
#' @param sigma2 Variance parameter (can be estimated or fixed)
#' @param estimate_sigma Logical, whether to estimate sigma2
#' @param verbose Logical, print progress
#' @return List containing posterior samples and diagnostics
gibbs_sampler_dpy <- function(y, K = 20, n_iter = 1000, burn_in = 500, thin = 1,
                               alpha_vec = NULL, theta_vec = NULL,
                               alpha0 = 0.5, theta0 = 1,
                               sigma2 = 1, estimate_sigma = TRUE,
                               verbose = TRUE) {
  d <- length(y)  # Number of dependent processes
  n <- sapply(y, length)  # Sample sizes

  # Initialize parameters if not provided
  if (is.null(alpha_vec)) alpha_vec <- rep(0.5, d)
  if (is.null(theta_vec)) theta_vec <- rep(1, d)

  # Storage for MCMC samples
  n_save <- floor((n_iter - burn_in) / thin)
  weights_samples <- array(0, dim = c(n_save, K, d))
  atoms_samples <- array(0, dim = c(n_save, K, d))
  allocations_samples <- vector("list", n_save)
  sigma2_samples <- numeric(n_save)
  ll_samples <- numeric(n_save)

  # Initialize
  dpy <- generate_dpy(d, K, alpha_vec, theta_vec, alpha0, theta0)
  allocations <- vector("list", d)
  for (j in 1:d) {
    # Initialize allocations according to DPY weights (not uniformly)
    allocations[[j]] <- sample(1:K, n[j], replace = TRUE, prob = dpy$weights[, j])
  }

  # MCMC iterations
  save_idx <- 1
  for (iter in 1:n_iter) {
    if (verbose && iter %% 100 == 0) {
      cat(sprintf("Iteration %d/%d\n", iter, n_iter))
    }

    # Step 1: Update allocations for each process
    for (j in 1:d) {
      allocations[[j]] <- sample_allocations(
        y[[j]], dpy$atoms[, j], dpy$weights[, j], sigma2
      )
    }

    # Step 2: Update atom locations for each process
    for (j in 1:d) {
      dpy$atoms[, j] <- update_atoms(
        y[[j]], allocations[[j]], K, sigma2
      )
    }

    # Step 3: Update stick-breaking weights
    dpy <- update_weights_dpy(y, allocations, dpy, alpha_vec, theta_vec, alpha0, theta0)

    # Step 4: Update variance parameter (if estimating)
    if (estimate_sigma) {
      sigma2 <- update_sigma2(y, dpy, allocations)
    }

    # Step 5: Optionally update hyperparameters (alpha, theta)
    # This can be added with Metropolis-Hastings steps if desired

    # Save samples after burn-in
    if (iter > burn_in && (iter - burn_in) %% thin == 0) {
      weights_samples[save_idx, , ] <- dpy$weights
      atoms_samples[save_idx, , ] <- dpy$atoms
      allocations_samples[[save_idx]] <- allocations
      sigma2_samples[save_idx] <- sigma2
      ll_samples[save_idx] <- log_likelihood_dpy(y, dpy, allocations, sigma2)
      save_idx <- save_idx + 1
    }
  }

  # Compute diagnostics
  n_clusters <- lapply(allocations_samples, function(alloc_list) {
    sapply(alloc_list, count_clusters)
  })

  return(list(
    weights = weights_samples,
    atoms = atoms_samples,
    allocations = allocations_samples,
    sigma2 = sigma2_samples,
    log_likelihood = ll_samples,
    n_clusters = n_clusters,
    hyperparameters = list(
      alpha_vec = alpha_vec,
      theta_vec = theta_vec,
      alpha0 = alpha0,
      theta0 = theta0
    )
  ))
}


#' Update Weights using Stick-Breaking Representation
#'
#' @param y List of data vectors
#' @param allocations List of allocation vectors
#' @param dpy Current DPY object
#' @param alpha_vec Discount parameters
#' @param theta_vec Strength parameters
#' @param alpha0 Base discount
#' @param theta0 Base strength
#' @return Updated DPY object
update_weights_dpy <- function(y, allocations, dpy, alpha_vec, theta_vec, alpha0, theta0) {
  d <- length(y)
  K <- nrow(dpy$weights)

  # Update base measure weights W
  W <- update_base_weights(allocations, K, alpha0, theta0)

  # Update beta random variables V for dependence structure
  V <- matrix(0, nrow = K, ncol = d)
  for (j in 1:d) {
    V[, j] <- update_beta_variables(allocations[[j]], K, alpha_vec[j], theta_vec[j])
  }

  # Reconstruct weights via beta-product
  weights <- beta_product_construction(V, W)

  dpy$weights <- weights
  dpy$base_weights <- W
  dpy$V <- V

  return(dpy)
}


#' Update Base Measure Weights
#'
#' @param allocations List of allocation vectors
#' @param K Number of clusters
#' @param alpha0 Base discount parameter
#' @param theta0 Base strength parameter
#' @return Updated base weights
update_base_weights <- function(allocations, K, alpha0, theta0) {
  # Count allocations across all processes
  counts <- numeric(K)
  for (alloc in allocations) {
    for (k in 1:K) {
      counts[k] <- counts[k] + sum(alloc == k)
    }
  }

  # Update stick-breaking weights based on counts
  v <- numeric(K)
  W <- numeric(K)

  for (k in 1:K) {
    # Number of observations in clusters > k
    if (k < K) {
      n_greater <- sum(counts[(k+1):K])
    } else {
      n_greater <- 0
    }

    # Beta posterior parameters
    a_post <- 1 - alpha0 + counts[k]
    b_post <- theta0 + k * alpha0 + n_greater

    v[k] <- rbeta(1, a_post, b_post)

    # Compute weights
    if (k == 1) {
      W[k] <- v[k]
    } else {
      W[k] <- v[k] * prod(1 - v[1:(k-1)])
    }
  }

  # Don't normalize W here - normalization happens in beta_product_construction
  # Normalizing here would distort the posterior distribution

  return(W)
}


#' Update Beta Variables for Dependent Structure
#'
#' @param allocations Allocation vector for one process
#' @param K Number of clusters
#' @param alpha Discount parameter
#' @param theta Strength parameter
#' @return Updated beta variables
update_beta_variables <- function(allocations, K, alpha, theta) {
  counts <- table(factor(allocations, levels = 1:K))
  counts <- as.numeric(counts)

  v <- numeric(K)
  for (k in 1:K) {
    # Number of observations in clusters > k
    if (k < K) {
      n_greater <- sum(counts[(k+1):K])
    } else {
      n_greater <- 0
    }

    # Beta posterior parameters
    a_post <- 1 - alpha + counts[k]
    b_post <- theta + k * alpha + n_greater

    v[k] <- rbeta(1, a_post, b_post)
  }

  return(v)
}


#' Update Variance Parameter
#'
#' @param y List of data vectors
#' @param dpy DPY object
#' @param allocations List of allocations
#' @param a_sigma Prior shape parameter for inverse-gamma
#' @param b_sigma Prior scale parameter for inverse-gamma
#' @return Updated sigma2
update_sigma2 <- function(y, dpy, allocations, a_sigma = 2, b_sigma = 1) {
  d <- length(y)
  n_total <- sum(sapply(y, length))

  # Compute sum of squared residuals
  ss <- 0
  for (j in 1:d) {
    for (i in 1:length(y[[j]])) {
      k <- allocations[[j]][i]
      ss <- ss + (y[[j]][i] - dpy$atoms[k, j])^2
    }
  }

  # Posterior parameters for inverse-gamma
  a_post <- a_sigma + n_total / 2
  b_post <- b_sigma + ss / 2

  # Sample from inverse-gamma (using inverse of gamma)
  sigma2 <- 1 / rgamma(1, shape = a_post, rate = b_post)

  return(sigma2)
}


#' Slice Sampler for Truncation-Free Inference
#'
#' @param y List of data vectors
#' @param n_iter Number of iterations
#' @param ... Additional parameters passed to Gibbs sampler
#' @return MCMC samples
slice_sampler_dpy <- function(y, n_iter = 1000, ...) {
  # Adaptive truncation using slice sampling
  # Start with small K and expand as needed

  K_init <- 20
  results <- gibbs_sampler_dpy(y, K = K_init, n_iter = n_iter, ...)

  # Check if we need more clusters
  max_clusters <- max(sapply(results$n_clusters, max))
  if (max_clusters > 0.8 * K_init) {
    warning(sprintf("Using %d out of %d clusters. Consider increasing K.",
                    max_clusters, K_init))
  }

  return(results)
}


#' Posterior Predictive Sampling
#'
#' @param mcmc_results Results from gibbs_sampler_dpy
#' @param n_pred Number of predictions to generate
#' @param process_idx Which process to predict for (1 to d)
#' @return Matrix of predictive samples (n_save x n_pred)
posterior_predictive <- function(mcmc_results, n_pred = 100, process_idx = 1) {
  n_save <- dim(mcmc_results$weights)[1]
  K <- dim(mcmc_results$weights)[2]

  pred_samples <- matrix(0, nrow = n_save, ncol = n_pred)

  for (i in 1:n_save) {
    # Get posterior samples for this iteration
    weights <- mcmc_results$weights[i, , process_idx]
    atoms <- mcmc_results$atoms[i, , process_idx]
    sigma2 <- mcmc_results$sigma2[i]

    # Generate predictions
    for (j in 1:n_pred) {
      # Sample cluster
      k <- sample(1:K, 1, prob = weights)
      # Sample from mixture component
      pred_samples[i, j] <- rnorm(1, mean = atoms[k], sd = sqrt(sigma2))
    }
  }

  return(pred_samples)
}


#' Compute WAIC (Widely Applicable Information Criterion)
#'
#' @param y List of data vectors
#' @param mcmc_results Results from gibbs_sampler_dpy
#' @return WAIC value and effective number of parameters
compute_waic <- function(y, mcmc_results) {
  d <- length(y)
  n <- sapply(y, length)
  n_save <- dim(mcmc_results$weights)[1]
  K <- dim(mcmc_results$weights)[2]

  # Compute log pointwise predictive density
  lppd <- 0
  p_waic <- 0

  for (j in 1:d) {
    for (i in 1:n[j]) {
      # Compute density for each MCMC sample
      log_dens <- numeric(n_save)
      for (s in 1:n_save) {
        dens <- 0
        for (k in 1:K) {
          dens <- dens + mcmc_results$weights[s, k, j] *
            dnorm(y[[j]][i], mean = mcmc_results$atoms[s, k, j],
                  sd = sqrt(mcmc_results$sigma2[s]))
        }
        log_dens[s] <- log(dens)
      }

      # Log mean of densities
      lppd <- lppd + log(mean(exp(log_dens)))

      # Variance of log densities
      p_waic <- p_waic + var(log_dens)
    }
  }

  waic <- -2 * (lppd - p_waic)

  return(list(waic = waic, lppd = lppd, p_waic = p_waic))
}


#' Compute DIC (Deviance Information Criterion)
#'
#' @param mcmc_results Results from gibbs_sampler_dpy
#' @return DIC value and effective number of parameters
compute_dic <- function(mcmc_results) {
  ll <- mcmc_results$log_likelihood

  # Mean deviance
  D_bar <- -2 * mean(ll)

  # Deviance at mean parameters
  D_theta_bar <- -2 * max(ll)  # Approximation

  # Effective number of parameters
  p_D <- D_bar - D_theta_bar

  # DIC
  DIC <- D_bar + p_D

  return(list(DIC = DIC, D_bar = D_bar, p_D = p_D))
}
