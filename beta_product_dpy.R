# Beta-Product Dependent Pitman-Yor Processes for Bayesian Inference
# Implementation based on Bassetti, Casarin, and Leisen (2014)
# Journal of Econometrics, Vol. 180(1), pages 49-72

# Required libraries
library(MASS)
library(MCMCpack)

#' Stick-Breaking Process for Pitman-Yor
#'
#' @param alpha Discount parameter (0 <= alpha < 1)
#' @param theta Strength parameter (theta > -alpha)
#' @param K Number of stick-breaking weights to generate
#' @return Vector of stick-breaking weights
stick_breaking_py <- function(alpha, theta, K) {
  if (alpha < 0 || alpha >= 1) {
    stop("alpha must be in [0, 1)")
  }
  if (theta <= -alpha) {
    stop("theta must be > -alpha")
  }

  v <- numeric(K)
  pi_weights <- numeric(K)

  # Generate stick-breaking weights
  for (k in 1:K) {
    # Beta distribution for stick-breaking
    v[k] <- rbeta(1, 1 - alpha, theta + k * alpha)

    # Compute weights
    if (k == 1) {
      pi_weights[k] <- v[k]
    } else {
      pi_weights[k] <- v[k] * prod(1 - v[1:(k-1)])
    }
  }

  return(pi_weights)
}


#' Beta-Product Construction for Dependent Random Measures
#'
#' @param V Matrix of beta random variables (K x d)
#' @param W Vector of weights from a common base measure
#' @return Matrix of dependent weights (K x d)
beta_product_construction <- function(V, W) {
  K <- nrow(V)
  d <- ncol(V)

  # Initialize weight matrix
  pi_matrix <- matrix(0, nrow = K, ncol = d)

  # Compute dependent weights for each dimension
  for (j in 1:d) {
    for (k in 1:K) {
      if (k == 1) {
        pi_matrix[k, j] <- W[k] * V[k, j]
      } else {
        pi_matrix[k, j] <- W[k] * V[k, j] * prod(1 - V[1:(k-1), j])
      }
    }
  }

  return(pi_matrix)
}


#' Generate Dependent Pitman-Yor Processes
#'
#' @param d Number of dependent processes
#' @param K Truncation level for stick-breaking
#' @param alpha_vec Vector of discount parameters for each process
#' @param theta_vec Vector of strength parameters for each process
#' @param alpha0 Discount parameter for base measure
#' @param theta0 Strength parameter for base measure
#' @return List containing dependent weights and atom locations
generate_dpy <- function(d, K, alpha_vec, theta_vec, alpha0, theta0) {
  # Generate base measure weights using Pitman-Yor
  W <- stick_breaking_py(alpha0, theta0, K)

  # Generate beta random variables for dependence structure
  V <- matrix(0, nrow = K, ncol = d)
  for (j in 1:d) {
    for (k in 1:K) {
      V[k, j] <- rbeta(1, 1 - alpha_vec[j], theta_vec[j] + k * alpha_vec[j])
    }
  }

  # Generate dependent weights via beta-product construction
  pi_matrix <- beta_product_construction(V, W)

  # Generate common atom locations
  # In practice, these would come from the base distribution G0
  atoms <- matrix(rnorm(K * d), nrow = K, ncol = d)

  return(list(
    weights = pi_matrix,
    atoms = atoms,
    base_weights = W,
    V = V
  ))
}


#' Mixture Model with Dependent Pitman-Yor Prior
#'
#' @param y List of d data vectors
#' @param K Truncation level
#' @param alpha_vec Discount parameters
#' @param theta_vec Strength parameters
#' @param alpha0 Base discount parameter
#' @param theta0 Base strength parameter
#' @param sigma2 Variance of the mixture components
#' @return List with mixture parameters
dpy_mixture_model <- function(y, K, alpha_vec, theta_vec, alpha0, theta0, sigma2) {
  d <- length(y)
  n <- sapply(y, length)

  # Generate DPY process
  dpy <- generate_dpy(d, K, alpha_vec, theta_vec, alpha0, theta0)

  # Allocate observations to mixture components
  allocations <- vector("list", d)
  for (j in 1:d) {
    allocations[[j]] <- sample(1:K, n[j], replace = TRUE, prob = dpy$weights[, j])
  }

  return(list(
    dpy = dpy,
    allocations = allocations,
    hyperparameters = list(
      alpha_vec = alpha_vec,
      theta_vec = theta_vec,
      alpha0 = alpha0,
      theta0 = theta0,
      sigma2 = sigma2
    )
  ))
}


#' Compute Predictive Distribution
#'
#' @param new_data Vector of new observations
#' @param dpy DPY object
#' @param j Index of the process (which time series)
#' @param sigma2 Variance parameter
#' @return Predictive density values
predictive_dpy <- function(new_data, dpy, j, sigma2) {
  K <- nrow(dpy$weights)
  n_new <- length(new_data)

  # Compute predictive density as mixture
  pred_density <- numeric(n_new)
  for (i in 1:n_new) {
    density_val <- 0
    for (k in 1:K) {
      # Normal mixture component
      density_val <- density_val +
        dpy$weights[k, j] * dnorm(new_data[i], mean = dpy$atoms[k, j], sd = sqrt(sigma2))
    }
    pred_density[i] <- density_val
  }

  return(pred_density)
}


#' Sample Cluster Allocations (Gibbs Step)
#'
#' @param y Data vector
#' @param atoms Vector of atom locations
#' @param weights Vector of weights
#' @param sigma2 Variance parameter
#' @return Vector of cluster allocations
sample_allocations <- function(y, atoms, weights, sigma2) {
  n <- length(y)
  K <- length(atoms)
  allocations <- numeric(n)

  for (i in 1:n) {
    # Compute probabilities for each cluster
    probs <- numeric(K)
    for (k in 1:K) {
      probs[k] <- weights[k] * dnorm(y[i], mean = atoms[k], sd = sqrt(sigma2))
    }
    probs <- probs / sum(probs)

    # Sample allocation
    allocations[i] <- sample(1:K, 1, prob = probs)
  }

  return(allocations)
}


#' Update Atom Locations (Gibbs Step)
#'
#' @param y Data vector
#' @param allocations Cluster allocations
#' @param K Number of clusters
#' @param sigma2 Data variance
#' @param mu0 Prior mean
#' @param tau2 Prior variance
#' @return Updated atom locations
update_atoms <- function(y, allocations, K, sigma2, mu0 = 0, tau2 = 10) {
  atoms <- numeric(K)

  for (k in 1:K) {
    # Get observations allocated to cluster k
    y_k <- y[allocations == k]
    n_k <- length(y_k)

    if (n_k > 0) {
      # Posterior parameters for normal-normal conjugate model
      tau2_post <- 1 / (1/tau2 + n_k/sigma2)
      mu_post <- tau2_post * (mu0/tau2 + sum(y_k)/sigma2)

      # Sample from posterior
      atoms[k] <- rnorm(1, mean = mu_post, sd = sqrt(tau2_post))
    } else {
      # Sample from prior if no observations
      atoms[k] <- rnorm(1, mean = mu0, sd = sqrt(tau2))
    }
  }

  return(atoms)
}


#' Log-Likelihood for DPY Mixture Model
#'
#' @param y List of data vectors
#' @param dpy DPY object
#' @param allocations List of allocation vectors
#' @param sigma2 Variance parameter
#' @return Log-likelihood value
log_likelihood_dpy <- function(y, dpy, allocations, sigma2) {
  d <- length(y)
  ll <- 0

  for (j in 1:d) {
    n_j <- length(y[[j]])
    for (i in 1:n_j) {
      k <- allocations[[j]][i]
      ll <- ll + dnorm(y[[j]][i], mean = dpy$atoms[k, j], sd = sqrt(sigma2), log = TRUE)
    }
  }

  return(ll)
}


#' Compute Effective Sample Size (ESS)
#'
#' @param weights Vector or matrix of weights
#' @return Effective sample size
compute_ess <- function(weights) {
  if (is.matrix(weights)) {
    ess <- apply(weights, 2, function(w) 1 / sum(w^2))
  } else {
    ess <- 1 / sum(weights^2)
  }
  return(ess)
}


#' Compute Number of Active Clusters
#'
#' @param allocations Vector of cluster allocations
#' @return Number of unique clusters
count_clusters <- function(allocations) {
  return(length(unique(allocations)))
}


#' Compute Pairwise Correlation of Weights
#'
#' @param weight_matrix Matrix of weights (K x d)
#' @return Correlation matrix (d x d)
weight_correlation <- function(weight_matrix) {
  d <- ncol(weight_matrix)
  cor_matrix <- cor(weight_matrix)
  return(cor_matrix)
}
