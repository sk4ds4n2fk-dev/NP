# ==============================================================================
# DEPENDENT PY PROCESS MIXTURE MODEL
# Implementation based on Bassetti, Casarin, and Leisen (2014)
# "Beta-product dependent Pitman-Yor processes for Bayesian inference"
# Journal of Econometrics, Vol. 180(1), pages 49-72
# ==============================================================================

library(ggplot2)
library(parallel)
library(mvtnorm)
library(dplyr)
library(lubridate)
library(zoo)
library(gridExtra)
library(viridis)
library(reshape2)
library(tidyverse)


# ==============================================================================
#   Section 1: Utility Functions
# ==============================================================================

rinvgamma <- function(n, shape, rate) {

  shape <- pmax(shape, 0.001) ## Protection from numerically unstable values
  rate  <- pmax(rate, 0.001)

  # Generate gamma
  vals <- tryCatch({
    rgamma(n, shape = shape, rate = rate)
  }, warning = function(w) {
    rep(1.0, n) # Fallback in case of extreme warning
  })

  # Protection: avoid division by zero or very small values
  vals <- pmax(vals, 1e-10)

  # Return inverse
  return(1 / vals)
}

logit <- function(p) log(p / (1 - p))
inv_logit <- function(x) exp(x) / (1 + exp(x))


# Compute weights from stick-breaking variables (Equation 2-3 in paper)
weights <- function(v_vec) {
  K <- length(v_vec)
  w <- numeric(K)
  resto <- 1
  for(k in 1:K) {
    w[k] <- v_vec[k] * resto
    resto <- resto * (1 - v_vec[k])
  }
  return(w)
}

# Log-Posterior for stick breaking (Equation 25 - Q_k function)
# Q_k(v_k|D, ψ) from equation after (25)
log_Qk <- function(v_vec, A_k, B_k, psi, k, tol = 1e-9) {
  v0 <- v_vec[1]; v1 <- v_vec[2]; v2 <- v_vec[3]
  a1 <- psi[1]; a2 <- psi[2]; l <- psi[3]

  # Numerical protection
  if(any(v_vec <= tol | v_vec >= (1 - tol))) return(-Inf)

  # Components from Q_k formula:
  # v_{0k}^{-l+A_{1k}+A_{2k}} (1-v_{0k})^{α₁-1}
  term_v0 <- (-l + A_k[1] + A_k[2]) * log(v0) + (a1 - 1) * log(1 - v0)

  # v_{ik}^{A_{ik}+α₁-l} (1-v_{ik})^{α₂+lk-1}
  term_v1 <- (A_k[1] + a1 - l) * log(v1) + (a2 + l*k - 1) * log(1 - v1)
  term_v2 <- (A_k[2] + a1 - l) * log(v2) + (a2 + l*k - 1) * log(1 - v2)

  # Constraints: (1-v_0*v_i)^{B_{ik}}
  arg_c1 <- 1 - v0 * v1
  arg_c2 <- 1 - v0 * v2
  if(arg_c1 <= 0 || arg_c2 <= 0) return(-Inf)

  val <- term_v0 + term_v1 + term_v2 + B_k[1] * log(arg_c1) + B_k[2] * log(arg_c2)

  if(is.na(val)) return(-Inf)
  return(val)
}

# ==============================================================================
# Section 2: Hyperparameter Configuration
# ==============================================================================

initialize_hyperparameters <- function(...) {

  configs <- list(
      # Atom variance priors
      s_sq_i = c(25, 25),
      s0_sq = 25,
      lambda = 3,
      eps = 3,

      # Stick-breaking (DP concentration) priors (for Gamma prior on α₁, α₂)
      psi_z11 = 10,  # ζ₁₁ in equation (50)
      psi_z21 = 10,  # ζ₂₁ in equation (50)
      psi_z12 = 10,  # ζ₁₂ in equation (50)
      psi_z22 = 10,  # ζ₂₂ in equation (50)

      # Initial values
      psi_init = c(1, 1, 0.25), # (α₁, α₂, l)

      # MCMC parameters
      scale_prop = 0.4,
      kappa_sq = 0.3,
      adapt_rate = 1.1,
      target_acc_low = 0.23,
      target_acc_high = 0.44,
      thin = 5
    )
  return(configs)
}

# ==============================================================================
# Section 3: MCMC Update Functions
# ==============================================================================

# Update Stick-Breaking Variables via Metropolis-Hastings (Equation 27)
# P{V* ∈ dv*|ψ̃, D} ∝ ∏_{k∈D*} Q_k(v_k|D, ψ̃)dv_k
update_sticks_MH <- function(V, D, psi, scale_prop = 0.3) {
  K_curr <- length(V$V0)
  new_V0 <- V$V0
  new_Vik <- V$Vik
  accepted <- 0

  for(k in 1:K_curr) {
    # Sufficient statistics
    A_k <- c(sum(D[[1]] == k), sum(D[[2]] == k))
    B_k <- c(sum(D[[1]] > k), sum(D[[2]] > k))

    # Current state
    v_curr <- c(V$V0[k], V$Vik[1, k], V$Vik[2, k])

    # Logit transformation for MH proposal
    logit_curr <- qlogis(v_curr)
    logit_curr[is.infinite(logit_curr) & logit_curr > 0] <- 20
    logit_curr[is.infinite(logit_curr) & logit_curr < 0] <- -20

    # Proposal
    logit_prop <- rnorm(3, mean = logit_curr, sd = scale_prop)
    v_prop <- plogis(logit_prop)

    # MH ratio computation
    log_q_curr <- log_Qk(v_curr, A_k, B_k, psi, k)
    log_q_prop <- log_Qk(v_prop, A_k, B_k, psi, k)
    log_jac_curr <- sum(log(v_curr) + log(1 - v_curr))
    log_jac_prop <- sum(log(v_prop) + log(1 - v_prop))

    mh_ratio <- if(is.infinite(log_q_curr) && log_q_curr < 0) {
      Inf
    } else {
      (log_q_prop + log_jac_prop) - (log_q_curr + log_jac_curr)
    }

    if(is.na(mh_ratio)) mh_ratio <- -Inf

    # Accept/Reject
    if(log(runif(1)) < mh_ratio) {
      new_V0[k] <- v_prop[1]
      new_Vik[1, k] <- v_prop[2]
      new_Vik[2, k] <- v_prop[3]
      accepted <- accepted + 1
    }
  }

  result <- list(V0 = new_V0, Vik = new_Vik)
  attr(result, "acc_rate") <- accepted / K_curr
  return(result)
}

# Update mixture atoms (Appendix B.1, Equations 44-49)
# Full conditionals for DPY(ψ, G₀) mixtures of Gaussians
update_atoms <- function(Y, D, mu_ik, sig_ik_sq, mu0, sig0_sq, hypers, K) {
  new_mu_ik <- mu_ik
  new_sig_ik_sq <- sig_ik_sq
  new_mu0 <- mu0
  new_sig0_sq <- sig0_sq

  s_sq_i <- hypers$s_sq_i
  s0_sq <- hypers$s0_sq
  lam <- hypers$lambda
  eps <- hypers$eps

  # Numerical protection
  MIN_VAR <- 1e-6

  for(k in 1:K) {
    # --- Component-specific part (Equations 44-45) ---
    for(i in 1:2) {
      idx <- which(D[[i]] == k)
      y_dat <- Y[[i]][idx]
      n_ik <- length(y_dat)

      # Use current variance with minimum protection
      curr_sig0 <- pmax(sig0_sq[k], MIN_VAR)
      curr_sigi <- pmax(sig_ik_sq[i, k], MIN_VAR)

      # Update μ_ik (Equation 44 - Normal posterior)
      prec_prior <- 1 / s_sq_i[i]
      prec_lik <- if(n_ik > 0) n_ik / (curr_sig0 * curr_sigi) else 0

      var_post <- 1 / (prec_prior + prec_lik)

      if(is.na(var_post) || var_post <= 0) var_post <- s_sq_i[i]

      if(n_ik > 0) {
        mean_post <- var_post * (sum(y_dat - mu0[k]) / (curr_sig0 * curr_sigi))
      } else {
        mean_post <- 0
      }

      if(is.na(mean_post)) mean_post <- 0

      new_mu_ik[i, k] <- rnorm(1, mean_post, sqrt(var_post))

      # Update σ²_ik (Equation 45 - Inverse Gamma posterior)
      eta2 <- if(n_ik > 0) sum((y_dat - (mu0[k] + new_mu_ik[i, k]))^2) else 0

      shape_p <- lam/2 + n_ik/2
      rate_p  <- lam/2 + eta2/(2*curr_sig0)

      new_sig_ik_sq[i, k] <- rinvgamma(1, shape = shape_p, rate = rate_p)
    }

    # --- Common part (Equations 47-49) ---
    prec_p0 <- 1 / s0_sq
    prec_lik_sum <- 0
    w_mean <- 0

    for(i in 1:2) {
      idx <- which(D[[i]] == k)
      y_dat <- Y[[i]][idx]
      n_ik <- length(y_dat)

      curr_sigi <- pmax(new_sig_ik_sq[i, k], MIN_VAR)
      curr_sig0 <- pmax(sig0_sq[k], MIN_VAR)

      if(n_ik > 0) {
        prec_lik_sum <- prec_lik_sum + n_ik / (curr_sig0 * curr_sigi)
        w_mean <- w_mean + sum(y_dat - new_mu_ik[i, k]) / (curr_sig0 * curr_sigi)
      }
    }

    var_p0 <- 1 / (prec_p0 + prec_lik_sum)
    if(is.na(var_p0) || var_p0 <= 0) var_p0 <- s0_sq

    mean_p0 <- var_p0 * w_mean
    if(is.na(mean_p0)) mean_p0 <- 0

    new_mu0[k] <- rnorm(1, mean_p0, sqrt(var_p0))

    # Update σ²_0k (Equation 49)
    rate_term <- 0
    nk_tot <- 0
    for(i in 1:2) {
      idx <- which(D[[i]] == k)
      y_dat <- Y[[i]][idx]
      nk_tot <- nk_tot + length(y_dat)

      curr_sigi <- pmax(new_sig_ik_sq[i, k], MIN_VAR)

      if(length(y_dat) > 0) {
        rate_term <- rate_term + sum((y_dat - (new_mu0[k] + new_mu_ik[i, k]))^2) / curr_sigi
      }
    }

    new_sig0_sq[k] <- rinvgamma(1,
                                shape = eps/2 + nk_tot/2,
                                rate = eps/2 + rate_term/2)
  }

  return(list(
    mu_ik = new_mu_ik,
    sig_ik_sq = new_sig_ik_sq,
    mu0 = new_mu0,
    sig0_sq = new_sig0_sq
  ))
}

# Update ψ = (α₁, α₂, l) (Equations 28, 50-51)
# P{ψ̃ ∈ (dα₁, dα₂, dl)|V*, D} from Equation (28) and (50)
update_psi <- function(psi, V, priors, kappa_sq) {

  # CORRECTED: This function computes log P(V|ψ) using PRIOR densities
  # from Equation (41), NOT the full conditional
  get_log_lik_psi_V <- function(p, v_obj) {
    a1 <- p[1]      # α₁
    a2 <- p[2]      # α₂
    alp <- p[3]     # l

    K_curr <- length(v_obj$V0)
    k_seq <- 1:K_curr

    # CORRECTED: Prior for V₀ (Equation 41)
    # P{V_{0k} ∈ dv|ψ} = v^{-l}(1-v)^{α₁-1} / B(1-l, α₁)
    log_prob_v0 <- -alp * log(v_obj$V0) +
      (a1 - 1) * log(1 - v_obj$V0) -  # FIXED: Changed from a1 to (a1 - 1)
      lbeta(1 - alp, a1)

    # CORRECTED: Prior for V_{ik} (Equation 41)
    # P{V_{ik} ∈ dv|k, ψ} = v^{α₁-l}(1-v)^{α₂+lk-1} / B(α₁+1-l, α₂+lk)
    log_norm_vi <- lbeta(a1 + 1 - alp, a2 + alp * k_seq)

    log_prob_v1 <- (a1 - alp) * log(v_obj$Vik[1,]) +
      (a2 + alp * k_seq - 1) * log(1 - v_obj$Vik[1,]) -  # FIXED: Added -1
      log_norm_vi

    log_prob_v2 <- (a1 - alp) * log(v_obj$Vik[2,]) +
      (a2 + alp * k_seq - 1) * log(1 - v_obj$Vik[2,]) -  # FIXED: Added -1
      log_norm_vi

    sum_val <- sum(log_prob_v0) + sum(log_prob_v1) + sum(log_prob_v2)
    if(is.na(sum_val) || is.nan(sum_val)) return(-Inf)
    return(sum_val)
  }

  # Proposal on transformed space (log for α's, logit for l)
  logit_l <- logit(psi[3])
  log_a1 <- log(psi[1])
  log_a2 <- log(psi[2])
  curr_trans <- c(log_a1, log_a2, logit_l)

  prop_trans <- as.vector(rmvnorm(1, mean = curr_trans, sigma = kappa_sq * diag(3)))
  psi_prop <- c(exp(prop_trans[1]), exp(prop_trans[2]), inv_logit(prop_trans[3]))

  # Compute Log-Posterior (Equation 50)

  # Log-Prior: Gamma(ζ₁₁, ζ₂₁) for α₁, Gamma(ζ₁₂, ζ₂₂) for α₂, Uniform for l
  lp_prior_curr <- (priors$z11 - 1) * log(psi[1]) - priors$z21 * psi[1] +
    (priors$z12 - 1) * log(psi[2]) - priors$z22 * psi[2]

  lp_prior_prop <- (priors$z11 - 1) * log(psi_prop[1]) - priors$z21 * psi_prop[1] +
    (priors$z12 - 1) * log(psi_prop[2]) - priors$z22 * psi_prop[2]

  # Log-Likelihood: P(V | ψ)
  ll_lik_curr <- get_log_lik_psi_V(psi, V)
  ll_lik_prop <- get_log_lik_psi_V(psi_prop, V)

  # Jacobian for transformation
  jac_curr <- sum(log(psi[1]), log(psi[2]), log(psi[3]) + log(1 - psi[3]))
  jac_prop <- sum(log(psi_prop[1]), log(psi_prop[2]), log(psi_prop[3]) + log(1 - psi_prop[3]))

  # Metropolis-Hastings Ratio
  log_prob_curr <- lp_prior_curr + ll_lik_curr + jac_curr
  log_prob_prop <- lp_prior_prop + ll_lik_prop + jac_prop

  mh_ratio <- log_prob_prop - log_prob_curr

  if(is.na(mh_ratio)) mh_ratio <- -Inf

  if(log(runif(1)) < mh_ratio) {
    return(psi_prop)
  } else {
    return(psi)
  }
}

# ==============================================================================
# Section 4: Main MCMC Sampler (Block Gibbs Sampler)
# Implementation of the algorithm described in Section 5
# ==============================================================================

run_DPY_sampler <- function(data,
                            hyperparams = initialize_hyperparameters(),
                            n_iter = 5000,
                            burnin = 2000,
                            adatta = 100,
                            thin = 5,
                            stampa = TRUE) {

  # Data setup
  Y <- list(data$Y1, data$Y2)
  N <- length(Y[[1]])

  # Hyperparameters
  hypers <- list(
    s_sq_i = hyperparams$s_sq_i,
    s0_sq = hyperparams$s0_sq,
    lambda = hyperparams$lambda,
    eps = hyperparams$eps
  )

  psi_priors <- list(
    z11 = hyperparams$psi_z11,
    z21 = hyperparams$psi_z21,
    z12 = hyperparams$psi_z12,
    z22 = hyperparams$psi_z22
  )

  # Initialize with enough components
  K <- 10

  # Initialize parameters
  psi <- hyperparams$psi_init
  range_y <- range(c(data$Y1, data$Y2))
  mu0 <- seq(range_y[1], range_y[2], length.out = K)
  sig0_sq <- rinvgamma(K, hypers$eps/2, hypers$eps/2)
  mu_ik <- matrix(rnorm(2*K, 0, sqrt(mean(hypers$s_sq_i))), 2, K)
  sig_ik_sq <- matrix(rinvgamma(2*K, hypers$lambda/2, hypers$lambda/2), 2, K)
  V <- list(V0 = rbeta(K, 1, 1), Vik = matrix(rbeta(2*K, 1, 1), 2, K))
  D <- list(sample(1:K, N, replace = TRUE), sample(1:K, N, replace = TRUE))

  # MCMC tuning parameters
  scale_prop <- hyperparams$scale_prop
  kappa_sq <- hyperparams$kappa_sq

  # Storage
  n_store <- floor((n_iter - burnin) / thin)
  trace_psi <- matrix(NA, n_iter, 3)
  trace_n_active <- matrix(NA, n_iter, 2)
  trace_atoms <- vector("list", n_store)
  trace_weights <- vector("list", n_store)
  trace_D <- vector("list", n_store)
  trace_acc_V <- numeric(n_iter)

  store_idx <- 1

  start_time <- Sys.time()

  # ===========================================================================
  # MCMC ITERATIONS
  # ===========================================================================

  for(iter in 1:n_iter) {

    # -------------------------------------------------------------------------
    # STEP 1: Update U (Equation 30)
    # P{U_{i,t} ∈ du|V, D} = 𝕀{u ≤ W_{i,D_{i,t}}} / W_{i,D_{i,t}} du
    # -------------------------------------------------------------------------

    W <- rbind(weights(V$Vik[1,]), weights(V$Vik[2,]))
    U <- list(runif(N, 0, W[1, D[[1]]]), runif(N, 0, W[2, D[[2]]]))

    # -------------------------------------------------------------------------
    # STEP 2: Adaptively expand K if needed (Slice Sampling - Walker 2007)
    # This implements the truncation scheme from Section 5.2
    # -------------------------------------------------------------------------

    min_u <- min(unlist(U))
    curr_resid <- 1 - rowSums(W)
    safety <- 0

    while(any(curr_resid > min_u) && safety < 50) {
      K <- K + 1
      safety <- safety + 1

      # Add new component
      mu0 <- c(mu0, rnorm(1, 0, sqrt(hypers$s0_sq)))
      sig0_sq <- c(sig0_sq, rinvgamma(1, hypers$eps/2, hypers$eps/2))
      mu_ik <- cbind(mu_ik, rnorm(2, 0, sqrt(hypers$s_sq_i)))
      sig_ik_sq <- cbind(sig_ik_sq, rinvgamma(2, hypers$lambda/2, hypers$lambda/2))
      V$V0 <- c(V$V0, rbeta(1, 1 - psi[3], psi[1] + 1))
      V$Vik <- cbind(V$Vik, rbeta(2, psi[1] + 1 - psi[3], psi[2] + psi[3]*K + 1))

      W <- rbind(weights(V$Vik[1,]), weights(V$Vik[2,]))
      curr_resid <- 1 - rowSums(W)
    }

    # -------------------------------------------------------------------------
    # STEP 3: Update D (Equation 31)
    # P{D_{i,t} = d|ϑ̃, V, U, ψ̃, Y} ∝ 𝒦_t(Y_{i,t}|ϑ̃_{id}, Z_t) 𝕀{U_{i,t} ≤ W_{i,d}}
    # -------------------------------------------------------------------------

    for(i in 1:2) {
      for(t in 1:N) {
        # Only consider clusters where U_{i,t} ≤ W_{i,k}
        active <- which(W[i, 1:K] > U[[i]][t])
        if(length(active) == 0) active <- which.max(W[i, 1:K])

        # Compute likelihood for each active cluster
        m <- mu0[active] + mu_ik[i, active]
        s <- sqrt(sig0_sq[active] * sig_ik_sq[i, active])

        lik <- dnorm(Y[[i]][t], m, s)

        lik[is.na(lik) | is.nan(lik)] <- 0

        if(sum(lik) == 0) {
          lik <- rep(1, length(lik))
        }

        # Sample allocation proportional to likelihood
        D[[i]][t] <- if(length(active) == 1) {
          active
        } else {
          active[sample.int(length(active), 1, prob = lik)]
        }
      }
    }

    # -------------------------------------------------------------------------
    # STEP 4: Update V via MH (Equation 27)
    # P{V* ∈ dv*|ψ̃, D} ∝ ∏_{k∈𝒟*} Q_k(v_k|D, ψ̃)dv_k
    # -------------------------------------------------------------------------

    V <- update_sticks_MH(V, D, psi, scale_prop)
    trace_acc_V[iter] <- attr(V, "acc_rate")

    # -------------------------------------------------------------------------
    # STEP 5: Update atoms ϑ̃ (Equation 24, Appendix B.1)
    # Gibbs sampling for Gaussian parameters
    # -------------------------------------------------------------------------

    atoms <- update_atoms(Y, D, mu_ik, sig_ik_sq, mu0, sig0_sq, hypers, K)
    mu_ik <- atoms$mu_ik
    sig_ik_sq <- atoms$sig_ik_sq
    mu0 <- atoms$mu0
    sig0_sq <- atoms$sig0_sq

    # -------------------------------------------------------------------------
    # STEP 6: Update ψ via MH (Equation 28, 50-51)
    # P{ψ̃ ∈ (dα₁, dα₂, dl)|V*, D}
    # -------------------------------------------------------------------------

    psi <- update_psi(psi, V, psi_priors, kappa_sq)

    # -------------------------------------------------------------------------
    # Adaptive tuning (only during burn-in)
    # -------------------------------------------------------------------------

    if(iter %% adatta == 0 && iter < burnin) {
      acc_rate <- mean(trace_acc_V[max(1, iter - adatta + 1):iter], na.rm = TRUE)

      if(acc_rate > hyperparams$target_acc_high) {
        scale_prop <- scale_prop * hyperparams$adapt_rate
      }
      if(acc_rate < hyperparams$target_acc_low) {
        scale_prop <- scale_prop / hyperparams$adapt_rate
      }

      scale_prop <- max(0.1, min(scale_prop, 2.0))
    }

    # -------------------------------------------------------------------------
    # Store results
    # -------------------------------------------------------------------------

    trace_psi[iter,] <- psi
    trace_n_active[iter,] <- c(length(unique(D[[1]])), length(unique(D[[2]])))

    if(iter > burnin && (iter - burnin) %% thin == 0) {
      trace_atoms[[store_idx]] <- list(
        mu0 = mu0,
        mu_ik = mu_ik,
        sig0_sq = sig0_sq,
        sig_ik_sq = sig_ik_sq
      )
      trace_weights[[store_idx]] <- list(W1 = W[1,], W2 = W[2,])
      trace_D[[store_idx]] <- D
      store_idx <- store_idx + 1
    }

    # Progress reporting
    if(stampa && iter %% 500 == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      cat(sprintf("Iteration %d/%d (%.1f min elapsed)\n", iter, n_iter, elapsed))
    }
  }

  return(list(
    psi = trace_psi,
    history_atoms = trace_atoms,
    history_weights = trace_weights,
    history_D = trace_D,
    n_active = trace_n_active,
    acc_rate = trace_acc_V,
    hyperparams = hyperparams,
    data = data
  ))
}

# ==============================================================================
# Section 5: Post-Processing and Diagnostics
# ==============================================================================

plot_diagnostics <- function(res) {
  n_iter <- nrow(res$psi)

  # PSI trace plots
  df_psi <- data.frame(
    Iteration = rep(1:n_iter, 2),
    Value = c(res$psi[,1], res$psi[,2]),
    Parameter = rep(c("alpha_1", "alpha_2"), each = n_iter)
  )

  p1 <- ggplot(df_psi, aes(x = Iteration, y = Value, color = Parameter)) +
    geom_line(alpha = 0.7) +
    facet_wrap(~Parameter, scales = "free_y", ncol = 1) +
    theme_minimal() +
    labs(title = "DPY Concentration Parameters (psi)", y = "Value") +
    theme(legend.position = "none")

  # Number of active components
  df_K <- data.frame(
    Iteration = rep(1:n_iter, 2),
    K = c(res$n_active[,1], res$n_active[,2]),
    Group = rep(c("Group 1", "Group 2"), each = n_iter)
  )

  p2 <- ggplot(df_K, aes(x = Iteration, y = K, color = Group)) +
    geom_line(alpha = 0.7) +
    theme_minimal() +
    labs(title = "Number of Active Components", y = "K") +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e"))

  # Acceptance rate
  df_acc <- data.frame(
    Iteration = 1:n_iter,
    AccRate = res$acc_rate
  )

  p3 <- ggplot(df_acc, aes(x = Iteration, y = AccRate)) +
    geom_line(alpha = 0.7, color = "darkgreen") +
    geom_hline(yintercept = c(0.23, 0.44), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "MH Acceptance Rate (Sticks)", y = "Acceptance Rate") +
    ylim(0, 1)

  grid.arrange(p1, p2, p3, ncol = 1)
}

# Dahl's method for point estimate (minimize posterior expected loss)
dahl <- function(res) {
  history <- res$history_D
  n_saved <- length(history)
  valid_idx <- 1:n_saved
  n_samples <- length(valid_idx)

  # Get dimensions
  N1 <- length(history[[1]][[1]])
  N2 <- length(history[[1]][[2]])
  N_tot <- N1 + N2

  # Allocation matrix (Rows = Iterations, Columns = Observations)
  alloc_mat <- matrix(NA, nrow = n_samples, ncol = N_tot)

  cat("Extracting allocations...\n")
  for(i in 1:n_samples) {
    iter <- valid_idx[i]
    d_list <- history[[iter]]
    alloc_mat[i, ] <- c(d_list[[1]], d_list[[2]])
  }

  # Compute pairwise probability matrix
  pairwise_prob <- matrix(0, N_tot, N_tot)

  for(i in 1:n_samples) {
    cl <- alloc_mat[i, ]
    adj <- outer(cl, cl, "==")
    pairwise_prob <- pairwise_prob + adj
  }
  pairwise_prob <- pairwise_prob / n_samples

  # Find iteration minimizing squared distance (Dahl's method)
  min_score <- Inf
  best_iter_idx <- -1

  for(i in 1:n_samples) {
    cl <- alloc_mat[i, ]
    current_adj <- outer(cl, cl, "==")
    score <- sum((current_adj - pairwise_prob)^2)

    if(score < min_score) {
      min_score <- score
      best_iter_idx <- valid_idx[i]
    }
  }

  # Extract parameters from best iteration
  best_atoms <- res$history_atoms[[best_iter_idx]]
  best_weights <- res$history_weights[[best_iter_idx]]
  best_D <- res$history_D[[best_iter_idx]]

  return(list(
    iter_idx = best_iter_idx,
    atoms = best_atoms,
    weights = best_weights,
    D = best_D,
    hyperparams = res$hyperparams
  ))
}

# Plot density estimates
plot_density <- function(best_model, res_original, threshold = 0.02, title = "Best Partition (Dahl)") {
  Y1 <- na.omit(res_original$data$Y1)
  Y2 <- na.omit(res_original$data$Y2)
  atoms <- best_model$atoms
  weights <- best_model$weights
  valid_idx <- which(weights$W1 > threshold)

  all_y <- c(Y1, Y2)
  x_seq <- seq(min(all_y), max(all_y), length.out = 500)

  calc_mixture <- function(x_vals, w_vec, mu_g, mu_l, sig_g, sig_l) {
    dens <- numeric(length(x_vals))
    active <- which(w_vec > threshold)
    for(k in active) {
      m <- mu_g[k] + mu_l[k]
      s <- sqrt(sig_g[k] * sig_l[k])
      d <- dnorm(x_vals, m, s)
      d[is.na(d)] <- 0
      dens <- dens + w_vec[k] * d
    }
    return(dens)
  }

  dens1 <- calc_mixture(x_seq, weights$W1, atoms$mu0, atoms$mu_ik[1,], atoms$sig0_sq, atoms$sig_ik_sq[1,])
  dens2 <- calc_mixture(x_seq, weights$W2, atoms$mu0, atoms$mu_ik[2,], atoms$sig0_sq, atoms$sig_ik_sq[2,])

  df_hist <- rbind(data.frame(Val = Y1, Group = "Group 1"), data.frame(Val = Y2, Group = "Group 2"))
  df_lines <- rbind(data.frame(x = x_seq, y = dens1, Group = "Group 1"), data.frame(x = x_seq, y = dens2, Group = "Group 2"))

  p <- ggplot() +
    geom_histogram(data = df_hist, aes(x = Val, y = after_stat(density)), fill="gray80", color="white", bins=50) +
    geom_line(data = df_lines, aes(x = x, y = y, color = Group), size = 1.2) +
    facet_wrap(~Group, scales = "free", ncol = 1) +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
    labs(title = title)

  print(p)
}
