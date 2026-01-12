# ==============================================================================
# DIAGNOSTIC AND IMPROVED VERSION OF DPY MIXTURE MODEL
# ==============================================================================

source("dpy_mixture_sampler.R")

# ==============================================================================
# PROBLEMA DIAGNOSTICATO:
# ==============================================================================
# La densità è troppo piatta perché:
# 1. Le varianze del modello sono il PRODOTTO σ²_0k * σ²_{ik}, che diventa molto grande
# 2. I prior sono troppo vaghi (s_sq_i = 25, s0_sq = 25)
# 3. Il modello ha μ_total = μ_0k + μ_{ik}, con troppa flessibilità
# 4. Threshold potrebbe eliminare componenti importanti

# ==============================================================================
# SOLUZIONE 1: Hyperparametri più informativi
# ==============================================================================

initialize_hyperparameters_improved <- function() {
  configs <- list(
      # REDUCED variance priors for better concentration
      s_sq_i = c(5, 5),      # Ridotto da 25
      s0_sq = 5,             # Ridotto da 25
      lambda = 5,            # Aumentato da 3 per prior più informativo
      eps = 5,               # Aumentato da 3

      # Stick-breaking priors - aumentiamo per favorire più cluster
      psi_z11 = 5,           # Ridotto da 10 per permettere più flessibilità
      psi_z21 = 1,           # Ridotto da 10
      psi_z12 = 5,
      psi_z22 = 1,

      # Initial values
      psi_init = c(0.5, 2, 0.1), # Modificato: α₁ più piccolo, α₂ più grande, l più piccolo

      # MCMC parameters
      scale_prop = 0.3,      # Ridotto per più stabilità
      kappa_sq = 0.2,        # Ridotto
      adapt_rate = 1.1,
      target_acc_low = 0.23,
      target_acc_high = 0.44,
      thin = 5
    )
  return(configs)
}

# ==============================================================================
# SOLUZIONE 2: Versione semplificata del modello
# ==============================================================================
# Invece di μ_total = μ_0k + μ_{ik} con var = σ²_0k * σ²_{ik}
# Usiamo solo μ_k con var = σ²_k (più semplice e interpretabile)

update_atoms_simplified <- function(Y, D, mu_k, sig_k_sq, hypers, K) {
  new_mu_k <- mu_k
  new_sig_k_sq <- sig_k_sq

  s0_sq <- hypers$s0_sq
  lam <- hypers$lambda

  MIN_VAR <- 1e-6

  for(k in 1:K) {
    # Pooled update across both groups
    y_all <- c()
    for(i in 1:2) {
      idx <- which(D[[i]] == k)
      if(length(idx) > 0) {
        y_all <- c(y_all, Y[[i]][idx])
      }
    }

    n_k <- length(y_all)

    if(n_k > 0) {
      # Update μ_k
      prec_prior <- 1 / s0_sq
      prec_lik <- n_k / pmax(sig_k_sq[k], MIN_VAR)
      var_post <- 1 / (prec_prior + prec_lik)
      mean_post <- var_post * (sum(y_all) / pmax(sig_k_sq[k], MIN_VAR))

      new_mu_k[k] <- rnorm(1, mean_post, sqrt(var_post))

      # Update σ²_k
      ss <- sum((y_all - new_mu_k[k])^2)
      new_sig_k_sq[k] <- rinvgamma(1,
                                   shape = lam/2 + n_k/2,
                                   rate = lam/2 + ss/2)
    } else {
      # Sample from prior
      new_mu_k[k] <- rnorm(1, 0, sqrt(s0_sq))
      new_sig_k_sq[k] <- rinvgamma(1, lam/2, lam/2)
    }
  }

  return(list(mu_k = new_mu_k, sig_k_sq = new_sig_k_sq))
}

# ==============================================================================
# SOLUZIONE 3: Run con diagnostici migliori
# ==============================================================================

run_DPY_sampler_improved <- function(data,
                                     hyperparams = initialize_hyperparameters_improved(),
                                     use_simplified_model = TRUE,
                                     n_iter = 10000,
                                     burnin = 5000,
                                     thin = 5,
                                     stampa = TRUE) {

  Y <- list(data$Y1, data$Y2)
  N <- length(Y[[1]])

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

  K <- 20  # Start with more components
  psi <- hyperparams$psi_init
  range_y <- range(c(data$Y1, data$Y2))

  if(use_simplified_model) {
    # Simplified model: single μ_k and σ²_k per cluster
    mu_k <- seq(range_y[1], range_y[2], length.out = K)
    sig_k_sq <- rinvgamma(K, hypers$lambda/2, hypers$lambda/2)
  } else {
    # Original complex model
    mu0 <- seq(range_y[1], range_y[2], length.out = K)
    sig0_sq <- rinvgamma(K, hypers$eps/2, hypers$eps/2)
    mu_ik <- matrix(rnorm(2*K, 0, sqrt(mean(hypers$s_sq_i))), 2, K)
    sig_ik_sq <- matrix(rinvgamma(2*K, hypers$lambda/2, hypers$lambda/2), 2, K)
  }

  V <- list(V0 = rbeta(K, 1, 1), Vik = matrix(rbeta(2*K, 1, 1), 2, K))
  D <- list(sample(1:K, N, replace = TRUE), sample(1:K, N, replace = TRUE))

  scale_prop <- hyperparams$scale_prop
  kappa_sq <- hyperparams$kappa_sq

  n_store <- floor((n_iter - burnin) / thin)
  trace_psi <- matrix(NA, n_iter, 3)
  trace_n_active <- matrix(NA, n_iter, 2)
  trace_atoms <- vector("list", n_store)
  trace_weights <- vector("list", n_store)
  trace_D <- vector("list", n_store)
  trace_acc_V <- numeric(n_iter)
  if(use_simplified_model) {
    trace_variances <- matrix(NA, n_store, K)
  }

  store_idx <- 1
  start_time <- Sys.time()

  for(iter in 1:n_iter) {

    # STEP 1: Update U
    W <- rbind(weights(V$Vik[1,]), weights(V$Vik[2,]))
    U <- list(runif(N, 0, W[1, D[[1]]]), runif(N, 0, W[2, D[[2]]]))

    # STEP 2: Expand K if needed
    min_u <- min(unlist(U))
    curr_resid <- 1 - rowSums(W)
    safety <- 0

    while(any(curr_resid > min_u) && safety < 50) {
      K <- K + 1
      safety <- safety + 1

      if(use_simplified_model) {
        mu_k <- c(mu_k, rnorm(1, 0, sqrt(hypers$s0_sq)))
        sig_k_sq <- c(sig_k_sq, rinvgamma(1, hypers$lambda/2, hypers$lambda/2))
      } else {
        mu0 <- c(mu0, rnorm(1, 0, sqrt(hypers$s0_sq)))
        sig0_sq <- c(sig0_sq, rinvgamma(1, hypers$eps/2, hypers$eps/2))
        mu_ik <- cbind(mu_ik, rnorm(2, 0, sqrt(hypers$s_sq_i)))
        sig_ik_sq <- cbind(sig_ik_sq, rinvgamma(2, hypers$lambda/2, hypers$lambda/2))
      }

      V$V0 <- c(V$V0, rbeta(1, 1 - psi[3], psi[1] + 1))
      V$Vik <- cbind(V$Vik, rbeta(2, psi[1] + 1 - psi[3], psi[2] + psi[3]*K + 1))

      W <- rbind(weights(V$Vik[1,]), weights(V$Vik[2,]))
      curr_resid <- 1 - rowSums(W)
    }

    # STEP 3: Update D
    for(i in 1:2) {
      for(t in 1:N) {
        active <- which(W[i, 1:K] > U[[i]][t])
        if(length(active) == 0) active <- which.max(W[i, 1:K])

        if(use_simplified_model) {
          m <- mu_k[active]
          s <- sqrt(sig_k_sq[active])
        } else {
          m <- mu0[active] + mu_ik[i, active]
          s <- sqrt(sig0_sq[active] * sig_ik_sq[i, active])
        }

        lik <- dnorm(Y[[i]][t], m, s)
        lik[is.na(lik) | is.nan(lik)] <- 0
        if(sum(lik) == 0) lik <- rep(1, length(lik))

        D[[i]][t] <- if(length(active) == 1) {
          active
        } else {
          active[sample.int(length(active), 1, prob = lik)]
        }
      }
    }

    # STEP 4: Update V
    V <- update_sticks_MH(V, D, psi, scale_prop)
    trace_acc_V[iter] <- attr(V, "acc_rate")

    # STEP 5: Update atoms
    if(use_simplified_model) {
      atoms <- update_atoms_simplified(Y, D, mu_k, sig_k_sq, hypers, K)
      mu_k <- atoms$mu_k
      sig_k_sq <- atoms$sig_k_sq
    } else {
      atoms <- update_atoms(Y, D, mu_ik, sig_ik_sq, mu0, sig0_sq, hypers, K)
      mu_ik <- atoms$mu_ik
      sig_ik_sq <- atoms$sig_ik_sq
      mu0 <- atoms$mu0
      sig0_sq <- atoms$sig0_sq
    }

    # STEP 6: Update psi
    psi <- update_psi(psi, V, psi_priors, kappa_sq)

    # Adaptive tuning
    if(iter %% 100 == 0 && iter < burnin) {
      acc_rate <- mean(trace_acc_V[max(1, iter - 99):iter], na.rm = TRUE)

      if(acc_rate > hyperparams$target_acc_high) {
        scale_prop <- scale_prop * hyperparams$adapt_rate
      }
      if(acc_rate < hyperparams$target_acc_low) {
        scale_prop <- scale_prop / hyperparams$adapt_rate
      }

      scale_prop <- max(0.1, min(scale_prop, 2.0))
    }

    # Store results
    trace_psi[iter,] <- psi
    trace_n_active[iter,] <- c(length(unique(D[[1]])), length(unique(D[[2]])))

    if(iter > burnin && (iter - burnin) %% thin == 0) {
      if(use_simplified_model) {
        trace_atoms[[store_idx]] <- list(mu_k = mu_k, sig_k_sq = sig_k_sq)
        trace_variances[store_idx,] <- sig_k_sq[1:ncol(trace_variances)]
      } else {
        trace_atoms[[store_idx]] <- list(
          mu0 = mu0, mu_ik = mu_ik,
          sig0_sq = sig0_sq, sig_ik_sq = sig_ik_sq
        )
      }
      trace_weights[[store_idx]] <- list(W1 = W[1,], W2 = W[2,])
      trace_D[[store_idx]] <- D
      store_idx <- store_idx + 1
    }

    if(stampa && iter %% 500 == 0) {
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      if(use_simplified_model) {
        avg_var <- mean(sig_k_sq[1:min(10, K)])
        cat(sprintf("Iter %d/%d (%.1f min) | Avg var: %.2f | K: %d | Acc: %.2f\n",
                    iter, n_iter, elapsed, avg_var, K, trace_acc_V[iter]))
      } else {
        cat(sprintf("Iter %d/%d (%.1f min) | K: %d | Acc: %.2f\n",
                    iter, n_iter, elapsed, K, trace_acc_V[iter]))
      }
    }
  }

  result <- list(
    psi = trace_psi,
    history_atoms = trace_atoms,
    history_weights = trace_weights,
    history_D = trace_D,
    n_active = trace_n_active,
    acc_rate = trace_acc_V,
    hyperparams = hyperparams,
    data = data,
    simplified = use_simplified_model
  )

  if(use_simplified_model) {
    result$trace_variances <- trace_variances
  }

  return(result)
}

# ==============================================================================
# Plot density function migliorata
# ==============================================================================

plot_density_improved <- function(best_model, res_original, threshold = 0.005,
                                  title = "Best Partition (Dahl)") {
  Y1 <- na.omit(res_original$data$Y1)
  Y2 <- na.omit(res_original$data$Y2)
  atoms <- best_model$atoms
  weights <- best_model$weights

  all_y <- c(Y1, Y2)
  x_seq <- seq(min(all_y), max(all_y), length.out = 500)

  if(res_original$simplified) {
    # Simplified model
    calc_mixture <- function(x_vals, w_vec, mu, sig_sq) {
      dens <- numeric(length(x_vals))
      active <- which(w_vec > threshold)
      cat(sprintf("Using %d components (threshold=%.4f)\n", length(active), threshold))

      for(k in active) {
        m <- mu[k]
        s <- sqrt(sig_sq[k])
        d <- dnorm(x_vals, m, s)
        d[is.na(d)] <- 0
        dens <- dens + w_vec[k] * d
      }
      return(dens)
    }

    dens1 <- calc_mixture(x_seq, weights$W1, atoms$mu_k, atoms$sig_k_sq)
    dens2 <- calc_mixture(x_seq, weights$W2, atoms$mu_k, atoms$sig_k_sq)

  } else {
    # Original complex model
    calc_mixture <- function(x_vals, w_vec, mu_g, mu_l, sig_g, sig_l, group_idx) {
      dens <- numeric(length(x_vals))
      active <- which(w_vec > threshold)
      cat(sprintf("Group %d: Using %d components (threshold=%.4f)\n",
                  group_idx, length(active), threshold))

      for(k in active) {
        m <- mu_g[k] + mu_l[k]
        s <- sqrt(sig_g[k] * sig_l[k])
        d <- dnorm(x_vals, m, s)
        d[is.na(d)] <- 0
        dens <- dens + w_vec[k] * d
      }
      return(dens)
    }

    dens1 <- calc_mixture(x_seq, weights$W1, atoms$mu0, atoms$mu_ik[1,],
                         atoms$sig0_sq, atoms$sig_ik_sq[1,], 1)
    dens2 <- calc_mixture(x_seq, weights$W2, atoms$mu0, atoms$mu_ik[2,],
                         atoms$sig0_sq, atoms$sig_ik_sq[2,], 2)
  }

  df_hist <- rbind(
    data.frame(Val = Y1, Group = "Group 1"),
    data.frame(Val = Y2, Group = "Group 2")
  )
  df_lines <- rbind(
    data.frame(x = x_seq, y = dens1, Group = "Group 1"),
    data.frame(x = x_seq, y = dens2, Group = "Group 2")
  )

  p <- ggplot() +
    geom_histogram(data = df_hist, aes(x = Val, y = after_stat(density)),
                   fill="gray80", color="white", bins=50) +
    geom_line(data = df_lines, aes(x = x, y = y, color = Group), size = 1.2) +
    facet_wrap(~Group, scales = "free", ncol = 1) +
    theme_minimal(base_size = 14) +
    scale_color_manual(values = c("#1f77b4", "#ff7f0e")) +
    labs(title = title)

  print(p)
}

# ==============================================================================
# Funzione per diagnostica completa
# ==============================================================================

diagnose_results <- function(res) {
  cat("\n========== DIAGNOSTIC REPORT ==========\n\n")

  # 1. Convergence diagnostics
  cat("1. CONVERGENCE DIAGNOSTICS\n")
  cat("--------------------------\n")
  cat(sprintf("Acceptance rate (V): %.3f (target: 0.23-0.44)\n",
              mean(res$acc_rate, na.rm=TRUE)))

  # 2. Number of active clusters
  cat("\n2. NUMBER OF ACTIVE CLUSTERS\n")
  cat("----------------------------\n")
  cat(sprintf("Group 1 - Mean: %.1f, Median: %d\n",
              mean(res$n_active[,1]), median(res$n_active[,1])))
  cat(sprintf("Group 2 - Mean: %.1f, Median: %d\n",
              mean(res$n_active[,2]), median(res$n_active[,2])))

  # 3. Concentration parameters
  cat("\n3. CONCENTRATION PARAMETERS (PSI)\n")
  cat("---------------------------------\n")
  psi_mean <- colMeans(res$psi, na.rm=TRUE)
  cat(sprintf("alpha_1: %.3f\n", psi_mean[1]))
  cat(sprintf("alpha_2: %.3f\n", psi_mean[2]))
  cat(sprintf("l: %.3f\n", psi_mean[3]))

  # 4. Variance diagnostics (if simplified model)
  if(!is.null(res$trace_variances)) {
    cat("\n4. VARIANCE DIAGNOSTICS\n")
    cat("-----------------------\n")
    avg_vars <- colMeans(res$trace_variances, na.rm=TRUE)
    cat(sprintf("Mean variance (1st 5 comp): %.2f\n", mean(avg_vars[1:5])))
    cat(sprintf("Std variance (1st 5 comp): %.2f\n", sd(avg_vars[1:5])))
  }

  cat("\n=======================================\n")
}
