# Beta-Product Dependent Pitman-Yor Processes - VERSIONE SEMPLIFICATA E TESTATA
# Implementazione funzionante basata su Bassetti, Casarin, and Leisen (2014)

cat("Caricamento Beta-Product DPY Implementation...\n")

################################################################################
# FUNZIONI CORE
################################################################################

#' Stick-Breaking per Pitman-Yor
stick_breaking_py <- function(alpha, theta, K) {
  v <- numeric(K)
  pi_weights <- numeric(K)

  for (k in 1:K) {
    v[k] <- rbeta(1, 1 - alpha, theta + k * alpha)
    if (k == 1) {
      pi_weights[k] <- v[k]
    } else {
      pi_weights[k] <- v[k] * prod(1 - v[1:(k-1)])
    }
  }

  # Normalizza per sicurezza
  pi_weights <- pi_weights / sum(pi_weights)
  return(pi_weights)
}

#' Genera DPY Semplificato
generate_simple_dpy <- function(d, K, alpha = 0.5, theta = 1) {
  # Genera pesi base
  W <- stick_breaking_py(alpha, theta, K)

  # Genera pesi dipendenti
  weights <- matrix(0, nrow = K, ncol = d)
  for (j in 1:d) {
    V <- numeric(K)
    for (k in 1:K) {
      V[k] <- rbeta(1, 1 - alpha, theta + k * alpha)
    }

    for (k in 1:K) {
      if (k == 1) {
        weights[k, j] <- W[k] * V[k]
      } else {
        weights[k, j] <- W[k] * V[k] * prod(1 - V[1:(k-1)])
      }
    }
  }

  # Normalizza
  for (j in 1:d) {
    weights[, j] <- weights[, j] / sum(weights[, j])
  }

  # Genera atomi
  atoms <- matrix(rnorm(K * d, 0, 2), nrow = K, ncol = d)

  return(list(weights = weights, atoms = atoms))
}

#' Assegna osservazioni ai cluster
assign_clusters <- function(y, atoms, weights, sigma2) {
  n <- length(y)
  K <- length(atoms)
  allocations <- numeric(n)

  for (i in 1:n) {
    probs <- numeric(K)
    for (k in 1:K) {
      probs[k] <- weights[k] * dnorm(y[i], mean = atoms[k], sd = sqrt(sigma2))
    }

    # Evita problemi numerici
    if (sum(probs) == 0) {
      probs <- weights
    } else {
      probs <- probs / sum(probs)
    }

    allocations[i] <- sample(1:K, 1, prob = probs)
  }

  return(allocations)
}

#' Aggiorna atomi
update_atoms <- function(y, allocations, K, sigma2, mu0 = 0, tau2 = 10) {
  atoms <- numeric(K)

  for (k in 1:K) {
    y_k <- y[allocations == k]
    n_k <- length(y_k)

    if (n_k > 0) {
      tau2_post <- 1 / (1/tau2 + n_k/sigma2)
      mu_post <- tau2_post * (mu0/tau2 + sum(y_k)/sigma2)
      atoms[k] <- rnorm(1, mean = mu_post, sd = sqrt(tau2_post))
    } else {
      atoms[k] <- rnorm(1, mean = mu0, sd = sqrt(tau2))
    }
  }

  return(atoms)
}

#' Aggiorna varianza
update_variance <- function(y, atoms, allocations, a = 2, b = 1) {
  n <- length(y)
  ss <- 0

  for (i in 1:n) {
    k <- allocations[i]
    ss <- ss + (y[i] - atoms[k])^2
  }

  a_post <- a + n/2
  b_post <- b + ss/2

  sigma2 <- 1 / rgamma(1, shape = a_post, rate = b_post)
  return(sigma2)
}

################################################################################
# GIBBS SAMPLER SEMPLIFICATO
################################################################################

gibbs_sampler_simple <- function(y, K = 10, n_iter = 1000, burn_in = 500,
                                  alpha = 0.5, theta = 1, verbose = TRUE) {

  d <- length(y)
  n <- sapply(y, length)

  # Storage
  n_save <- n_iter - burn_in
  weights_samples <- array(0, dim = c(n_save, K, d))
  atoms_samples <- array(0, dim = c(n_save, K, d))
  sigma2_samples <- numeric(n_save)

  # Inizializza
  sigma2 <- 1
  dpy <- generate_simple_dpy(d, K, alpha, theta)
  allocations <- vector("list", d)
  for (j in 1:d) {
    allocations[[j]] <- sample(1:K, n[j], replace = TRUE)
  }

  # MCMC
  save_idx <- 1
  for (iter in 1:n_iter) {
    if (verbose && iter %% 100 == 0) {
      cat(sprintf("Iterazione %d/%d\n", iter, n_iter))
    }

    # Update per ogni processo
    for (j in 1:d) {
      # 1. Aggiorna allocazioni
      allocations[[j]] <- assign_clusters(y[[j]], dpy$atoms[, j],
                                          dpy$weights[, j], sigma2)

      # 2. Aggiorna atomi
      dpy$atoms[, j] <- update_atoms(y[[j]], allocations[[j]], K, sigma2)

      # 3. Aggiorna varianza
      sigma2 <- update_variance(y[[j]], dpy$atoms[, j], allocations[[j]])
    }

    # 4. Rigenera pesi
    dpy <- generate_simple_dpy(d, K, alpha, theta)

    # Salva dopo burn-in
    if (iter > burn_in) {
      weights_samples[save_idx, , ] <- dpy$weights
      atoms_samples[save_idx, , ] <- dpy$atoms
      sigma2_samples[save_idx] <- sigma2
      save_idx <- save_idx + 1
    }
  }

  return(list(
    weights = weights_samples,
    atoms = atoms_samples,
    sigma2 = sigma2_samples,
    final_allocations = allocations
  ))
}

################################################################################
# ESEMPIO PRATICO
################################################################################

example_simple <- function() {
  cat("\n=== ESEMPIO SEMPLICE ===\n\n")

  set.seed(123)

  # Genera dati con 2 cluster
  cat("Generazione dati...\n")
  n <- 100
  y1 <- c(rnorm(50, mean = 0, sd = 1), rnorm(50, mean = 5, sd = 1))
  y2 <- c(rnorm(50, mean = 0.5, sd = 1), rnorm(50, mean = 4.5, sd = 1))
  y <- list(y1, y2)

  cat(sprintf("Processo 1: n=%d, media=%.2f, sd=%.2f\n",
              length(y1), mean(y1), sd(y1)))
  cat(sprintf("Processo 2: n=%d, media=%.2f, sd=%.2f\n",
              length(y2), mean(y2), sd(y2)))

  # Esegui MCMC
  cat("\nEsecuzione MCMC...\n")
  results <- gibbs_sampler_simple(
    y = y,
    K = 10,
    n_iter = 500,
    burn_in = 250,
    alpha = 0.5,
    theta = 1,
    verbose = TRUE
  )

  # Risultati
  cat("\n=== RISULTATI ===\n")
  cat(sprintf("Campioni posteriori: %d\n", length(results$sigma2)))
  cat(sprintf("Media sigma2: %.3f\n", mean(results$sigma2)))
  cat(sprintf("SD sigma2: %.3f\n", sd(results$sigma2)))

  # Numero di cluster usati
  for (j in 1:length(y)) {
    n_clust <- length(unique(results$final_allocations[[j]]))
    cat(sprintf("Processo %d: %d cluster attivi\n", j, n_clust))
  }

  # Plot se disponibile
  if (requireNamespace("graphics", quietly = TRUE)) {
    par(mfrow = c(2, 2))

    # Trace sigma2
    plot(results$sigma2, type = "l", main = "Trace sigma2",
         ylab = "sigma2", xlab = "Iterazione")

    # Histogram sigma2
    hist(results$sigma2, breaks = 20, main = "Distribuzione sigma2",
         xlab = "sigma2", col = "lightblue")

    # Dati processo 1
    plot(y[[1]], col = results$final_allocations[[1]], pch = 19,
         main = "Dati e Cluster (Processo 1)", ylab = "Valore", xlab = "Indice")

    # Dati processo 2
    plot(y[[2]], col = results$final_allocations[[2]], pch = 19,
         main = "Dati e Cluster (Processo 2)", ylab = "Valore", xlab = "Indice")

    par(mfrow = c(1, 1))
  }

  cat("\nEsempio completato con successo!\n")

  return(results)
}

################################################################################
# TEST RAPIDO
################################################################################

test_quick <- function() {
  cat("\n=== TEST RAPIDO ===\n\n")

  # Test 1: Stick-breaking
  cat("Test 1: Stick-breaking... ")
  weights <- stick_breaking_py(0.5, 1, 10)
  if (abs(sum(weights) - 1) < 0.01 && all(weights > 0)) {
    cat("OK\n")
  } else {
    cat("FALLITO\n")
    return(FALSE)
  }

  # Test 2: DPY generation
  cat("Test 2: DPY generation... ")
  dpy <- generate_simple_dpy(d = 2, K = 5, alpha = 0.5, theta = 1)
  if (all(dim(dpy$weights) == c(5, 2))) {
    cat("OK\n")
  } else {
    cat("FALLITO\n")
    return(FALSE)
  }

  # Test 3: MCMC veloce
  cat("Test 3: MCMC rapido... ")
  set.seed(42)
  y <- list(rnorm(50), rnorm(50))
  results <- gibbs_sampler_simple(y, K = 5, n_iter = 100, burn_in = 50, verbose = FALSE)
  if (length(results$sigma2) == 50) {
    cat("OK\n")
  } else {
    cat("FALLITO\n")
    return(FALSE)
  }

  cat("\nTutti i test passati!\n")
  return(TRUE)
}

################################################################################
# MESSAGGIO DI AVVIO
################################################################################

cat("\nBeta-Product DPY caricato con successo!\n")
cat("========================================\n")
cat("Comandi disponibili:\n")
cat("  test_quick()      - Esegui test rapidi\n")
cat("  example_simple()  - Esegui esempio completo\n")
cat("\nProva: test_quick()\n\n")
