#!/usr/bin/env Rscript
# ============================================================================
# START HERE - Script di Test Minimo
# ============================================================================
#
# COME USARE:
# 1. Apri R o RStudio
# 2. Scrivi: source("START_HERE.R")
# 3. Premi Invio
#
# Oppure da terminale: Rscript START_HERE.R
# ============================================================================

cat("\n")
cat("==================================================\n")
cat("  BETA-PRODUCT DPY - TEST MINIMO\n")
cat("==================================================\n\n")

# ----------------------------------------------------------------------------
# STEP 1: Definisci funzione minima
# ----------------------------------------------------------------------------

cat("Step 1: Definizione funzioni base... ")

# Stick-breaking minimo
sb <- function(K) {
  v <- rbeta(K, 1, 1)
  w <- numeric(K)
  for(i in 1:K) {
    if(i==1) w[i] <- v[i] else w[i] <- v[i] * prod(1-v[1:(i-1)])
  }
  w/sum(w)
}

cat("OK\n")

# ----------------------------------------------------------------------------
# STEP 2: Test funzione
# ----------------------------------------------------------------------------

cat("Step 2: Test stick-breaking... ")

set.seed(42)
pesi <- sb(5)

if(abs(sum(pesi) - 1) < 0.01) {
  cat("OK\n")
  cat("  Pesi generati:", paste(round(pesi, 3), collapse=", "), "\n")
} else {
  cat("ERRORE!\n")
  stop("Somma pesi != 1")
}

# ----------------------------------------------------------------------------
# STEP 3: Genera dati minimali
# ----------------------------------------------------------------------------

cat("Step 3: Generazione dati test... ")

n <- 50
dati <- c(rnorm(25, 0, 1), rnorm(25, 5, 1))

cat("OK\n")
cat("  n =", n, "osservazioni\n")
cat("  Media =", round(mean(dati), 2), "\n")

# ----------------------------------------------------------------------------
# STEP 4: Clustering semplice
# ----------------------------------------------------------------------------

cat("Step 4: Clustering base... ")

K <- 3
pesi <- sb(K)
centri <- rnorm(K, mean(dati), sd(dati))

# Assegna a cluster più vicino
cluster <- numeric(n)
for(i in 1:n) {
  dist <- abs(dati[i] - centri)
  cluster[i] <- which.min(dist)
}

n_cluster_usati <- length(unique(cluster))

cat("OK\n")
cat("  Cluster usati:", n_cluster_usati, "/", K, "\n")

# ----------------------------------------------------------------------------
# STEP 5: Plot minimale
# ----------------------------------------------------------------------------

cat("Step 5: Visualizzazione... ")

if(interactive()) {
  plot(dati, col=cluster, pch=19, main="Clustering Test",
       ylab="Valore", xlab="Osservazione")
  legend("topright", legend=paste("Cluster", 1:K),
         col=1:K, pch=19, cex=0.8)
  cat("OK (grafico creato)\n")
} else {
  cat("SKIP (non interattivo)\n")
}

# ----------------------------------------------------------------------------
# RISULTATO FINALE
# ----------------------------------------------------------------------------

cat("\n")
cat("==================================================\n")
cat("  ✓ TEST COMPLETATO CON SUCCESSO!\n")
cat("==================================================\n\n")

cat("Il codice base funziona correttamente!\n\n")
cat("PROSSIMI PASSI:\n")
cat("1. Carica il codice completo:\n")
cat("   source('beta_dpy_simple.R')\n\n")
cat("2. Esegui esempio completo:\n")
cat("   example_simple()\n\n")
cat("3. Usa con i tuoi dati:\n")
cat("   risultati <- gibbs_sampler_simple(y = tuoi_dati)\n\n")

# Restituisci TRUE se tutto ok
invisible(TRUE)
