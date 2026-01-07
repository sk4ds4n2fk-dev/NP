# COVID-19 Case Data Analysis with Beta-Product DPY
# This script shows how to properly analyze COVID case count data

source("beta_product_dpy.R")
source("mcmc_sampler.R")
source("plot_clusters.R")

set.seed(42)

# --- DATA GENERATION (Your code) ---
giorni <- 1:100

simula_covid <- function(n_giorni, inizio_boom, durata_salita, durata_plateau, picco_max) {
  trend <- numeric(n_giorni)

  for (i in 1:n_giorni) {
    if (i < inizio_boom) {
      trend[i] <- 1
    } else if (i < (inizio_boom + durata_salita)) {
      progress <- (i - inizio_boom) / durata_salita
      trend[i] <- 1 + (picco_max * progress)
    } else if (i < (inizio_boom + durata_salita + durata_plateau)) {
      trend[i] <- picco_max
    } else {
      passi_dal_plateau <- i - (inizio_boom + durata_salita + durata_plateau)
      trend[i] <- max(0, picco_max - (picco_max / 30 * passi_dal_plateau))
    }
  }

  casi_reali <- rpois(n_giorni, lambda = trend)
  return(casi_reali)
}

# Generate data
y1 <- simula_covid(n_giorni = 100, inizio_boom = 5, durata_salita = 15,
                   durata_plateau = 20, picco_max = 800)
y2 <- simula_covid(n_giorni = 100, inizio_boom = 25, durata_salita = 25,
                   durata_plateau = 15, picco_max = 400)

# --- OPTION 1: LOG TRANSFORM (Recommended for count data) ---
cat("=== Option 1: Log-transformed data ===\n\n")

# Add 1 to avoid log(0)
y_log <- list(
  Regione1 = log(y1 + 1),
  Regione2 = log(y2 + 1)
)

results_log <- gibbs_sampler_dpy(
  y = y_log,
  K = 20,
  n_iter = 3000,
  burn_in = 1500,
  alpha_vec = c(0.25, 0.25),
  theta_vec = c(5, 5),
  alpha0 = 0.25,
  theta0 = 5,
  sigma2 = 1,
  estimate_sigma = TRUE,
  verbose = TRUE
)

cat(sprintf("\nLog-transformed results:\n"))
cat(sprintf("Mean clusters (Regione 1): %.1f\n",
            mean(sapply(results_log$n_clusters, function(x) x[1]))))
cat(sprintf("Mean clusters (Regione 2): %.1f\n",
            mean(sapply(results_log$n_clusters, function(x) x[2]))))

# Plot
plot_both_processes(y_log, results_log)


# --- OPTION 2: FIRST DIFFERENCES (For time series patterns) ---
cat("\n\n=== Option 2: First differences ===\n\n")

# Compute daily changes
diff_y1 <- diff(y1)
diff_y2 <- diff(y2)

y_diff <- list(
  Regione1 = diff_y1,
  Regione2 = diff_y2
)

results_diff <- gibbs_sampler_dpy(
  y = y_diff,
  K = 20,
  n_iter = 3000,
  burn_in = 1500,
  alpha_vec = c(0.25, 0.25),
  theta_vec = c(5, 5),
  alpha0 = 0.25,
  theta0 = 5,
  sigma2 = 50,  # Larger variance for differences
  estimate_sigma = TRUE,
  verbose = TRUE
)

cat(sprintf("\nFirst differences results:\n"))
cat(sprintf("Mean clusters (Regione 1): %.1f\n",
            mean(sapply(results_diff$n_clusters, function(x) x[1]))))
cat(sprintf("Mean clusters (Regione 2): %.1f\n",
            mean(sapply(results_diff$n_clusters, function(x) x[2]))))


# --- OPTION 3: STANDARDIZED DATA ---
cat("\n\n=== Option 3: Standardized data ===\n\n")

# Standardize to mean=0, sd=1
y_std <- list(
  Regione1 = scale(y1)[,1],
  Regione2 = scale(y2)[,1]
)

results_std <- gibbs_sampler_dpy(
  y = y_std,
  K = 20,
  n_iter = 3000,
  burn_in = 1500,
  alpha_vec = c(0.25, 0.25),
  theta_vec = c(5, 5),
  alpha0 = 0.25,
  theta0 = 5,
  sigma2 = 1,
  estimate_sigma = TRUE,
  verbose = TRUE
)

cat(sprintf("\nStandardized results:\n"))
cat(sprintf("Mean clusters (Regione 1): %.1f\n",
            mean(sapply(results_std$n_clusters, function(x) x[1]))))
cat(sprintf("Mean clusters (Regione 2): %.1f\n",
            mean(sapply(results_std$n_clusters, function(x) x[2]))))


# --- COMPREHENSIVE COMPARISON ---
cat("\n\n=== COMPARISON OF METHODS ===\n")
cat("Log transform:  ",
    sprintf("R1=%.1f clusters, R2=%.1f clusters\n",
            mean(sapply(results_log$n_clusters, function(x) x[1])),
            mean(sapply(results_log$n_clusters, function(x) x[2]))))
cat("First diff:     ",
    sprintf("R1=%.1f clusters, R2=%.1f clusters\n",
            mean(sapply(results_diff$n_clusters, function(x) x[1])),
            mean(sapply(results_diff$n_clusters, function(x) x[2]))))
cat("Standardized:   ",
    sprintf("R1=%.1f clusters, R2=%.1f clusters\n",
            mean(sapply(results_std$n_clusters, function(x) x[1])),
            mean(sapply(results_std$n_clusters, function(x) x[2]))))


# --- VISUALIZATION ---
par(mfrow = c(2, 2))

# Original data
plot(giorni, y1, type = "l", col = "red", lwd = 2,
     main = "Original COVID Data", xlab = "Days", ylab = "Cases")
lines(giorni, y2, col = "blue", lwd = 2)
legend("topright", legend = c("Regione 1", "Regione 2"),
       col = c("red", "blue"), lty = 1, lwd = 2)

# Log-transformed
plot(giorni, y_log$Regione1, type = "l", col = "red", lwd = 2,
     main = "Log-transformed", xlab = "Days", ylab = "log(Cases + 1)")
lines(giorni, y_log$Regione2, col = "blue", lwd = 2)

# First differences
plot(giorni[-1], diff_y1, type = "l", col = "red", lwd = 2,
     main = "First Differences (Daily Change)", xlab = "Days", ylab = "Change")
lines(giorni[-1], diff_y2, col = "blue", lwd = 2)
abline(h = 0, lty = 2, col = "gray")

# Standardized
plot(giorni, y_std$Regione1, type = "l", col = "red", lwd = 2,
     main = "Standardized", xlab = "Days", ylab = "Z-score")
lines(giorni, y_std$Regione2, col = "blue", lwd = 2)
abline(h = 0, lty = 2, col = "gray")

par(mfrow = c(1, 1))


cat("\n\n=== RECOMMENDATION ===\n")
cat("For COVID case count data, LOG TRANSFORMATION is usually best because:\n")
cat("1. It handles the large scale differences (0 to 800)\n")
cat("2. It stabilizes variance across the range\n")
cat("3. It makes the data more Gaussian-like\n")
cat("4. It preserves the relative patterns of growth/decline\n")
cat("\nUse FIRST DIFFERENCES if you want to identify phases:\n")
cat("- Negative differences = declining cases\n")
cat("- Positive differences = increasing cases\n")
cat("- Near-zero differences = plateau\n")
