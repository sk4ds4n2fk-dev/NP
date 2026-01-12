# Fixing Poor Density Fit in DPY Mixture Model

## Problema Osservato

Il plot mostra densità stimate che non catturano bene la struttura dei dati:
- **Group 1**: Densità troppo piatta/uniforme invece di multimodale
- **Group 2**: Cattura i picchi principali ma perde dettagli e variabilità

## Cause Identificate

### 1. **Varianze Troppo Grandi**
Il modello originale usa:
```r
Y_it ~ N(μ_0k + μ_ik, σ²_0k * σ²_ik)
```

La varianza finale è il **PRODOTTO** di due varianze:
- `σ²_0k` = varianza comune
- `σ²_ik` = varianza specifica per gruppo

Con prior vaghi (`s_sq_i = 25, s0_sq = 25`), questo prodotto diventa molto grande, causando:
- Sovrapposizione eccessiva tra componenti
- Densità troppo smooth/piatta
- Perdita di struttura multimodale

### 2. **Hyperparametri Troppo Vaghi**
```r
s_sq_i = c(25, 25)    # Troppo grande
s0_sq = 25            # Troppo grande
lambda = 3            # Prior debole per varianze
eps = 3               # Prior debole per varianze
```

### 3. **Threshold Troppo Alto**
```r
threshold = 0.02
```
Elimina componenti con peso < 2%, che potrebbero essere importanti per catturare le code.

### 4. **Modello Troppo Complesso**
Il modello a due livelli (`μ_0k + μ_ik`) aggiunge flessibilità non necessaria per molti casi.

---

## Soluzioni Proposte

### ✅ SOLUZIONE 1: Modello Semplificato (RACCOMANDATO)

Usa il file `dpy_mixture_improved.R` con `use_simplified_model = TRUE`:

```r
source("dpy_mixture_improved.R")

res <- run_DPY_sampler_improved(
  data = data,
  hyperparams = initialize_hyperparameters_improved(),
  use_simplified_model = TRUE,  # ← IMPORTANTE
  n_iter = 10000,
  burnin = 5000,
  thin = 5
)

best <- dahl(res)
plot_density_improved(best, res, threshold = 0.005)
```

**Cambiamenti:**
- Modello: `Y_it ~ N(μ_k, σ²_k)` invece di `N(μ_0k + μ_ik, σ²_0k * σ²_ik)`
- Varianza singola invece di prodotto
- Più interpretabile e stabile

---

### ✅ SOLUZIONE 2: Hyperparametri Migliorati

```r
initialize_hyperparameters_improved <- function() {
  list(
    s_sq_i = c(5, 5),      # ← Ridotto da 25
    s0_sq = 5,             # ← Ridotto da 25
    lambda = 5,            # ← Aumentato da 3
    eps = 5,               # ← Aumentato da 3

    psi_z11 = 5,           # ← Ridotto da 10
    psi_z21 = 1,           # ← Ridotto da 10
    psi_z12 = 5,
    psi_z22 = 1,

    psi_init = c(0.5, 2, 0.1),  # ← α₁ più piccolo, α₂ più grande

    scale_prop = 0.3,      # ← Ridotto per stabilità
    kappa_sq = 0.2,        # ← Ridotto
    thin = 5
  )
}
```

**Effetti:**
- Prior più informativi sulle varianze → componenti più concentrate
- Più flessibilità sui parametri di concentrazione
- Migliore stabilità numerica

---

### ✅ SOLUZIONE 3: Threshold Più Basso

```r
plot_density_improved(best, res, threshold = 0.005)  # Invece di 0.02
```

Include più componenti nella densità finale.

---

### ✅ SOLUZIONE 4: Più Iterazioni

```r
res <- run_DPY_sampler_improved(
  ...,
  n_iter = 15000,   # ← Aumenta da 5000
  burnin = 10000,   # ← Aumenta da 2000
  thin = 5
)
```

Migliore convergenza, specialmente per modelli complessi.

---

## Come Usare

### Quick Start (Modello Semplificato)

```r
# 1. Source the improved version
source("dpy_mixture_improved.R")

# 2. Prepare your data
data <- list(Y1 = your_data_group1, Y2 = your_data_group2)

# 3. Run sampler
res <- run_DPY_sampler_improved(
  data = data,
  use_simplified_model = TRUE,
  n_iter = 10000,
  burnin = 5000
)

# 4. Get best partition and plot
best <- dahl(res)
plot_density_improved(best, res, threshold = 0.005)

# 5. Check diagnostics
diagnose_results(res)
plot_diagnostics(res)
```

### Full Example

Vedi `example_improved_usage.R` per un esempio completo con:
- Confronto tra modello semplificato, complesso e originale
- Diagnostic plots
- Raccomandazioni

```r
source("example_improved_usage.R")
```

---

## Diagnostics da Controllare

### 1. Acceptance Rate
```
Target: 0.23 - 0.44
Se < 0.23: Aumenta scale_prop
Se > 0.44: Diminuisci scale_prop
```

### 2. Numero di Cluster Attivi
```
Se troppo pochi:
- Riduci psi_init (α₁, α₂)
- Aumenta K iniziale
- Riduci threshold nel plot

Se troppi:
- Aumenta psi_init
- Prior più informativi sulle varianze
```

### 3. Convergenza di ψ
```
Controlla trace plots di (α₁, α₂, l)
Devono mostrare mixing e stazionarietà
```

### 4. Varianze
```
Nel modello semplificato, controlla trace_variances
Varianze troppo grandi → densità piatta
Soluzione: Prior più informativi (lambda, eps)
```

---

## Confronto Modelli

| Feature | Original Complex | Simplified |
|---------|-----------------|------------|
| Parametri | μ_0k + μ_ik, σ²_0k * σ²_ik | μ_k, σ²_k |
| Interpretabilità | Bassa | Alta |
| Varianza | Prodotto (instabile) | Singola (stabile) |
| Convergenza | Lenta | Veloce |
| Use case | Struttura gerarchica complessa | Misture standard |

---

## Troubleshooting

### Problema: Densità ancora troppo piatta

**Soluzioni:**
1. Riduci ulteriormente `s_sq_i` e `s0_sq` (prova 2-3)
2. Aumenta `lambda` ed `eps` (prova 10)
3. Aumenta `n_iter` a 20000+
4. Riduci `threshold` a 0.001

### Problema: Troppi cluster

**Soluzioni:**
1. Aumenta `psi_init[1]` e `psi_init[2]`
2. Prior più restrittivi: `psi_z11 = 15, psi_z12 = 15`
3. Aumenta varianza prior: `s0_sq = 10`

### Problema: Nessuna convergenza

**Soluzioni:**
1. Controlla accettazione rate
2. Adatta `scale_prop` e `kappa_sq`
3. Inizializza meglio i parametri
4. Aumenta burnin

### Problema: Code pesanti non catturate

**Soluzioni:**
1. Usa `threshold = 0.001` invece di 0.005
2. Aumenta `K` iniziale a 30-50
3. Prior meno informativi su varianze (lambda = 3)

---

## File Aggiuntivi

- `dpy_mixture_sampler.R` - Implementazione base corretta
- `dpy_mixture_improved.R` - Versione migliorata con modello semplificato
- `example_improved_usage.R` - Esempi completi
- `CORRECTIONS_DOCUMENTATION.md` - Documentazione correzioni matematiche

---

## Riferimenti

- Modello originale: `dpy_mixture_sampler.R`
- Paper: Bassetti, Casarin, and Leisen (2014) "Beta-product dependent Pitman-Yor processes for Bayesian inference", Journal of Econometrics

---

## Contatti e Contributi

Per problemi o domande:
1. Controlla prima i diagnostics con `diagnose_results(res)`
2. Prova il modello semplificato
3. Aumenta iterazioni prima di modificare altro
