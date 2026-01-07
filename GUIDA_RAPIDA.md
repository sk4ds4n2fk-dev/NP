# 🚀 GUIDA RAPIDA - Come Eseguire il Codice

## ✅ PASSO 1: Verifica R

Apri il terminale/prompt dei comandi e verifica che R sia installato:

```bash
R --version
```

Se non è installato:
- **Windows/Mac**: [Scarica R da qui](https://cran.r-project.org/)
- **Linux**: `sudo apt-get install r-base`

## ✅ PASSO 2: Scarica i File

### Opzione A: Se hai Git configurato
```bash
git clone <url-repository>
cd NP
git checkout claude/implement-paper-code-r-agiCf
```

### Opzione B: Copia manuale
Scarica questi file nella stessa cartella:
- `beta_dpy_simple.R` ⭐ (QUESTO È IL FILE PRINCIPALE)
- Oppure tutti i file se vuoi la versione completa

## ✅ PASSO 3: Avvia R

### Da Terminale:
```bash
cd /percorso/alla/cartella/NP
R
```

### Da RStudio:
1. Apri RStudio
2. File → Open File → Seleziona `beta_dpy_simple.R`
3. Session → Set Working Directory → To Source File Location

## ✅ PASSO 4: Esegui il Codice

### TEST RAPIDO (30 secondi):
```r
# Carica il file
source("beta_dpy_simple.R")

# Esegui test veloce
test_quick()
```

**Output atteso:**
```
Test 1: Stick-breaking... OK
Test 2: DPY generation... OK
Test 3: MCMC rapido... OK

Tutti i test passati!
```

### ESEMPIO COMPLETO (2-3 minuti):
```r
# Esegui esempio con grafici
risultati <- example_simple()
```

**Output atteso:**
```
=== ESEMPIO SEMPLICE ===

Generazione dati...
Processo 1: n=100, media=2.52, sd=2.45
Processo 2: n=100, media=2.27, sd=2.35

Esecuzione MCMC...
Iterazione 100/500
Iterazione 200/500
...

=== RISULTATI ===
Campioni posteriori: 250
Media sigma2: 1.234
Processo 1: 3 cluster attivi
Processo 2: 3 cluster attivi
```

## ✅ PASSO 5: Usa con i TUOI Dati

```r
# Carica il codice
source("beta_dpy_simple.R")

# I tuoi dati (lista di vettori numerici)
miei_dati <- list(
  serie1 = c(1.2, 3.4, 5.6, 2.1, ...),
  serie2 = c(2.3, 4.5, 3.2, 1.8, ...)
)

# Esegui l'analisi
risultati <- gibbs_sampler_simple(
  y = miei_dati,
  K = 10,           # Max numero di cluster
  n_iter = 1000,    # Iterazioni totali
  burn_in = 500,    # Burn-in
  alpha = 0.5,      # Parametro discount
  theta = 1,        # Parametro strength
  verbose = TRUE
)

# Guarda i risultati
summary(risultati$sigma2)
plot(risultati$sigma2, type = "l")
```

## 🔧 Risoluzione Problemi

### ❌ Errore: "cannot open file"
**Soluzione:**
```r
# Controlla dove sei
getwd()

# Cambia directory
setwd("/percorso/corretto/alla/cartella")

# Verifica che il file esista
file.exists("beta_dpy_simple.R")
```

### ❌ Errore: "could not find function"
**Soluzione:**
```r
# Ricarica il file
source("beta_dpy_simple.R")
```

### ❌ Il codice è troppo lento
**Soluzione:**
```r
# Usa meno iterazioni
risultati <- gibbs_sampler_simple(
  y = miei_dati,
  K = 5,          # Meno cluster
  n_iter = 200,   # Meno iterazioni
  burn_in = 100
)
```

### ❌ "object 'y' not found"
**Soluzione:**
```r
# I dati devono essere una LISTA di vettori
y <- list(c(1,2,3,4), c(5,6,7,8))  # ✓ Corretto
y <- c(1,2,3,4)                     # ✗ Sbagliato
```

## 📊 Interpretazione Risultati

```r
# Dopo aver eseguito gibbs_sampler_simple()

# 1. Varianza posteriore
mean(risultati$sigma2)           # Media posteriore
quantile(risultati$sigma2, c(0.025, 0.975))  # Intervallo 95%

# 2. Numero di cluster
for(j in 1:length(miei_dati)) {
  n_clust <- length(unique(risultati$final_allocations[[j]]))
  cat("Serie", j, ":", n_clust, "cluster\n")
}

# 3. Grafico trace
plot(risultati$sigma2, type="l", main="Convergenza")

# 4. Istogramma posteriore
hist(risultati$sigma2, main="Distribuzione Posteriore")
```

## 📧 Serve Aiuto?

Se hai problemi:
1. Copia l'errore esatto che vedi
2. Copia il comando che hai usato
3. Chiedi aiuto specificando questi dettagli

## ⚡ Comandi Rapidi

```r
# Setup completo in 3 righe
source("beta_dpy_simple.R")
test_quick()                    # Verifica funziona
risultati <- example_simple()   # Esegui esempio

# Con i tuoi dati in 2 righe
source("beta_dpy_simple.R")
risultati <- gibbs_sampler_simple(y = lista_tuoi_dati)
```

---

**🎯 INIZIA DA QUI:** Apri R e scrivi: `source("beta_dpy_simple.R")` poi `test_quick()`
