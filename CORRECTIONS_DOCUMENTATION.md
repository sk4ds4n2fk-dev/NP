# Correzioni al Codice DPY Mixture Model

## Data: 2026-01-12

Questo documento descrive le correzioni apportate all'implementazione del modello Dependent Pitman-Yor Process Mixture per allinearlo esattamente alle equazioni del paper "Beta-product dependent Pitman-Yor processes for Bayesian inference" di Bassetti, Casarin, and Leisen (2014).

---

## Errori Identificati e Correzioni

### 1. **Errore Critico nella Funzione `update_psi`**

#### Localizzazione
Funzione `get_log_lik_psi_V` all'interno di `update_psi` (circa linea 157-177 nel codice originale)

#### Problema
L'implementazione calcolava la log-likelihood di V dato ψ usando esponenti **sbagliati** per il termine $(1-V)$:

**Codice Errato (Originale):**
```r
# Per V0
log_prob_v0 <- -alp * log(v_obj$V0) +
  a1 * log(1 - v_obj$V0) -              # ERRORE: dovrebbe essere (a1 - 1)
  lbeta(1 - alp, a1)

# Per Vik
log_prob_v1 <- (a1 - alp) * log(v_obj$Vik[1,]) +
  (a2 + alp * k_seq) * log(1 - v_obj$Vik[1,]) -  # ERRORE: dovrebbe essere (a2 + alp*k - 1)
  log_norm_vi
```

**Codice Corretto:**
```r
# Per V0 - Equazione (41) del paper
# P{V_{0k} ∈ dv|ψ} = v^{-l}(1-v)^{α₁-1} / B(1-l, α₁)
log_prob_v0 <- -alp * log(v_obj$V0) +
  (a1 - 1) * log(1 - v_obj$V0) -        # CORRETTO
  lbeta(1 - alp, a1)

# Per Vik - Equazione (41) del paper
# P{V_{ik} ∈ dv|k, ψ} = v^{α₁-l}(1-v)^{α₂+lk-1} / B(α₁+1-l, α₂+lk)
log_prob_v1 <- (a1 - alp) * log(v_obj$Vik[1,]) +
  (a2 + alp * k_seq - 1) * log(1 - v_obj$Vik[1,]) -  # CORRETTO
  log_norm_vi
```

#### Riferimento al Paper
- **Equazione (41)** nella terza immagine fornita dal paper
- Questa equazione specifica le densità prior per le variabili stick-breaking

#### Impatto
Questo errore causava:
- Sampling scorretto dei parametri di concentrazione ψ = (α₁, α₂, l)
- Distribuzione posterior incorretta per i parametri del processo
- Potenziali problemi di convergenza dell'MCMC

---

### 2. **Verifica Funzione `log_Qk`**

#### Stato
✅ **CORRETTA** - Nessuna modifica necessaria

La funzione implementa correttamente l'equazione per $Q_k(v_k|D, \tilde{\psi})$ come definita dopo l'equazione (25):

$$Q_k(v_k|D, \tilde{\psi}) = v_{0k}^{-l+A_{1k}+A_{2k}}(1-v_{0k})^{\alpha_1-1} \prod_{i=1,2} v_{ik}^{A_{ik}+\alpha_1-l}(1-v_{ik})^{\alpha_2+lk-1}(1-v_0v_i)^{B_{ik}}$$

Dove:
- $A_{ik} = |\{t \in \{1,...,T_i\} : D_{i,t} = k\}|$ (numero di osservazioni nel cluster k)
- $B_{ik} = |\{t \in \{1,...,T_i\} : D_{i,t} > k\}|$ (numero di osservazioni in cluster superiori a k)

---

### 3. **Verifica Altre Funzioni**

#### `update_sticks_MH` (Equazione 27)
✅ **CORRETTA** - Implementa correttamente:
$$P\{V^* \in dv^*|\tilde{\psi}, D\} \propto \prod_{k \in \mathcal{D}^*} Q_k(v_k|D, \tilde{\psi})dv_k$$

#### `update_atoms` (Equazioni 44-49, Appendice B.1)
✅ **CORRETTA** - Implementa correttamente le full conditionals per:
- μ_ik (parametri component-specific)
- σ²_ik (varianze component-specific)
- μ_0k (parametri comuni)
- σ²_0k (varianze comuni)

#### Update di D (Equazione 31)
✅ **CORRETTA** - Implementa correttamente:
$$P\{D_{i,t} = d|\tilde{\vartheta}, V, U, \tilde{\psi}, Y\} \propto \mathcal{K}_t(Y_{i,t}|\tilde{\vartheta}_{id}, Z_t) \mathbb{I}\{U_{i,t} \leq W_{i,d}\}$$

#### Update di U (Equazione 30)
✅ **CORRETTA** - Implementa correttamente:
$$P\{U_{i,t} \in du|V, Y, \tilde{\vartheta}, D\} = \frac{\mathbb{I}\{u \leq W_{i,D_{i,t}}\}}{W_{i,D_{i,t}}}du$$

---

## Struttura del Block Gibbs Sampler

Il sampler implementa correttamente la struttura a blocchi descritta nel paper (Section 5):

### Iterazione MCMC (ogni iterazione campiona in ordine):

1. **U | V, D** (Equazione 30)
   - Variabili ausiliarie per slice sampling

2. **Espansione adattiva di K**
   - Segue Walker (2007) e Kalli et al. (2011)
   - Aggiunge componenti finché $\sum_{k=1}^K W_{ik} > \min_t U_{i,t}$

3. **D | Y, V, U, ψ, ϑ** (Equazione 31)
   - Allocazioni cluster via Gibbs sampling

4. **V | D, ψ** (Equazione 27)
   - Stick-breaking variables via Metropolis-Hastings
   - Usa trasformazione logit per migliorare mixing

5. **ϑ | D, Y** (Equazione 24, Appendice B.1)
   - Atomi della mistura via Gibbs sampling
   - Normale-Inverse Gamma coniugate

6. **ψ | V, D** (Equazioni 28, 50-51)
   - Parametri di concentrazione via Metropolis-Hastings
   - Prior: Gamma(ζ₁₁, ζ₂₁) per α₁, Gamma(ζ₁₂, ζ₂₂) per α₂, Uniform per l

---

## Miglioramenti Aggiuntivi

### 1. Documentazione Dettagliata
- Aggiunto riferimenti espliciti alle equazioni del paper
- Commenti che spiegano ogni passo dell'algoritmo
- Struttura chiara con sezioni ben definite

### 2. Commenti sulle Equazioni
Ogni funzione ora include:
- Numero dell'equazione di riferimento dal paper
- Formula matematica in forma LaTeX nei commenti
- Spiegazione dei parametri

### 3. Nomi delle Variabili
Allineati alla notazione del paper:
- ψ = (α₁, α₂, l) invece di generici nomi
- ϑ per gli atomi
- D per le allocazioni
- V per stick-breaking variables

---

## Testing e Validazione

### Test Consigliati

1. **Test di Convergenza**
   ```r
   # Verificare che i trace plots di ψ mostrino convergenza
   res <- run_DPY_sampler(data, n_iter=10000, burnin=5000)
   plot_diagnostics(res)
   ```

2. **Test di Acceptance Rate**
   ```r
   # Verificare che acceptance rate per V sia tra 0.23 e 0.44
   mean(res$acc_rate[5000:10000])
   ```

3. **Test di Prior Predictive**
   ```r
   # Generare dati dal prior e verificare che il sampler recuperi i parametri
   ```

---

## Confronto Versione Originale vs Corretta

| Componente | Versione Originale | Versione Corretta | Riferimento |
|------------|-------------------|-------------------|-------------|
| V0 prior exponent | α₁ | α₁ - 1 | Eq. (41) |
| Vik prior exponent | α₂ + lk | α₂ + lk - 1 | Eq. (41) |
| log_Qk | ✓ Corretta | ✓ Corretta | Eq. (25) |
| update_atoms | ✓ Corretta | ✓ Corretta | Eq. (44-49) |
| Slice sampling | ✓ Corretta | ✓ Corretta | Sec. 5.2 |

---

## Note Importanti

### Interpretazione della Densità Beta

Per una distribuzione Beta(a, b), la densità è:
$$f(v) = \frac{v^{a-1}(1-v)^{b-1}}{B(a,b)}$$

Quindi se dall'equazione (41) abbiamo:
$$P\{V_{0k} \in dv|\psi\} = \frac{v^{-l}(1-v)^{\alpha_1-1}}{B(1-l, \alpha_1)}$$

I parametri Beta sono:
- Primo parametro: $a = 1 - l$ (perché $a - 1 = -l$)
- Secondo parametro: $b = \alpha_1$ (perché $b - 1 = \alpha_1 - 1$)

Questo è **esattamente** ciò che è implementato nella versione corretta.

---

## Conclusioni

Le correzioni apportate assicurano che l'implementazione sia **matematicamente equivalente** alle equazioni del paper. Gli errori corretti erano critici per la correttezza statistica del metodo, in particolare per l'inferenza sui parametri di concentrazione del processo Pitman-Yor.

La versione corretta ora implementa fedelmente:
- Block Gibbs sampler (Section 5)
- Full conditionals esatte (Equations 24, 27, 28, 30, 31)
- Prior corretti per stick-breaking variables (Equation 41)
- Full conditionals per modelli Gaussiani (Appendix B.1)

---

## File Prodotti

1. `DPY_mixture_model_corrected.R` - Implementazione corretta completa
2. `CORRECTIONS_DOCUMENTATION.md` - Questo documento

---

## Riferimenti

Bassetti, F., Casarin, R., & Leisen, F. (2014). Beta-product dependent Pitman-Yor processes for Bayesian inference. Journal of Econometrics, 180(1), 49-72.

Walker, S. G. (2007). Sampling the Dirichlet mixture model with slices. Communications in Statistics—Simulation and Computation, 36(1), 45-54.

Kalli, M., Griffin, J. E., & Walker, S. G. (2011). Slice sampling mixture models. Statistics and Computing, 21(1), 93-105.
