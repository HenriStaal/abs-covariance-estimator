Adaptive Beta Shrinkage (ABS) estimator

Implements steps 2–4 of the ABS algorithm from Staal & Flint (manuscript under review):
* **Step 1 (volatility adjustment)** is *not* included here, to stay consistent with the application in Staal & Flint (manuscript under review). If you would like to implement this step, apply it to the returns matrix before calling this function.

**Args:**
* **R**  
  t × n matrix of returns (t time points, n assets). May be raw or already volatility-adjusted.  
* **lambda**  
  (optional) Ridge penalty for step 3. Default: `max(log(n / t), 0)`.  
* **V**  
  (optional) Length-t vector of observation weights. Converted internally to `W = diag(V)`. Default: `rep(1, t)`.  
* **B_prior**  
  Shrinkage prior for β:  
  * `"0"` – zero prior (default, as used in Staal & Flint (manuscript under review))  
  * `"MLE"` – uses the empirical-Bayes MLE of β as prior  
  * *numeric* – fixed prior value specified by the user  

**Output:**
* **P_hat**  
  n × n shrunk correlation matrix.  

**Notation** (in code vs. in Staal & Flint (2025)):
- **P_sample** – ρₛ (sample correlation matrix)  
- **S** – S (standardised returns matrix)  
- **B_hat** – β̂ (ridge‐estimated pairwise beta)  
- **B_prior** – βₚᵣᵢₒᵣ (prior mean for beta)  
- **sigma2_hat_B** – σ̂²ᵦ (sampling variance of β̂)  
- **sigma2_B_hat** – σ²ᵦ̂ (cross‐sectional variance of β̂)  
- **sigma2_prior** – σ²ₚᵣᵢₒᵣ (prior variance used in Vasicek shrinkage)  
- **B_post** – βₚₒₛₜ (posterior mean, i.e., shrunken beta)  
- **P_hat** – P̂ (final shrunken correlation matrix)  
*Note:* An “off” suffix indicates that the variable is a vector that stores the respective values for all off-diagonal entries (i.e. for all i, j such that i ≠ j).

**Requires:**
* `library(Matrix)` for `nearPD()` function
