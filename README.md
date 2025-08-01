The R script abs-covariance-estimator.R provides the ABS() function, which implements the Adaptive Beta Shrinkage estimator for covariance matrices (Staal & Flint, manuscript under review).

The function implements steps 2–4, as detailed in Staal & Flint (manuscript under review). **Step 1 (volatility adjustment)** is *not* included here, to stay consistent with the application in Staal & Flint (manuscript under review). To include Step 1, apply a volatility adjustment to the returns matrix before calling this function.

**Usage:**
```r
# Install dependency
install.packages("Matrix")

# Load the ABS function
source("abs-covariance-estimator.R")

# Example call
# R is a t × n returns matrix
P_hat <- ABS(R, lambda  = 0.1, V = rep(1, nrow(R)), B_prior = "0")
```

## Quick Start

### 1. Source directly from GitHub
```r
install.packages("Matrix")  
source("https://raw.githubusercontent.com/HenriStaal/abs-covariance-estimator/main/abs-covariance-estimator.R")  
P_hat <- ABS(R)

### 2. Clone locally
```bash
git clone https://github.com/HenriStaal/abs-covariance-estimator.git  
cd abs-covariance-estimator

install.packages("Matrix")  
source("abs-covariance-estimator.R")  
P_hat <- ABS(R)



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

**Notation** (in code vs. that used in Staal & Flint (manuscript under review)):

- **P_sample** – ρ<sub>s</sub> (sample correlation matrix)  
- **S** – S (standardised returns matrix)  
- **B_hat** – β<sup>̂</sup> (ridge-estimated pairwise beta)  
- **B_prior** – β<sub>prior</sub> (prior mean for beta)  
- **sigma2_hat_B** – σ<sup>2</sup><sub>β̂</sub> (sampling variance of β̂)  
- **sigma2_B_hat** – σ<sup>2</sup><sub>β̂</sub> (cross-sectional variance of β̂)  
- **sigma2_prior** – σ<sup>2</sup><sub>prior</sub> (prior variance used in Vasicek shrinkage)  
- **B_post** – β<sub>post</sub> (posterior mean, i.e., shrunken beta)  
- **P_hat** – P<sup>̂</sup> (final shrunken correlation matrix) 

*Note:* An “off” suffix indicates that the variable is a vector that stores the respective values for all off-diagonal entries (i.e. for all i, j such that i ≠ j).

**Requires:**
* `library(Matrix)` for `nearPD()` function

**Citation**
Staal, H. (2025). ABS Covariance Estimator [R script]. GitHub.
https://github.com/HenriStaal/abs-covariance-estimator
