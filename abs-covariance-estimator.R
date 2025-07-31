library(Matrix)

ABS <- function(R, lambda = NULL, V = NULL, B_prior = "0") {
  
  # Finding the dimensions n and t
  n <- ncol(R)
  t <- nrow(R)
  
  # Setting lambda
  if (is.null(lambda)) {
    lambda <- max(log(n / t), 0)
  }
  
  # Setting V
  if (is.null(V)) {
    V <- rep(1, t)
  } else if (length(V) != t) {
    stop("the weight vector V must be of length t")
  }
  
  # Sample correlation matrix
  P_sample <- cor(R)
  
  #STEP (2) of algorithm: standardisation
  
    # Standardize returns
    S <- scale(R, center = TRUE, scale = TRUE)
  
  
  #STEP (3) of algorithm: weighted ridge regression
  
    # Create the weight matrix, W
    W <- diag(V)
    
    # Let M = s'Ws. This is needed in order to vectorise the ridge regression computations.
    M <- t(S) %*% W %*% S
    
    # Off‐diagonal indices
    offdiag_idx <- which(row(M) != col(M))
    i_idx       <- row(M)[offdiag_idx]
    j_idx       <- col(M)[offdiag_idx]
    
    # Elements for ridge regression
    M_jj <- M[cbind(j_idx, j_idx)]
    M_ji <- M[cbind(j_idx, i_idx)]
    M_ii <- M[cbind(i_idx, i_idx)]
    
    # Ridge‐estimate of beta for off-diagonal entries
    B_hat_off <- (M_ji) / (M_jj + lambda)
  
  
  # STEP (4) of algorithm: beta shrinkage via the Vasicek (1973) adjustment
  
    # Create a vector to populate with the B_prior value for the off-diagonal entries, so that vector operations can be performed in the Vasicek adjustment
    B_prior_off <- numeric(length(offdiag_idx))
    
    if (B_prior == "0") {
      B_prior_off[] <- 0
    } else if (B_prior == "MLE") {
      B_prior_off[] <- mean(B_hat_off, na.rm = TRUE)
    }
    
    # Calculate the residual variances and the sampling variance, i.e. sigma2_hat_B, for each off-diagonal asset pair
    sum_of_squared_residuals <- M_ii - 2 * B_hat_off * M_ji + (B_hat_off^2) * M_jj
    regression_residual_variance <- sum_of_squared_residuals / (t - 1)
    sigma2_hat_B_off <- regression_residual_variance / (M_jj + lambda)
    
    # 2) Compute sigma2_B_hat, i.e. the variance of all the beta estimates
    sigma2_B_hat <- var(B_hat_off, na.rm = TRUE)
    
    # 3) Calculate pairwise prior variance
    sigma2_prior_off <- pmax(sigma2_B_hat - sigma2_hat_B_off, 1e-6)
    
    # Perform the Vasicek adjustment
    B_post_off <- ((B_prior_off / sigma2_prior_off) + (B_hat_off / sigma2_hat_B_off)) / 
      ((1 / sigma2_prior_off) + (1 / sigma2_hat_B_off))
    
    # Populate the final, shrunk estimate of the correlation matrix
    P_hat <- matrix(1, n, n)
    P_hat[cbind(i_idx, j_idx)] <- B_post_off
  
  # Average the two pairwise betas (though this is only necessary if non-unit weights are used in W)
  P_hat <- (P_hat + t(P_hat)) / 2
  
  # Test whether the shrunk correlation matrix is positive definite. If not, find the nearest positive definite matrix using the nearPD function. 
  if(any(eigen(P_hat, only.values = TRUE)$values <= 0)) {
    P_hat <- as.matrix(nearPD(P_hat)$mat)
  }
  
  return(P_hat)
}
