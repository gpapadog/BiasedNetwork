#' Updating the coefficients in the traits model.
#' 
#' @param obs_cov Matrix of observed values for the covariates. Rows
#' correspond to units and columns to covariates.
#' @param latfac Latent factor positions for all subjects. Matrix of rows
#' equal to the observations and columns equal to the number of factors.
#' @param resid_var Value of the residual variance for continuous traits.
#' @param shr_var Value of the variances across the H factor components.
#' @param extra_var Values of the additional variances of all coefficients.
#' Matrix of dimensions equal to the number of covariates and number of
#' factors.
#' @param curr_coefs Matrix of current coefficients. Dimensions correspond to
#' covariate and intercept/factor. They are used to draw Polya-Gamma latent
#' variables for the update of coefficients of binary models.
#' @param num_covs Vector including two values: number of continuous and binary
#' covariates.
#' @param prior_mu0 Mean of the prior distribution on intercepts.
#' @param prior_sigmasq0 Variance of the prior distribution on intercepts.
#' 
UpdTraitCoef <- function(obs_cov, latfac, resid_var, shr_var, extra_var,
                         curr_coefs, num_covs, prior_mu0, prior_sigmasq0)  {
  
  des_mat <- cbind(1, latfac)
  des_mat_sq <- t(des_mat) %*% des_mat
  num_obs <- nrow(obs_cov)
  Hval <- dim(latfac)[2]
  
  r <- matrix(NA, nrow = sum(num_covs), ncol = Hval + 1)
  om <- matrix(NA, nrow = nrow(obs_cov), ncol = num_covs[2])
  
  for (jj in 1 : sum(num_covs)) {
    
    prior_m <- matrix(c(prior_mu0, rep(0, Hval)), ncol = 1)
    prior_S_inv <- diag(1 / c(prior_sigmasq0, shr_var * extra_var[jj, ]))
    
    if (jj <= num_covs[1]) {  # Continuous covariates:
      
      new_S <- solve(des_mat_sq / resid_var[jj] + prior_S_inv)
      new_m <- t(des_mat) %*% matrix(obs_cov[, jj], ncol = 1) / resid_var[jj]

    } else {  # Binary covariates:
      
      # Drawing the Polya-Gamma variables.
      omegas <- BayesLogit::rpg(num_obs, 1, des_mat %*% matrix(curr_coefs[jj, ], ncol = 1))
      om[, jj - num_covs[1]] <- omegas
      
      new_S <- solve(t(des_mat) %*% diag(omegas) %*% des_mat + prior_S_inv)
      new_m <- t(des_mat) %*% matrix(obs_cov[, jj] - 1 / 2, ncol = 1)

    }
    
    new_m <- new_S %*% (new_m + prior_S_inv %*% prior_m)
    r[jj, ] <- mvnfast::rmvn(1, mu = new_m, sigma = new_S)
    
  }
  
  return(list(coefs = r, omegas = om))
  
}

