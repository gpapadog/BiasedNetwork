#' Updating the model's latent factors
#' 
#' MCMC sampling of latent factors for either set of species.
#' 
#' @param latfac Current values of the latent factors to be updated. Matrix
#' with number of rows equal to the number of species and number of columns
#' equal to the number of factors.
#' @param latfac_others Current values of the latent factors fot the other set
#' of species, not to be updated. Matrix. Rows correspond to species and
#' columns correspond to latent factors.
#' @param probobs Vector including the values for the probability of observing
#' a given species interaction within a study. Should be of length equal to the
#' number of species among this set.
#' @param coefs_probobs Coefficients of the latent factors for the normal model
#' on the logit of probability of observing a given species. If these
#' coefficients are all NA, then this implies that no bias correction for
#' geographical and taxonomical bias is performed, so this part of the model
#' does not inform the updates of the latent factors. Should be of length equal
#' to the number of latent factors used plus 1 for the intercept.
#' @param var_probobs Numeric. Variance in the model for the probability of
#' detecting an interaction. Used only if performing bias correction.
#' @param obs_covs  Matrix of observed traits. Rows correspond to species and
#' columns correspond to covariates. Continuous covariates should be included
#' first.
#' @param omega_obs_covs The Polya-Gamma omegas from the models of binary
#' covariates. Matrix. Rows correspond to units, columns to binary traits.
#' @param num_covs Vector of length 2 including the number of continuous and
#' binary traits.
#' @param coefs_covs Matrix of latent factor coefficients in the traits models.
#' Rows correspond to traits and columns to factors.
#' @param var_covs Vector of length equal to the number of continuous traits
#' including the residual variance of the corresponding model.
#' @param curr_inter Matrix with rows corresponding to the set of species whose
#' latent factors we are updating and columns corresponding to the other set of
#' species. Entries are 0 or 1.
#' @param coefs_inter Vector including the coefficients in the network model.
#' @param omega_inter The omegas generated from the Polya-Gamma for the model
#' of the true interactions L. The rows should correspond to the units whose
#' latent factors are updated using this function.
#' @param prior_S_inv The inverse of the current correlation matrix for the
#' latent factors that are being updated.
#' 
UpdLatFac <- function(latfac, latfac_others,
                      probobs, coefs_probobs, var_probobs,
                      obs_covs, omega_obs_covs, num_covs, coefs_covs,
                      var_covs,
                      curr_inter, coefs_inter, omega_inter,
                      prior_S_inv) {
  
  
  ret <- latfac
  Hval <- dim(latfac)[2]
  num_obs <- dim(latfac)[1]
  num_obs_other <- dim(latfac_others)[1]
  
  if (any(dim(omega_inter) != c(num_obs, num_obs_other))) {
    stop('Wrong dimensions of omega_inter.')
  }
  
  # If the coefficients for the observation model are NA, it implies that bias
  # correction is not performed.
  bias_cor <- !(all(is.na(coefs_probobs)))

  
  for (hh in 1 : Hval) {
    
    # Attaching the column of 1s to the latent factors:
    des_mat <- cbind(1, ret[, - hh])
    
    
    # Part 1(A): Combining the prior and the likelihood from continuous traits.
    # It leads to an updated normal distribution with parameters mpart, Spart.

  
    # Partial covariance matrix.
    coef_diag <- sum(coefs_covs[1 : num_covs[1], hh + 1] ^ 2 / var_covs)
    Spart_inv <- coef_diag * diag(num_obs) + prior_S_inv

    # Partial mean vector.
    mpart <- rep(0, num_obs)
    for (jj in 1 : num_covs[1]) {
      pred <- des_mat %*% matrix(coefs_covs[jj, - (hh + 1)], ncol = 1)
      resid <- obs_covs[, jj] - pred
      mpart <- mpart + coefs_covs[jj, hh + 1] * resid / var_covs[jj]
    }
    
    
    # Part 1(B): If bias correction is performed, I update the prior and
    # continuous traits to includethe model for the probability of observing.
    
    if (bias_cor) {
      
      coef_diag <- coefs_probobs[hh + 1] ^ 2 / var_probobs
      Spart_inv <- Spart_inv + coef_diag * diag(num_obs)
      
      pred <- des_mat %*% matrix(coefs_probobs[- (hh + 1)], ncol = 1)
      resid <- logit(probobs) - pred
      mpart <- mpart + coefs_probobs[hh + 1] * resid / var_probobs
      
    }
    
    
    # NOTE: mpart and Spart are not truly means and variances because
    # I do not invert Spart_inv and I won't multiple mpart with Spart since
    # I don't need the standardized ones later, and it takes time.
    
    
    
    # Part 2: Updating the variance to include the Polya-Gamma components:
    # model for the true interactions L, and models for binary covariates.
    
    # Starting from what we have (the partial information):
    new_S <- Spart_inv
    
    # Adding the terms for interaction model:
    omega_latfac <- sweep(omega_inter, 2, FUN = '*', latfac_others[, hh] ^ 2)
    add_var <- coefs_inter[hh + 1] ^ 2 * diag(apply(omega_latfac, 1, sum))
    new_S <- new_S + add_var
    
    # Adding the terms for the binary covariate models:
    for (jj in 1 : num_covs[2]) {
      use_coef <- coefs_covs[num_covs[1] + jj, hh + 1]
      new_S <- new_S + use_coef ^ 2 * diag(omega_obs_covs[, jj])
    }
    
    # Inverting:
    new_S <- chol2inv(chol(new_S))
    
    
    
    # Part 3: Updating the mean of the conditional posterior to include the
    # Polya-Gamma components for themodel for the true interactions L, and
    # the models for binary covariates.
    
    
    # Starting from what we have (the partial information)
    new_m <- mpart
    
    # ----- Adding the part corresponding to the L model ------ #
    
    # U*V array (to avoid looping over j)
    ext_UV <- array(rep(ret[, - hh], num_obs_other),
                    dim = c(num_obs, Hval - 1, num_obs_other))
    ext_UV <- sweep(ext_UV, c(2, 3), FUN = '*', t(latfac_others[, - hh]))
    
    # Predictions excluding factor h (using sweep for computational gains):
    pred <- sweep(ext_UV, 2, FUN = '*', coefs_inter[- c(1, hh + 1)])
    pred <- coefs_inter[1] + apply(pred, c(1, 3), sum)
    resid <- (curr_inter - 1 / 2) / omega_inter - pred
    scaled_resid <- sweep(omega_inter * resid, 2, FUN = '*', latfac_others[, hh])
    
    # Updating the mean:
    new_m <- new_m + coefs_inter[hh + 1] * apply(scaled_resid, 1, sum)
    
    # ----- Adding the part corresponding to the binary covariates ------ #
    
    for (jj in 1 : num_covs[2]) {
      use_coef <- coefs_covs[num_covs[1] + jj, hh + 1]
      pred <- des_mat %*% coefs_covs[num_covs[1] + jj, - (hh + 1)]
      resid <- (obs_covs[, num_covs[1] + jj] - 1 / 2) / omega_obs_covs[, jj] - pred
      new_m <- new_m +  use_coef * resid * omega_obs_covs[, jj]
    }
    
    # Multiplying with the variance:
    new_m <- new_S %*% matrix(new_m, ncol = 1)
    
    
    # Part 4: Sampling.
    
    ret[, hh] <- mvnfast::rmvn(1, new_m, new_S)
    
  }
  
  return(ret)
  
}


