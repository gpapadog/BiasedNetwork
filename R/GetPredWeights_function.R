#' Importance sampling weights for a set of units
#' 
#' Corresponds to the likelihood of the species covariates given the generated
#' latent factors and current parameters.
#' 
#' @param pred_covs The covariates of the species we predict on. Can include
#' NAs for missing values.
#' @param pred_latfac The generated latent factors for this species.
#' @param mcmc_coefs Posterior samples for the coefficients in the covariate
#' models.
#' @param mcmc_vars Posterior samples for the residual variances in the models
#' for the continuous covariates.
#' 
GetPredWeights <- function(pred_covs, pred_latfac, mcmc_coefs, mcmc_vars) {
  
  # Continuous and binary covariates:
  num_covs <- c(dim(mcmc_vars)[2], dim(mcmc_coefs)[2] - dim(mcmc_vars)[2])
  Nsims <- dim(mcmc_coefs)[1]
  
  # Density of covariates at current iteration:
  pred_den <- matrix(NA, Nsims, sum(num_covs))
  
  for (jj in 1 : sum(num_covs)) {
    if (!is.na(pred_covs[jj])) {
      # The linear predictor of the model - vector of length Nsims.
      lin_pred <- apply(mcmc_coefs[, jj, ] * cbind(1, pred_latfac), 1, sum)
      if (jj <= num_covs[1]) {
        pred_den[, jj] <- dnorm(pred_covs[jj], lin_pred, sqrt(mcmc_vars[, jj]))
      } else {
        pred_den[, jj] <- dbinom(pred_covs[jj], 1, prob = expit(lin_pred))
      }
    }
  }
  
  r <- apply(pred_den, 1, prod, na.rm = TRUE)
  return(r)

}

