#' Trait matching for species interactions
#' 
#' Fitting multivariate regression of the logit of fitted interaction
#' probability values on all covariates and using permutation to get their
#' distribution under the null of no association.
#' 
#' @param B Numeric. Number of times to perform the resampling. Default is 500.
#' @param mod_pL1s Posterior samples for the fitted probabilities of the
#' interaction model. Returned by the function MCMC.
#' @param Xs Posterior samples of imputed values for the covariates of the
#' first set of species. Returned by the function MCMC.
#' @param Ws Posterior samples of imputed values for the covariates of the
#' second set of species. Returned by the function MCMC.
#' @param obs_X The design matrix of covariates for the first set of species.
#' Number of rows is number of species and number of columns is number of
#' covariates.
#' @param obs_W The design matrix of covariates for the second set of species.
#' Number of rows is number of species and number of columns is number of
#' covariates.
#' 
#' @export
#' 
TraitMatching3 <- function(B = 500, mod_pL1s, Xs, Ws, obs_X, obs_W) {
  
  Nsims <- dim(mod_pL1s)[1]
  nB <- dim(mod_pL1s)[2]
  nP <- dim(mod_pL1s)[3]
  sum_pB <- ncol(obs_X)
  sum_pP <- ncol(obs_W)
  
  logit_mod_pL1s <- logit(mod_pL1s)
  
  # Where to save R squared:
  coefs_resampling_X <- array(NA, dim = c(B, Nsims, sum_pB))
  coefs_resampling_W <- array(NA, dim = c(B, Nsims, sum_pP))
  
  coefs_obs_X <- array(NA, dim = c(Nsims, sum_pB))
  coefs_obs_W <- array(NA, dim = c(Nsims, sum_pP))
  
  
  cat('Covariates of first set of species. ')
  
  # For each posterior sample
  for (ss in 1 : Nsims) {
    
    if (ss %% 100 == 0) print(ss)
    
    this_response <- logit_mod_pL1s[ss, , ]
    these_pred <- obs_X
    for (mm in 1 : sum_pB) {
      wh_na <- which(is.na(these_pred[, mm]))
      if (length(wh_na) > 0) {
        these_pred[wh_na, mm] <- Xs[[mm]][ss, ]
      }
    }
    
    # Calculate the abs_coefs for each plant species, and take average
    abs_coefs <- matrix(NA, nP, sum_pB)
    for (jj in 1 : nP) {
      abs_coefs[jj, ] <- abs(lm(this_response[, jj] ~ these_pred)$coef[- 1])
    }
    coefs_obs_X[ss, ] <- apply(abs_coefs, 2, mean)
    
    # Reorder the covariates and do the same.  
    for (bb in 1 : B) {
      these_pred <- these_pred[sample(1 : nB, nB, replace = FALSE), ]
      abs_coefs <- matrix(NA, nP, sum_pB)
      for (jj in 1 : nP) {
        abs_coefs[jj, ] <- abs(lm(this_response[, jj] ~ these_pred)$coef[- 1])
      }
      coefs_resampling_X[bb, ss, ] <- apply(abs_coefs, 2, mean)
    }
  }
  
  
  cat('Covariates of second set of species. ')
  
  
  # For each posterior sample
  for (ss in 1 : Nsims) {
    
    if (ss %% 100 == 0) print(ss)
    
    this_response <- logit_mod_pL1s[ss, , ]
    these_pred <- obs_W
    for (ll in 1 : sum_pP) {
      wh_na <- which(is.na(these_pred[, ll]))
      if (length(wh_na) > 0) {
        these_pred[wh_na, ll] <- Ws[[ll]][ss, ]
      }
    }
    
    # Calculate the abs_coefs for each plant species, and take average
    abs_coefs <- matrix(NA, nB, sum_pP)
    for (ii in 1 : nB) {
      abs_coefs[ii, ] <- abs(lm(this_response[ii, ] ~ these_pred)$coef[- 1])
    }
    coefs_obs_W[ss, ] <- apply(abs_coefs, 2, mean)
    
    # Reorder the covariates and do the same.  
    for (bb in 1 : B) {
      these_pred <- these_pred[sample(1 : nP, nP, replace = FALSE), ]
      abs_coefs <- matrix(NA, nB, sum_pP)
      for (ii in 1 : nB) {
        abs_coefs[ii, ] <- abs(lm(this_response[ii, ] ~ these_pred)$coef[- 1])
      }
      coefs_resampling_W[bb, ss, ] <- apply(abs_coefs, 2, mean)
    }
  }
  
  coefs_resampling_X <- apply(coefs_resampling_X, c(1, 3), mean)
  coefs_resampling_W <- apply(coefs_resampling_W, c(1, 3), mean)
  coefs_obs_X <- apply(coefs_obs_X, 2, mean)
  coefs_obs_W <- apply(coefs_obs_W, 2, mean)
  
  dimnames(coefs_resampling_X) <- list(boot = 1 : B, cov = 1 : sum_pB)
  names(dimnames(coefs_resampling_X)) <- c('boot', 'covariate')
  
  dimnames(coefs_resampling_W) <- list(boot = 1 : B, cov = 1 : sum_pP)
  names(dimnames(coefs_resampling_W)) <- c('boot', 'covariate')
  
  names(coefs_obs_X) <- 1 : sum_pB
  names(coefs_obs_W) <- 1 : sum_pP
  
  return(list(coefs_resampling_X = coefs_resampling_X,
              coefs_resampling_W = coefs_resampling_W,
              coefs_obs_X = coefs_obs_X,
              coefs_obs_W = coefs_obs_W))
  
}