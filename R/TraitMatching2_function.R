#' Trait matching for species interactions
#' 
#' @param B Numeric. Number of times to perform the resampling. Default is 500.
#' @param obs_only Logical. If set to TRUE only the observations with observed
#' covariate values will be used. Defaults to FALSE.
#' 
TraitMatching2 <- function(B = 500, mod_pL1s, Xs, Ws, obs_X, obs_W,
                           obs_only = FALSE) {
  
  Nsims <- dim(mod_pL1s)[1]
  nB <- dim(mod_pL1s)[2]
  nP <- dim(mod_pL1s)[3]
  sum_pB <- ncol(obs_X)
  sum_pP <- ncol(obs_W)
  
  logit_mod_pL1s <- logit(mod_pL1s)
  
  # Where to save R squared:
  rsq_resampling_X <- array(NA, dim = c(B, Nsims, sum_pB))
  rsq_resampling_W <- array(NA, dim = c(B, Nsims, sum_pP))
  
  rsq_obs_X <- array(NA, dim = c(Nsims, sum_pB))
  rsq_obs_W <- array(NA, dim = c(Nsims, sum_pP))
  
  
  cat('Covariates of first set of species. ')
  
  # For each covariate
  for (mm in 1 : sum_pB) {
    cat(mm, ' ')
    
    # For each posterior sample
    for (ss in 1 : Nsims) {
      
      this_response <- logit_mod_pL1s[ss, , ]
      
      this_cov <- obs_X[, mm]
      wh_na <- which(is.na(this_cov))
      
      if (length(wh_na) > 0) {
        if (obs_only) {
          this_cov <- this_cov[- wh_na]
          this_response <- this_response[- wh_na, ]
        } else {
          this_cov[wh_na] <- Xs[[mm]][ss, ]
        }
      }
      
      # Calculate the rsquared for each plant species, and take average
      rsq <- rep(NA, nP)
      for (jj in 1 : nP) {
        rsq[jj] <- cor(this_response[, jj], this_cov) ^ 2
      }
      rsq_obs_X[ss, mm] <- mean(rsq)
      
      # Reorder the covariates and do the same.  
      for (bb in 1 : B) {
        this_cov <- this_cov[sample(1 : length(this_cov), length(this_cov), replace = FALSE)]
        rsq <- rep(NA, nP)
        for (jj in 1 : nP) {
          rsq[jj] <- cor(this_response[, jj], this_cov) ^ 2
        }
        rsq_resampling_X[bb, ss, mm] <- mean(rsq)
      }
    }
  }
  
  cat('Covariates of second set of species. ')
  
  # For each covariate
  for (ll in 1 : sum_pP) {
    cat(ll, ' ')
    
    # For each posterior sample
    for (ss in 1 : Nsims) {
      
      this_response <- logit_mod_pL1s[ss, , ]
      
      this_cov <- obs_W[, ll]
      wh_na <- which(is.na(this_cov))
      if (length(wh_na) > 0) {
        if (obs_only) {
          this_cov <- this_cov[- wh_na]
          this_response <- this_response[, - wh_na]
        } else {
          this_cov[wh_na] <- Ws[[ll]][ss, ]
        }
      }
      
      # Calculate the rsquared for each plant species, and take average
      rsq <- rep(NA, nB)
      for (ii in 1 : nB) {
        rsq[ii] <- cor(this_response[ii, ], this_cov) ^ 2
      }
      rsq_obs_W[ss, ll] <- mean(rsq)
      
      # Reorder the covariates and do the same.  
      for (bb in 1 : B) {
        this_cov <- this_cov[sample(1 : length(this_cov), length(this_cov), replace = FALSE)]
        rsq <- rep(NA, nB)
        for (jj in 1 : nB) {
          rsq[jj] <- cor(this_response[ii, ], this_cov) ^ 2
        }
        rsq_resampling_W[bb, ss, ll] <- mean(rsq)
      }
    }
  }
  
  rsq_resampling_X <- apply(rsq_resampling_X, c(1, 3), mean)
  rsq_resampling_W <- apply(rsq_resampling_W, c(1, 3), mean)
  rsq_obs_X <- apply(rsq_obs_X, 2, mean)
  rsq_obs_W <- apply(rsq_obs_W, 2, mean)
  
  dimnames(rsq_resampling_X) <- list(boot = 1 : B, cov = 1 : sum_pB)
  names(dimnames(rsq_resampling_X)) <- c('boot', 'covariate')
  
  dimnames(rsq_resampling_W) <- list(boot = 1 : B, cov = 1 : sum_pP)
  names(dimnames(rsq_resampling_W)) <- c('boot', 'covariate')
  
  names(rsq_obs_X) <- 1 : sum_pB
  names(rsq_obs_W) <- 1 : sum_pP
  
  return(list(rsq_resampling_X = rsq_resampling_X,
              rsq_resampling_W = rsq_resampling_W,
              rsq_obs_X = rsq_obs_X,
              rsq_obs_W = rsq_obs_W))
  
}