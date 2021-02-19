#' Trait matching for species interactions
#' 
#' Using the posterior samples for covariate-informed latent factor network
#' model to infer trait imprtance in forming interactions.
#' 
#' @param B Numeric. Number of times to perform the resampling. Default is 500.
#' @param h_choice1 List of functions to be used for choosing h. These
#' functions should include functions that lead to randomness in the choice of
#' the component. If left NULL two choices will be considered: equal weight,
#' and proportional to the absolute value.
#' @param h_choice2 List of functions to be used for choosing h. These
#' functions should include functions that lead to a single component. If left
#' NULL, the maximum absolute element will be chosen.
#' @param lambdas Matrix including the coefficients of the network model. The 
#' rows correspond to posterior samples and the columns to intercept and latent
#' factors components.
#' 
#' @param betas 
#' @param gammas
#' @param Us
#' @param Vs
#' @param Xs
#' @param Ws
#' 
TraitMatching <- function(B = 500, h_choice1 = NULL, h_choice2 = NULL,
                          lambdas, betas, gammas, Us, Vs, Xs, Ws) {
  
  H <- dim(lambdas)[2] - 1
  Nsims <- nrow(lambdas)
  
  # Number of continuous and binary covariates:
  entries_X <- apply(Xs, 2, function(x) length(unique(x[!is.na(x)])))
  pB <- c(sum(entries_X > 2), sum(entries_X == 2))
  entries_W <- apply(Ws, 2, function(x) length(unique(x[!is.na(x)])))
  pP <- c(sum(entries_W > 2), sum(entries_W == 2))
  if (any(c(entries_X[1 : pB[1]] == 2, entries_W[1 : pP[1]] == 2))) {
    stop('Reorder covariates')
  }
  
  # Specifying how the latent factor component will be chosen if left NULL.
  if (is.null(h_choice1)) {
    h_choice1 <- list(f1 = function(x) rep(1, length(x)),
                      f2 = function(x) abs(x))
  }
  if (is.null(h_choice2)) {
    h_choice2 <- list(f1 = function(x) abs(x) * (abs(x) == max(abs(x))))
  }
  
  # Calculating the rescaled lambda by multiplying the value of the coefficient
  # with the standard deviation of the predictor:
  rescaled_lambda <- lambdas[, - 1]
  for (ss in 1 : Nsims) {
    for (hh in 1 : H) {
      x <- outer(Us[ss, , hh], Vs[ss, , hh])
      rescaled_lambda[ss, hh] <- rescaled_lambda[ss, hh] * sd(as.numeric(x))
    }
  }
  
  
  # Number of resampling of components:
  n_h_choice <- length(h_choice1) + length(h_choice2)
  
  # Where to save mean squared residuals:
  resampling_X <- array(NA, dim = c(B, Nsims, sum(pB), n_h_choice))
  resampling_W <- array(NA, dim = c(B, Nsims, sum(pP), n_h_choice))
  
  for (bb in 1 : B) {
    
    if (bb %% 50 == 0) print(bb)
    
    for (ss in 1 : Nsims) {
      
      # Which latent factor will be used:
      wh_h <- NULL
      for (ll in 1 : length(h_choice1)) {
        wh_h <- c(wh_h, sample(1 : H, 1, prob = h_choice1[[ll]](rescaled_lambda[ss, ])))
      }
      for (ll in 1 : length(h_choice2)) {
        wh_h <- c(wh_h, sample(1 : H, 1, prob = h_choice2[[ll]](rescaled_lambda[ss, ])))
      }
      
      for (mm in 1 : sum(pB)) {
        wh_obs <- which(!is.na(Xs[, mm]))
        this_response <- Xs[wh_obs, mm]
        
        for (hh in 1 : n_h_choice) {
          if (!(bb > 1 & hh %in% seq(length(h_choice1) + 1, n_h_choice))) {
            
            # Which latent factor component is used:
            this_component <- wh_h[hh]
            # The intercept
            this_inter <- betas[ss, mm, 1]
            # What is the coefficient of the latent factor at this iteration:
            this_coef <- betas[ss, mm, this_component + 1]
            # Calculating coefficient times latent factor as the linear predictor:
            this_lin_pred <- this_inter + this_coef * Us[ss, wh_obs, this_component]
            if (mm <= pB[1]) { # If continous, sum of the squared residuals:
              resampling_X[bb, ss, mm, hh] <- sum((this_response - this_lin_pred) ^ 2)
            } else { # If binary, classification accuracy (AUROC)
              resampling_X[bb, ss, mm, hh] <- as.numeric(pROC::roc(response = this_response,
                                                                   predictor = this_lin_pred,
                                                                   quiet = TRUE)$auc)
            }
          }
        }
      }
      
      for (ll in 1 : sum(pP)) {
        wh_obs <- which(!is.na(Ws[, ll]))
        this_response <- Ws[wh_obs, ll]
        
        for (hh in 1 : n_h_choice) {
          if (!(bb > 1 & hh %in% seq(length(h_choice1) + 1, n_h_choice))) {
            
            this_component <- wh_h[hh]
            this_inter <- gammas[ss, ll, 1]
            this_coef <- gammas[ss, ll, this_component + 1]
            this_lin_pred <- this_inter + this_coef * Vs[ss, wh_obs, this_component]
            if (ll <= pP[1]) {
              resampling_W[bb, ss, ll, hh] <- sum((this_response - this_lin_pred) ^ 2)
            } else {
              resampling_W[bb, ss, ll, hh] <- as.numeric(pROC::roc(response = this_response,
                                                                   predictor = this_lin_pred,
                                                                   quiet = TRUE)$auc)
            }
          }
        }
      }
    }
  }
  
  return(list(SSR_X = resampling_X, SSR_W = resampling_W))
  
}