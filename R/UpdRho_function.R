#' Updating the correlation coefficient for latent factors.
#' 
#' @param curr_r Current value of the correlation coefficient.
#' @param curr_S Current value of the correlation matrix for the latent
#' factors.
#' @param corr_C Phylogenetic correlation matrix.
#' @param latfac Current values of the latent factors. Matrix with rows
#' corresponding to the species and columns to factors.
#' @param mh_n Parameter n in the Beta proposal for updating the 
#' correlation coefficient for latent factors.
#' @param prior_rho Hyperparameters of the beta prior for the weight rho in 
#' the correlation matrix for latent factors. Defaults to (5, 5).
#' 
UpdRho <- function(curr_r, curr_S, corr_C, latfac, mh_n, prior_rho) {
  
  r <- list(new_value_r = curr_r, new_value_S = curr_S, accepted = 0)
  num_obs <- nrow(latfac)
  
  # Generating the proposal:
  r_prop <- rbeta(1, shape1 = mh_n * curr_r, shape2 = mh_n * (1 - curr_r))
  S_prop <- r_prop * corr_C + (1 - r_prop) * diag(num_obs)
  
  # Calculating the likelihood, prior and proposal for proposed value:
  log_prior_prop <- dbeta(r_prop, prior_rho[1], prior_rho[2], log = TRUE)
  log_lik_prop <- sum(mvnfast::dmvn(t(latfac), mu = rep(0, num_obs),
                                    sigma = S_prop, log = TRUE))
  log_prop <- dbeta(r_prop, shape1 =  mh_n * curr_r,
                    shape2 = mh_n * (1 - curr_r), log = TRUE)
  
  # Calculating the likelihood, prior and proposal for current value:
  log_prior_curr <- dbeta(curr_r, prior_rho[1], prior_rho[2], log = TRUE)
  log_lik_curr <- sum(mvnfast::dmvn(t(latfac), mu = rep(0, num_obs),
                                    sigma = curr_S, log = TRUE))
  log_curr <- dbeta(curr_r, shape1 = mh_n * r_prop,
                    shape2 = mh_n * (1 - r_prop), log = TRUE)
  
  AP <- log_lik_prop + log_prior_prop + log_curr
  AP <- AP - (log_lik_curr + log_prior_curr + log_prop)
  
  if (log(runif(1)) < AP) {
    r <- list(new_value_r = r_prop, new_value_S = S_prop, accepted = 1)
  }
  
  return(r)
  
}