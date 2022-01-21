#' Updating the species probability of being observed in a study
#' 
#' A species probabiliy of being observed, pi, is updated based on the observed
#' interactions, latent traits, and current true interactions.
#' 
#' @param probobs_curr Vector of current values for the probability of
#' observation for the set of species we are updating.
#' @param probobs_others Vector of current values for the probability of
#' observation for the other set of species, the ones that are not being
#' updated.
#' @param curr_inter Matrix with 0,1 entries for the current posterior samples
#' of possible interactions. The rows correspond to the species we are updating
#' and the columns correspond to the other set of species.
#' @param obs_inter Matrix with 0,1 entries for the detected interactions. The
#' rows correspond to the species we are updating and the columns correspond to
#' the other set of species.
#' @param mh_n Numeric. Parameter n in the Beta proposal distribution.
#' @param obs_in_poss Matrix with rows corresponding to the set of species we
#' are updating and rows corresponding to the other set. The matrix includes
#' a count for the number of studies that could have recorded the interaction
#' and did record it.
#' @param unobs_in_poss Matrix with rows corresponding to the set of species we
#' are updating and rows corresponding to the other set. The matrix includes
#' a count for the number of studies that could have recorded the interaction
#' but they did not record it.
#' @param latfac Current values of the latent factors for the set of species
#' whose probability of detection is to be updated. Matrix with rows
#' corresponding to the species and columns corresponding to the factors.
#' @param coefs_probobs Coefficients of the latent factors for the normal model
#' on the logit of probability of observing a given species. If these
#' coefficients are all NA, then this implies that no bias correction for
#' geographical and taxonomical bias is performed, so this part of the model
#' does not inform the updates of the latent factors. Should be of length equal
#' to the number of latent factors used plus 1 for the intercept.
#' @param var_probobs Numeric. Variance in the model for the probability of
#' detecting an interaction. Used only if performing bias correction.
#' 
#' 
UpdProbObs <- function(probobs_curr, probobs_others, curr_inter, obs_inter,
                       mh_n, obs_in_poss, unobs_in_poss, latfac, coefs_probobs,
                       var_probobs) {
  
  # Acquiring quantities that will be used later:
  num_obs <- nrow(latfac)
  prior_mean <- cbind(1, latfac) %*% matrix(coefs_probobs, ncol = 1)
  
  # Generating the proposed values:
  p_prop <- rbeta(num_obs, mh_n * probobs_curr, mh_n * (1 - probobs_curr))
  
  # Likelihood part:
  pipj_prop <- outer(p_prop, probobs_others)
  pipj_curr <- outer(probobs_curr, probobs_others)
  
  lik_p1 <- (pipj_prop / pipj_curr) ^ (curr_inter * obs_in_poss)
  lik_p2 <- ((1 - pipj_prop) / (1 - pipj_curr)) ^ (curr_inter * unobs_in_poss)

  # Prior part:
  prior_lik <- (dnorm(logit(p_prop), prior_mean, sqrt(var_probobs)) /
                  dnorm(logit(probobs_curr), prior_mean, sqrt(var_probobs)))
  
  # Proposal part:
  prop_lik <- (dbeta(probobs_curr, mh_n * p_prop, mh_n * (1 - p_prop)) /
                 dbeta(p_prop, mh_n * probobs_curr, mh_n * (1 - probobs_curr)))
  
  # Acceptance probability.
  AP <- apply(lik_p1 * lik_p2, 1, prod) * prior_lik * prop_lik
  
  u <- runif(num_obs)
  
  acc <- (u <= AP)
  r <- probobs_curr
  r[acc] <- p_prop[acc]
  return(list(new_values = r, accepted = acc, model_values = prior_mean))
  
}