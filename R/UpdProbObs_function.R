#' Updating the species probability of being observed in a study
#' 
#' A species probabiliy of being observed, pi, is updated based on the observed
#' interactions, latent traits, and current true interactions.
#' 
#' 
UpdProbObs <- function(probobs_curr, probobs_others, curr_inter, obs_inter,
                       mh_n, n_studies, latfac, coefs_probobs, var_probobs) {
  
  # Acquiring quantities that will be used later:
  num_obs <- nrow(latfac)
  prior_mean <- cbind(1, latfac) %*% matrix(coefs_probobs, ncol = 1)
  
  # Generating the proposed values:
  p_prop <- rbeta(num_obs, mh_n * probobs_curr, mh_n * (1 - probobs_curr))
  
  # Likelihood part:
  pipj_prop <- outer(p_prop, probobs_others)
  pipj_curr <- outer(probobs_curr, probobs_others)
  
  lik_part1 <- ((1 - (1 - pipj_prop) ^ n_studies) /
                  (1 - (1 - pipj_curr) ^ n_studies))
  lik_part2 <- (1 - pipj_prop) / (1 - pipj_curr)
  lik <- (lik_part1 ^ (obs_inter * curr_inter) *
            lik_part2 ^ (n_studies * (1 - obs_inter) * curr_inter))
  
  # Prior part:
  prior_lik <- (dnorm(logit(p_prop), prior_mean, sqrt(var_probobs)) /
                  dnorm(logit(probobs_curr), prior_mean, sqrt(var_probobs)))
  
  # Proposal part:
  prop_lik <- (dbeta(probobs_curr, mh_n * p_prop, mh_n * (1 - p_prop)) /
                 dbeta(p_prop, mh_n * probobs_curr, mh_n * (1 - probobs_curr)))
  
  # Acceptance probability.
  AP <- apply(lik, 1, prod) * prior_lik * prop_lik
  
  u <- runif(num_obs)
  
  acc <- (u <= AP)
  r <- probobs_curr
  r[acc] <- p_prop[acc]
  return(list(new_values = r, accepted = acc, model_values = prior_mean))
  
}