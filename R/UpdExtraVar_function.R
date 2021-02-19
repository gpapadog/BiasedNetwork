#' Updating the additional variance component for coefficients.
#' 
#' Function that updates the extra variance, and can be used for coefficients
#' across multiple models.
#' 
#' @param mod_coef Values for coefficients whose extra variance will be
#' updated. These can be the coefficients of the physical traits (matrix format
#' with rows corresponding to the trait, and columns to the latent factors),
#' coefficients of the network model (vector with length equal to the number of
#' latent factors), coefficients of the models for probability of observing a
#' certain species (vector).
#' @param shr_var Value of the variances across the H factor components.
#' @param prior_spec Values of parameters in the prior specification.
#' 
UpdExtraVar <- function(mod_coef, shr_var, prior_spec) {
  
  # New shape parameter:
  new_a <- (prior_spec[1] + 1) / 2
  
  # New rate parameter:
  if (class(mod_coef) == 'numeric') {
    new_b_add <- mod_coef ^ 2 / shr_var
    
  } else {  # Class is matrix:
    new_b_add <- as.numeric(sweep(mod_coef ^ 2, 2, FUN = '/', shr_var))
  }
  
  new_b <- (prior_spec[2] + new_b_add) / 2
  
  # Drawing values:
  r <- 1 / rgamma(n = length(mod_coef), shape = new_a, rate = new_b)
  
  if (class(mod_coef) == 'matrix') {
    r <- matrix(r, nrow = nrow(mod_coef), ncol = ncol(mod_coef))
  }
  
  return(r)
}