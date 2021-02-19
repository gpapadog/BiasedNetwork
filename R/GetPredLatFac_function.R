#' Samples from the distribution of the new latent factor.
#' 
#' @param mcmc_latfac Posterior samples from the distribution of latent factors
#' for species in the observed data.
#' @param mcmc_rho Posterior distribution of the weight in the latent factor
#' correlation matrix.
#' @param C Correlation matrix on species including the observed data species
#' and the one we wish to predict for.
#' 
#' @return Matrix of with number of rows equal to the number of posterior
#' samples and number of columns equal to the number of latent factors.
#' 
GetPredLatFac <- function(mcmc_latfac, mcmc_rho, C) {
  
  # Simply drawing from a conditional of a multivariate normal distribution.
  
  Nsims <- dim(mcmc_latfac)[1]
  use_H <- dim(mcmc_latfac)[3]
  n <- dim(C)[1] - 1
  
  # Where to save results:
  r <- matrix(NA, nrow = Nsims, ncol = use_H)
  
  for (ss in 1 : Nsims) {
    
    # Current correlation matrix of latent factors:
    this_C <- mcmc_rho[ss] * C + (1 - mcmc_rho[ss]) * diag(n + 1)
    # Converting these to the conditional variance, same for all h.
    this_S11 <- this_C[n + 1, n + 1]  # This is always equal to 1.
    this_S12 <- matrix(this_C[n + 1, 1 : n], nrow = 1)
    this_S22 <- this_C[1 : n, 1 : n]
    this_S22_inv <- chol2inv(chol(this_S22))
    this_ssq <- this_S11 - this_S12 %*% this_S22_inv %*% t(this_S12)
    
    for (hh in 1 : use_H) {
      this_latfac <- matrix(mcmc_latfac[ss, , hh], ncol = 1)
      this_mean <- this_S12 %*% this_S22_inv %*% this_latfac
      r[ss, hh] <- rnorm(1, mean = this_mean, sd = sqrt(this_ssq))
    }
  }
  
  return(r)
  
}
