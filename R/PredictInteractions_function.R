#' Predicting the probability of species interaction
#' 
#' Using samples from the posterior distribution of model parameters to predict
#' the probability of interaction for in-sample co-occuring, in-sample
#' neglected, half-in-sample, and out-of-sample pairs of species.
#' 
#' @param pred_X Matrix of covariates for the first set of species. Number of
#' rows are the number of measured covariates, and number of columns is the
#' number of species. The matrix can include missing data.
#' @param pred_W Matrix of covariates for the second set of species. Number of
#' rows are the number of measured covariates, and number of columns is the
#' number of species. The matrix can include missing data.
#' @param mcmc List of posterior samples as returned from the MCMC function.
#' @param corr_mat List of two elements. Each element corresponds to each set
#' of species, and it includes the species phylogenetic correlation matrix
#' extended to include the species on which we predict. For each species, we
#' only use the sub-matrix including the observed data and the target species.
#' If the species is included in the original data, it still has to be included
#' here, but its entries can be arbitrary, since they are not used.
#' 
#' @return List including four arrays of dimensions corresponding to the
#' number of the two sets of species and posterior samples. The first array
#' includes draws from the probability of an interaction between the two
#' species based on the model, these probabilities are not in any way adjusted
#' by importance sampling weights. The second array includes draws from the
#' unweighted distribution without importance sampling reweighting. The third
#' array includes the importance sampling weights. The fourth array includes
#' samples from the posterior distribution of interactions after re-weighting
#' according to the importance sampling weights.
#' 
#' @export
#' 
PredictInteractions <- function(pred_X, pred_W, mcmc, corr_mat) {
  
  # ---------- STEP 0: Get quantities and transform data ----------- #
  
  Nsims <- dim(mcmc$Ls)[1]
  nB <- dim(mcmc$Ls)[2]
  nP <- dim(mcmc$Ls)[3]
  use_H <- dim(mcmc$Us)[3]
  
  # If there is only one species, turn vector to matrix:
  if (class(pred_X) == 'numeric') {
    pred_X <- matrix(pred_X, nrow = 1)
  }
  if (class(pred_W) == 'numeric') {
    pred_W <- matrix(pred_W, nrow = 1)
  }
  
  # Number of species to predict.
  total_nB <- nrow(pred_X)
  total_nP <- nrow(pred_W)
  
  
  # ---------- STEP 1: Where to save results ----------- #
  
  # Where to save results:
  unweigh_prob <- array(NA, dim = c(Nsims, total_nB, total_nP))
  unweigh_pred <- array(NA, dim = c(Nsims, total_nB, total_nP))
  imp_sam_weights <- array(NA, dim = c(Nsims, total_nB, total_nP))
  weigh_prob <- array(NA, dim = c(total_nB, total_nP))
  weigh_pred <- array(NA, dim = c(Nsims, total_nB, total_nP))
  
  
  # ---------- STEP 2: Latent factors and weights ----------- #
  
  # For each set of species, we acquire their latent factors and importance
  # sampling weights in advance, to avoid re-calculating them multiple times.
  #
  # For the species not in the original sample, we have to generate
  # latent factors and calculate importance sampling weights.
  # For units in the original sample, the latent factors are acquired from
  # the MCMC and their weights are 1.
  
  cat('Predicting latent factors and calculating importance sampling weights',
      'for the first set of species.', fill = TRUE)
  
  weights_X <- matrix(1, nrow = Nsims, ncol = total_nB)
  pred_U <- array(NA, dim = c(Nsims, total_nB, use_H))
  
  for (ii in 1 : total_nB) {
    
    if (ii <= nB) {  # In sample.
      pred_U[, ii, ] <- mcmc$Us[, ii, ]
    } else {  # Out-of-sample.
      use_C <- corr_mat[[1]][c(1 : nB, ii), c(1 : nB, ii)]
      pred_U[, ii, ] <- GetPredLatFac(mcmc_latfac = mcmc$Us,
                                      mcmc_rho = mcmc$rU, C = use_C)
      weights_X[, ii] <- GetPredWeights(pred_covs = pred_X[ii, ],
                                        pred_latfac = pred_U[, ii, ],
                                        mcmc_coefs = mcmc$betas,
                                        mcmc_vars = mcmc$sigmasq_m)
    }
  }
  
  
  cat('Predicting latent factors and calculating importance sampling weights',
      'for the second set of species.', fill = TRUE)
  
  weights_W <- matrix(1, nrow = Nsims, ncol = total_nP)
  pred_V <- array(NA, dim = c(Nsims, total_nP, use_H))
  
  for (jj in 1 : total_nP) {
    
    if (jj <= nP) {
      pred_V[, jj, ] <- mcmc$Vs[, jj, ]
    } else {
      use_C <- corr_mat[[2]][c(1 : nP, jj), c(1 : nP, jj)]
      pred_V[, jj, ] <- GetPredLatFac(mcmc_latfac = mcmc$Vs,
                                      mcmc_rho = mcmc$rV, C = use_C)
      weights_W[, jj] <- GetPredWeights(pred_covs = pred_W[jj, ],
                                        pred_latfac = pred_V[, jj, ],
                                        mcmc_coefs = mcmc$gammas,
                                        mcmc_vars = mcmc$sigmasq_l)
    }
  }
  
  
  # ---------- STEP 3: Performing the predictions ----------- #
  
  # What type of pair is it? In-sample, out-of-sample, half-in-sample?
  # Prediction will proceed separately for each type of pair.
  # Prediction for in-sample pairs is based on the samples from the MCMC.
  
  cat('Predicting the interactions.', fill = TRUE)
  
  for (ii in 1 : total_nB) {
    for (jj in 1 : total_nP) {
      
      # Assigning the importance sampling weights as the product of the weights
      # from the two species.
      imp_sam_weights[, ii, jj] <- weights_X[, ii] * weights_W[, jj]
      
      # If this is an in-sample pair, mcmc gives the posterior of interaction:
      if (ii <= nB & jj <= nP) {
        unweigh_pred[, ii, jj] <- mcmc$Ls[, ii, jj]
        weigh_pred[, ii, jj] <- unweigh_pred[, ii, jj]
      }
      
      # If the pair is half-in-sample or out of sample, we have to generate L.
      if (ii > nB | jj > nP) {
        
        # Calculating probability of interaction according to the model:
        lin_pred <- mcmc$lambdas[, - 1] * pred_U[, ii, ] * pred_V[, jj, ]
        lin_pred <- mcmc$lambdas[, 1] + apply(lin_pred, 1, sum)
        unweigh_prob[, ii, jj] <- expit(lin_pred)
        unweigh_pred[, ii, jj] <- rbinom(Nsims, 1, unweigh_prob[, ii, jj])
        
        # Drawing from the probability of interaction according to the
        # importance sampling weights:
        wh_0 <- which(unweigh_pred[, ii, jj] == 0)
        prob_0 <- sum(imp_sam_weights[wh_0, ii, jj])
        prob_1 <- sum(imp_sam_weights[- wh_0, ii, jj])
        prob_1 <- prob_1 / (prob_0 + prob_1)
        weigh_prob[ii, jj] <- prob_1
        weigh_pred[, ii, jj] <- rbinom(Nsims, 1, prob = prob_1)
        
      }
    }
  }
  
  return(list(unweigh_prob = unweigh_prob, unweigh_pred = unweigh_pred,
              imp_sam_weights = imp_sam_weights, weigh_prob = weigh_prob,
              weigh_pred = weigh_pred, weights_X = weights_X,
              weights_W = weights_W))
  
}








