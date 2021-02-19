#' MCMC for bipartite network model with trait information and unrecorded
#' interactions
#' 
#' Perform MCMC to acquire samples from the posterior distribution of model
#' parameters for a model for the true interaction matrix among different sets
#' of species with available trait information.
#' 
#' @param obs_A Observed interaction matrix. Contains values of 0 and 1.
#' @param obs_n Matrix with the same dimensions as obs_A. Includes the number
#' of studies in which both species were observed.
#' @param obs_X Matrix of observed covariates for the first set of species.
#' Rows correspond to species, columns to covariates. Continuous covariates
#' should be first, and binary covariates should follow.
#' @param obs_W Same structure as obs_X but for the second set of species.
#' @param Cu Phylogenetic correlation matrix for the first set of species.
#' @param Cv Phylogenetic correlation matrix for the second set of species.
#' @param Nsims Number of posterior samples to be kept.
#' @param burn Number of samples to be burnt in the beginning of the MCMC.
#' @param thin MCMC sampling thinning.
#' @param use_H Number of latent factors to be used. Defaults to 10.
#' @param bias_cor Logical. Whether the model should aim to correct for
#' geographical and taxonomical bias. If set to FALSE, the MCMC is simplied
#' since the observed matrix of interactions is considered the same as the
#' true matrix of interactions. Defaults to TRUE.
#' @param theta_inf Variance value theta infinity in the spike part of the
#' prior distribution for increasing shrinkage. Defaults to 0.01.
#' @param mh_n_pis Parameter n in the Beta proposal for updating pi. Defaults
#' to 100.
#' @param mh_n_pjs Same as mh_n_pis but for the update of pj.
#' @param mh_n_rho Same as mh_n_pis but for the update of rho.
#' @param stick_alpha Alpha value in the increasing shrinkage prior. Defaults
#' to 5.
#' @param prior_theta Hyperparameters of the inverse gamma distribution in the
#' slab part of the increasing shrinkage prior. Defaults to (1, 1).
#' @param pruir_tau Hyperarameters of the inverse gamma prior for the
#' additional variance parameter in the coefficients. Defaults to (5, 5).
#' @param prior_rho Hyperparameters of the beta prior for the weight rho in 
#' the correlation matrix for latent factors. Defaults to (5, 5).
#' @param prior_mu0 Mean of the normal prior distribution for all intercepts.
#' Defaults to 0.
#' @param prior_sigmasq0 Variance of the normal prior distribution for all
#' intercepts. Defaults to 10.
#' @param prior_sigmasq Hyperparameters of the inverse gamma prior on the
#' variance terms of continuous traits. Defaults to (1, 1).
#' @param start_values List that can include starting values for any subset of
#' parameters. Starting values for parameters that are not specified in
#' start_values will be sampled. Defaults to NULL.
#' @param sampling List specifying which parameters should be sampled by
#' setting the value to TRUE, and which should not and should be kept at
#' their starting value by setting the corresponding list element to FALSE.
#' If sampling is set to FALSE for a parameter, it is recommended that the
#' corresponding start_value is specified to be set to the parameter's true
#' value. Defaults to NULL, and when set to NULL all parameters are sampled.
#' 
#' @export
#' 
MCMC <- function(obs_A, obs_n, obs_X, obs_W, Cu, Cv,
                 Nsims, burn, thin, use_H = 10, bias_cor = TRUE,
                 theta_inf = 0.01,
                 mh_n_pis = 100, mh_n_pjs = 100, mh_n_rho = 100,
                 stick_alpha = 5, prior_theta = c(1, 1), prior_tau = c(5, 5),
                 prior_rho = c(5, 5), prior_mu0 = 0, prior_sigmasq0 = 10,
                 prior_sigmasq = c(1, 1), start_values = NULL,
                 sampling = NULL) {
  
  
  # -------------------- PART 0 ------------------- #
  # ------- Specifying what will be sampled. ------- #
  
  if (is.null(sampling)) {
    sampling <- list(L = TRUE, lambda = TRUE, tau = TRUE, beta = TRUE,
                     gamma = TRUE, sigmasq = TRUE, sigmasq_p = TRUE,
                     delta = TRUE, zeta = TRUE, U = TRUE, V = TRUE, v = TRUE,
                     z = TRUE, theta = TRUE, pis = TRUE, pjs = TRUE, rU = TRUE,
                     rV = TRUE, miss_X = TRUE, miss_W = TRUE)
  }
  
  if (!bias_cor) {  # Without bias adjustment, some parameters are not updated:
    cat('Without bias correction a number of parameters will not be sampled.', fill = TRUE)
    sampling$L <- FALSE
    sampling$sigmasq_p <- FALSE
    sampling$delta <- FALSE
    sampling$zeta <- FALSE
    sampling$pis <- FALSE
    sampling$pjs <- FALSE
  }
  
  
  # ---------------- PART 1 ------------- #
  # Getting the parameters that we use throughout:
  
  # Sample size
  nB <- nrow(obs_A)
  nP <- ncol(obs_A)
  
  cat('MCMC on', nB, 'x', nP, 'number of species.', fill = TRUE)
  
  # Which A indices are equal to 0.
  index_A0 <- which(obs_A == 0)
  quant_A0 <- length(index_A0)
  
  # Continuous and binary covariates.
  entries_X <- apply(obs_X, 2, function(x) length(unique(x[!is.na(x)])))
  pB <- c(sum(entries_X > 2), sum(entries_X == 2))
  entries_W <- apply(obs_W, 2, function(x) length(unique(x[!is.na(x)])))
  pP <- c(sum(entries_W > 2), sum(entries_W == 2))
  # Making sure covariates are ordered as continuous first.
  if (any(c(entries_X[1 : pB[1]] == 2, entries_W[1 : pP[1]] == 2))) {
    stop('Reorder covariates')
  }
  
  # If the covariates have missing values, we check for it here.
  any_X_miss <- any(apply(obs_X, 2, function(x) sum(is.na(x))) > 0)
  any_W_miss <- any(apply(obs_W, 2, function(x) sum(is.na(x))) > 0)
  
  # Indices with missing values:
  if (any_X_miss) {
    miss_X_ind <- apply(obs_X, 2, function(x) which(is.na(x)))
  }
  if (any_W_miss) {
    miss_W_ind <- apply(obs_W, 2, function(x) which(is.na(x)))
  }
  
  
  
  # ---------------- PART 2 ------------- #
  # Creating the elements where the MCMC samples will be saved.
  
  Ls <- array(NA, dim = c(Nsims, nB, nP))  # True interaction matrix.
  mod_pL1s <- array(NA, dim = c(Nsims, nB, nP))  # Probability of interaction without correction.
  pL1s <- array(NA, dim = c(Nsims, nB, nP))  # Probability of interaction.
  Us <- array(NA, dim = c(Nsims, nB, use_H))  # Bird latent factors.
  Vs <- array(NA, dim = c(Nsims, nP, use_H))  # Plant latent factors.
  lambdas <- array(NA, dim = c(Nsims, use_H + 1))  # Coefficients, network model.
  taus_lambda <- array(NA, dim = c(Nsims, use_H))  # Extra variance, network model.
  betas <- array(NA, dim = c(Nsims, sum(pB), use_H + 1))  # Coefficients, bird traits.
  gammas <- array(NA, dim = c(Nsims, sum(pP), use_H + 1))  # Coefficients, plant traits.
  taus_beta <- array(NA, dim = c(Nsims, sum(pB), use_H))  # Extra variance, bird traits.
  taus_gamma <- array(NA, dim = c(Nsims, sum(pP), use_H))  # Extra variance, plant traits.
  sigmasq_m <- array(NA, dim = c(Nsims, pB[1]))  # Residual variance, bird traits.
  sigmasq_l <- array(NA, dim = c(Nsims, pP[1]))  # Residual variance, plant traits.
  deltas <- array(NA, dim = c(Nsims, use_H + 1))  # Coefficients, observing bird.
  taus_delta <- array(NA, dim = c(Nsims, use_H))  # Extra variance, observing bird.
  zetas <- array(NA, dim = c(Nsims, use_H + 1))  # Coefficients, observing plant.
  taus_zeta <- array(NA, dim = c(Nsims, use_H))  # Extra variance, observing plant.
  thetas <- array(NA, dim = c(Nsims, use_H))  # Shrinking variance for loadings.
  vs <- array(NA, dim = c(Nsims, use_H))  # Stick breaking prior for shrinking variance.
  omegas <- array(NA, dim = c(Nsims, use_H))  # Length of broken sticks.
  zs <- array(NA, dim = c(Nsims, use_H))  # Which stick break the H components belong to.
  pis <- array(NA, dim = c(Nsims, nB))  # Probabilities of observing bird.
  pi_accepted <- rep(0, nB)  # Number of MCMC iterations for which the bird proposal was accepted.
  pjs <- array(NA, dim = c(Nsims, nP))  # Probabilities of observing plant.
  pj_accepted <- rep(0, nP)  # Number of MCMC iterations for which the plant proposal was accepted.
  sigmasq_pB <- rep(NA, Nsims)  # Residual variance, observing bird.
  sigmasq_pP <- rep(NA, Nsims)  # Residual variance, observing plant.
  rU <- rep(NA, Nsims)  # Correlation, latent factors for birds.
  ru_accepted <- 0
  rV <- rep(NA, Nsims)  # Correlation, latent factors for plants.
  rv_accepted <- 0
  
  if (any_X_miss) {
    Xs <- lapply(miss_X_ind, function(x) matrix(NA, nrow = Nsims, ncol = length(x)))
  }
  if (any_W_miss) {
    Ws <- lapply(miss_W_ind, function(x) matrix(NA, nrow = Nsims, ncol = length(x)))
  }
  
  
  # ---------------- PART 3 ------------- #
  # Setting starting values for the parameters.
  
  this_L <- matrix(rbinom(nB * nP, 1, 1 / 2), nB, nP)
  this_L[obs_A == 1] <- 1
  this_U <- matrix(rnorm(nB * use_H), nB, use_H)
  this_V <- matrix(rnorm(nP * use_H), nP, use_H)
  this_lambda <- rnorm(use_H + 1, 0, 1)
  this_tau_lambda <- 1 / rgamma(use_H, 10, 10)
  this_beta <- matrix(rnorm(sum(pB) * (use_H + 1)), sum(pB), use_H + 1)
  this_gamma <- matrix(rnorm(sum(pP) * (use_H + 1)), sum(pP), use_H + 1)
  this_tau_beta <- matrix(1 / rgamma(sum(pB) * use_H, 10, 10), sum(pB), use_H)
  this_tau_gamma <- matrix(1 / rgamma(sum(pP) * use_H, 10, 10), sum(pP), use_H)
  this_sigmasq_m <- 1 / rgamma(pB[1], 10, 10)
  this_sigmasq_l <- 1 / rgamma(pP[1], 10, 10)
  this_delta <- rnorm(use_H + 1, 0, 1)
  this_tau_delta <- 1 / rgamma(use_H, 10, 10)
  this_zeta <- rnorm(use_H + 1, 0, 1)
  this_tau_zeta <- 1 / rgamma(use_H, 10, 10)
  this_theta <- sort(1 / rgamma(use_H, 5, 5), decreasing = TRUE)
  this_v <- c(rbeta(use_H - 1, 1, 3), 1)
  this_sigmasq_pB <- 1 / rgamma(1, 10, 10)
  this_sigmasq_pP <- 1 / rgamma(1, 10, 10)
  this_pi <- runif(nB, 0, 1)
  this_pj <- runif(nP, 0, 1)
  this_ru <- rbeta(1, 5, 5)
  this_rv <- rbeta(1, 5, 5)
  this_X <- obs_X
  this_W <- obs_W
  this_mod_pL1 <- matrix(NA, nrow = nB, ncol = nP)
  this_pL1 <- matrix(NA, nrow = nB, ncol = nP)
  
  # Draw starting values for missing covariates from observed distribution.
  if (any_X_miss) {
    for (jj in 1 : pB[1]) {
      use_stats <- c(mean(obs_X[, jj], na.rm = TRUE), sd(obs_X[, jj], na.rm = TRUE))
      this_X[miss_X_ind[[jj]], jj] <- rnorm(length(miss_X_ind[[jj]]),
                                            use_stats[1], use_stats[2])
    }
    for (jj in (pB[1] + 1) : sum(pB)) {
      use_stats <- mean(obs_X[, jj], na.rm = TRUE)
      this_X[miss_X_ind[[jj]], jj] <- rbinom(length(miss_X_ind[[jj]]), 1, use_stats)
    }
  }
  
  if (any_W_miss) {
    for (jj in 1 : pP[1]) {
      use_stats <- c(mean(obs_W[, jj], na.rm = TRUE), sd(obs_W[, jj], na.rm = TRUE))
      this_W[miss_W_ind[[jj]], jj] <- rnorm(length(miss_W_ind[[jj]]),
                                            use_stats[1], use_stats[2])
    }
    for (jj in (pP[1] + 1) : sum(pP)) {
      use_stats <- mean(obs_W[, jj], na.rm = TRUE)
      this_W[miss_W_ind[[jj]], jj] <- rbinom(length(miss_W_ind[[jj]]), 1, use_stats)
    }
  }

  # If the starting values have been specified in start_values, use those.
  if (!is.null(start_values)) {
    for (pp in 1 : length(start_values)) {
      assign(x = names(start_values)[pp], start_values[[pp]])
    }
  }
  
  this_omega <- OmegaFromV(v_val = this_v)
  this_Su <- this_ru * Cu + (1 - this_ru) * diag(nB)
  this_Sv <- this_rv * Cv + (1 - this_rv) * diag(nP)
  this_z <- sample(c(1 : use_H), use_H, replace = TRUE, prob = this_omega)
  
  
  # If we do not perform bias correction, change the starting values:
  if (!bias_cor) {
    this_L <- obs_A
    this_delta <- rep(NA, use_H + 1)
    this_tau_delta <- rep(NA, use_H)
    this_zeta <- rep(NA, use_H + 1)
    this_tau_zeta <- rep(NA, use_H)
    this_sigmasq_pB <- NA
    this_sigmasq_pP <- NA
    this_pi <- rep(1, nB)
    this_pj <- rep(1, nP)
    obs_n <- matrix(NA, nB, nP)  # Does not inform anything.
  }
  
  
  # ---------------- PART 4 ------------- #
  # Performing the MCMC.
  
  # Which index we are saving at. Increased by 1 every thin iterations.
  keep_index <- 1
  
  cat('Total number of iterations:', Nsims * thin + burn, fill = TRUE)
  
  for (ss in 1 : (Nsims * thin + burn)) {
    
    if (ss %% 100 == 0) print(ss)
    
    
    # ------ L: Update true interactions. -------- #
    
    if (sampling$L | sampling$lambda | sampling$U | sampling$V) {
      
      # Probability of L = 1 under the network model.
      logit_pijL <- matrix(this_lambda[1], nB, nP)
      for (hh in 1 : use_H) {
        lat_prod <- matrix(this_U[, hh], ncol = 1) %*% matrix(this_V[, hh], nrow = 1)
        logit_pijL <- logit_pijL + this_lambda[hh + 1] * lat_prod
      }
      this_mod_pL1 <- expit(logit_pijL)
      
    }
    
    if (sampling$L) {
      
      # Probability of observing interaction (i, j).
      pipj <- outer(this_pi, this_pj)
      
      # Probability of L = 1, or L = 0, when A = 0.
      pL1 <- this_mod_pL1 * (1 - pipj) ^ obs_n
      pL0 <- 1 - this_mod_pL1
      this_pL1 <- pL1 / (pL0 + pL1)  # Standardize.
      
      # Update L.
      this_L <- matrix(1, nB, nP)
      this_L[index_A0] <- rbinom(n = quant_A0, 1, prob = this_pL1[index_A0])
      
    }
    
    
    
    # ------- lambda: Update the coefficients in network model ------- #
    
    if (sampling$lambda | sampling$U) {
      # Sample the Polya-Gamma random variables:
      omega_L <- BayesLogit::rpg(nB * nP, h = rep(1, nB * nP), z = as.numeric(logit_pijL))
    }
    
    if (sampling$lambda) {
      
      des_mat <- matrix(1, nB * nP)
      for (hh in 1 : use_H) {
        lat_prod <- matrix(this_U[, hh], ncol = 1) %*% matrix(this_V[, hh], nrow = 1)
        des_mat <- cbind(des_mat, as.numeric(lat_prod))
      }
      
      prior_m <- c(prior_mu0, rep(0, use_H))
      prior_S_inv <- diag(1 / c(prior_sigmasq0, this_tau_lambda * this_theta))
      
      # Following line is faster but identical to
      # t(des_mat) %*% diag(omega_L) %*% des_mat
      new_S <- sweep(t(des_mat), 2, FUN = '*', omega_L) %*% des_mat
      new_S <- chol2inv(chol(new_S + prior_S_inv))
      new_m <- t(des_mat) %*% matrix(as.numeric(this_L - 1 / 2), ncol = 1)
      new_m <- new_m + prior_S_inv %*% prior_m
      new_m <- new_S %*% new_m
      
      # Draw new values for lambda:
      this_lambda <- mvnfast::rmvn(1, new_m, new_S)
      
    }
    
    
    # ----------- tau: Update the extra variances ------------ #
    
    if (sampling$tau) {
      
      this_tau_beta <- UpdExtraVar(mod_coef = this_beta[, - 1],
                                   shr_var = this_theta, prior_spec = prior_tau)
      
      this_tau_gamma <- UpdExtraVar(mod_coef = this_gamma[, - 1],
                                    shr_var = this_theta, prior_spec = prior_tau)
      
      this_tau_lambda <- UpdExtraVar(mod_coef = this_lambda[- 1],
                                     shr_var = this_theta, prior_spec = prior_tau)
      
    }
    
    # Without bias correction, the following updates will not be performed:
    
    if (sampling$tau & (sampling$delta | sampling$sigmasq_p | sampling$pis)) {
      this_tau_delta <- UpdExtraVar(mod_coef = this_delta[- 1],
                                    shr_var = this_theta, prior_spec = prior_tau)
    }
    
    if (sampling$tau & (sampling$zeta | sampling$sigmasq_p | sampling$pjs)) {
      this_tau_zeta <- UpdExtraVar(mod_coef = this_zeta[- 1], 
                                   shr_var = this_theta, prior_spec = prior_tau)
    }
    
    
    # ------- beta: Coefficients of bird physical traits --------- #
    
    if (sampling$beta) {
      r <- UpdTraitCoef(obs_cov = this_X, latfac = this_U,
                        resid_var = this_sigmasq_m, shr_var = this_theta,
                        extra_var = this_tau_beta, curr_coefs = this_beta,
                        num_covs = pB, prior_mu0 = prior_mu0,
                        prior_sigmasq0 = prior_sigmasq0)
      this_beta <- r$coefs
      omega_obsX <- r$omegas
    }
    
    
    # -------- gamma: Coefficients of plant physical traits ---------- #
    
    if (sampling$gamma) {
      r <- UpdTraitCoef(obs_cov = this_W, latfac = this_V,
                        resid_var = this_sigmasq_l, shr_var = this_theta,
                        extra_var = this_tau_gamma, curr_coefs = this_gamma,
                        num_covs = pP, prior_mu0 = prior_mu0,
                        prior_sigmasq0 = prior_sigmasq0)
      this_gamma <- r$coefs
      omega_obsW <- r$omegas  
    }
    
    
    # ------- sigmasq: Residual variance for continuous traits -------- #
    
    if (sampling$sigmasq) {
      # For bird continuous traits:  
      new_a <- prior_sigmasq[1] + nB / 2
      for (jj in 1 : pB[1]) {
        resid <- this_X[, jj] - cbind(1, this_U) %*% matrix(this_beta[jj, ], ncol = 1)
        new_b <- prior_sigmasq[2] + sum(resid ^ 2) / 2
        this_sigmasq_m[jj] <- 1 / rgamma(1, shape = new_a, rate = new_b)
      }
      
      # For plant continuous traits:
      new_a <- prior_sigmasq[1] + nP / 2
      for (jj in 1 : pP[1]) {
        resid <- this_W[, jj] - cbind(1, this_V) %*% matrix(this_gamma[jj, ], ncol = 1)
        new_b <- prior_sigmasq[2] + sum(resid ^ 2) / 2
        this_sigmasq_l[jj] <- 1 / rgamma(1, shape = new_a, rate = new_b)
      }  
    }
    
    
    # ----- sigmasq_p: Residual variance for probability of observing ------- #
    
    if (sampling$sigmasq_p) {
      new_a <- prior_sigmasq[1] + nB / 2
      resid <- logit(this_pi) - cbind(1, this_U) %*% matrix(this_delta, ncol = 1)
      new_b <- prior_sigmasq[2] + sum(resid ^ 2) / 2
      this_sigmasq_pB <- 1 / rgamma(1, shape = new_a, rate = new_b)
      
      new_a <- prior_sigmasq[1] + nP / 2
      resid <- logit(this_pj) - cbind(1, this_V) %*% matrix(this_zeta, ncol = 1)
      new_b <- prior_sigmasq[2] + sum(resid ^ 2) / 2
      this_sigmasq_pP <- 1 / rgamma(1, shape = new_a, rate = new_b)  
    }
    
    
    # ------- deltas, zetas: Coefficients of probability of observing ------- #
    
    if (sampling$delta) {
      # deltas.
      des_mat <- cbind(1, this_U)
      prior_m <- matrix(c(prior_mu0, rep(0, use_H)), ncol = 1)
      prior_S_inv <- diag(1 / c(prior_sigmasq0, this_tau_delta * this_theta))
      
      new_S <- chol2inv(chol(t(des_mat) %*% des_mat / this_sigmasq_pB + prior_S_inv))
      new_m <- t(des_mat) %*% matrix(logit(this_pi), ncol = 1) / this_sigmasq_pB
      new_m <- new_S %*% (new_m + prior_S_inv %*% prior_m)
      this_delta <- mvnfast::rmvn(1, new_m, new_S)
    }
    
    if (sampling$zeta) {
      # zetas:
      des_mat <- cbind(1, this_V)
      prior_m <- matrix(c(prior_mu0, rep(0, use_H)), ncol = 1)
      prior_S_inv <- diag(1 / c(prior_sigmasq0, this_tau_zeta * this_theta))
      
      new_S <- chol2inv(chol(t(des_mat) %*% des_mat / this_sigmasq_pP + prior_S_inv))
      new_m <- t(des_mat) %*% matrix(logit(this_pj), ncol = 1) / this_sigmasq_pP
      new_m <- new_S %*% (new_m + prior_S_inv %*% prior_m)
      this_zeta <- mvnfast::rmvn(1, new_m, new_S)
    }
    
    
    # ------- U(h): Latent factors for the birds.
    
    if (sampling$U | sampling$V) {
      omega_L <- matrix(omega_L, nrow = nB, ncol = nP)
    }
    
    if (sampling$U) {
      this_Su_inv <- chol2inv(chol(this_Su))
      this_U <- UpdLatFac(latfac = this_U, latfac_others = this_V, probobs = this_pi,
                          coefs_probobs = this_delta, var_probobs = this_sigmasq_pB,
                          obs_covs = this_X, omega_obs_covs = omega_obsX,
                          num_covs = pB, coefs_covs = this_beta,
                          var_covs = this_sigmasq_m,
                          curr_inter = this_L, coefs_inter = this_lambda,
                          omega_inter = omega_L,
                          prior_S_inv = this_Su_inv)  
    }
    
    
    # ------- V(h): Latent factors for the plants.
    
    if (sampling$V) {
      this_Sv_inv <- chol2inv(chol(this_Sv))
      this_V <- UpdLatFac(latfac = this_V, latfac_others = this_U, probobs = this_pj,
                          coefs_probobs = this_zeta, var_probobs = this_sigmasq_pP,
                          obs_covs = this_W, omega_obs_covs = omega_obsW,
                          num_covs = pP, coefs_covs = this_gamma,
                          var_covs = this_sigmasq_l,
                          curr_inter = t(this_L), coefs_inter = this_lambda,
                          omega_inter = t(omega_L),
                          prior_S_inv = this_Sv_inv)  
    }
    
    
    # ------- v, omega: Stick breaking weights.
    
    if (sampling$v) {
      new_a <- 1 + sapply(1 : use_H, function(x) sum(this_z == x))
      new_b <- stick_alpha + sapply(1 : use_H, function(x) sum(this_z > x))
      this_v <- rbeta(use_H, new_a, new_b)
      this_omega <- OmegaFromV(v_val = this_v)
    }
    
    
    # ------- z: Latent variables for which of H components.
    
    if (sampling$z) {
      
      all_est <- rbind(this_beta, this_gamma, this_lambda)
      all_var <- rbind(this_tau_beta, this_tau_gamma, this_tau_lambda)
      
      if (bias_cor) {
        all_est <- rbind(all_est, this_delta, this_zeta)
        all_var <- rbind(all_var, this_tau_delta, this_tau_zeta)
      }
      
      all_est <- all_est[, - 1]
      all_var <- all_var * prior_theta[2] / prior_theta[1]
      

      ind_norm_lik <- dnorm(all_est, sd = sqrt(all_var * theta_inf), log = TRUE)
      norm_lik <- apply(ind_norm_lik, 2, sum)
      
      tau_lik <- rep(NA, use_H)
      for (hh in 1 : use_H) {
        tau_lik[hh] <- mvtnorm::dmvt(x = all_est[, hh], df = 2 * prior_theta[1],
                                     sigma = diag(all_var[, hh]), log = TRUE)
      }
      
      for (hh in 1 : use_H) {
        
        log_prob <- log(this_omega[1 : hh]) + norm_lik[1 : hh]
        if (hh < use_H) {
          add_on <- log(this_omega[(hh + 1) : use_H]) + tau_lik[(hh + 1) : use_H]
          log_prob <- c(log_prob, add_on)
        }
        
        probs <- exp(log_prob) / sum(exp(log_prob))
        
        this_z[hh] <- sample(1 : use_H, 1, prob = probs)
        
      }
    }
    
    
    # ------- theta: Shrinking variable of parameters at dimension h.
    if (sampling$theta) {
      new_a <- prior_theta[1] + sum(c(pB, pP, 1)) / 2
      add_on <- (apply(this_beta[, - 1] ^ 2 / this_tau_beta, 2, sum) +
                   apply(this_gamma[, - 1] ^ 2 / this_tau_gamma, 2, sum) +
                   this_lambda[- 1] ^ 2 / this_tau_lambda)
      if (bias_cor) {
        new_a <- new_a + 1
        add_on <- (add_on + this_delta[- 1] ^ 2 / this_tau_delta +
                     this_zeta[- 1] ^ 2 / this_tau_zeta)
      }
      new_b <- prior_theta[2] + add_on / 2
      theta_vals <- 1 / rgamma(use_H, shape = new_a, rate = new_b)
      this_theta <- ifelse(this_z <= 1 : use_H, theta_inf, theta_vals)
    }
    
    
    # ------- pis, pjs: Probability of observing for birds and plants.
    
    if (sampling$pis) {
      upd_prob_obs <- UpdProbObs(probobs_curr = this_pi,
                                 probobs_others = this_pj,
                                 curr_inter = this_L,
                                 obs_inter = obs_A,
                                 mh_n = mh_n_pis,
                                 n_studies = obs_n,
                                 latfac = this_U,
                                 coefs_probobs = this_delta,
                                 var_probobs = this_sigmasq_pB)
      this_pi <- upd_prob_obs$new_values
      pi_accepted <- pi_accepted + upd_prob_obs$accepted
    }
    
    if (sampling$pjs) {
      upd_prob_obs <- UpdProbObs(probobs_curr = this_pj,
                                 probobs_others = this_pi,
                                 curr_inter = t(this_L),
                                 obs_inter = t(obs_A),
                                 mh_n = mh_n_pjs,
                                 n_studies = t(obs_n),
                                 latfac = this_V,
                                 coefs_probobs = this_zeta,
                                 var_probobs = this_sigmasq_pP)
      this_pj <- upd_prob_obs$new_values
      pj_accepted <- pj_accepted + upd_prob_obs$accepted
    }
    
    
    # ------ rho: Weight of phylogenetic information in latent factor covariance.
    
    if (sampling$rU) {
      upd_rho <- UpdRho(curr_r = this_ru, curr_S = this_Su, corr_C = Cu,
                        latfac = this_U, mh_n = mh_n_rho, prior_rho = prior_rho)
      this_ru <- upd_rho$new_value_r
      this_Su <- upd_rho$new_value_S
      ru_accepted <- ru_accepted + upd_rho$accepted
    }
    
    if (sampling$rV) {
      upd_rho <- UpdRho(curr_r = this_rv, curr_S = this_Sv, corr_C = Cv,
                        latfac = this_V, mh_n = mh_n_rho, prior_rho = prior_rho)
      this_rv <- upd_rho$new_value_r
      this_Sv <- upd_rho$new_value_S
      rv_accepted <- rv_accepted + upd_rho$accepted
    }
    
    
    # -------- Missing covariate values --------- #
    
    if (any_X_miss & sampling$miss_X) {
      for (jj in 1 : sum(pB)) {
        this_miss <- miss_X_ind[[jj]]
        if (length(this_miss) > 0) {
          X_mean <- cbind(1, this_U[this_miss, , drop = FALSE]) %*% matrix(this_beta[jj, ], ncol = 1)
          if (jj <= pB[1]) {
            this_X[this_miss, jj] <- rnorm(length(this_miss), mean = X_mean,
                                           sd = sqrt(this_sigmasq_m[jj]))
          } else {
            X_prob <- expit(X_mean)
            this_X[this_miss, jj] <- rbinom(length(this_miss), 1, prob = X_prob)
          }
        }
      }
    }
      
    if (any_W_miss & sampling$miss_W) {
      for (jj in 1 : sum(pP)) {
        this_miss <- miss_W_ind[[jj]]
        if (length(this_miss) > 0) {
          W_mean <- cbind(1, this_V[this_miss, ]) %*% matrix(this_gamma[jj, ], ncol = 1)
          if (jj <= pP[1]) {
            this_W[this_miss, jj] <- rnorm(length(this_miss), mean = W_mean,
                                           sd = sqrt(this_sigmasq_l[jj]))
          } else {
            W_prob <- expit(W_mean)
            this_W[this_miss, jj] <- rbinom(length(this_miss), 1, prob = W_prob)
          }
        }
      }
    }
    
    
    # -------------- END OF MCMC UPDATES ------------- #
      
    
    # ------ Saving the results every thin iteration after burn in:
    if ((ss - burn) %% thin == 0 & ss > burn) {
      
      Ls[keep_index, , ] <- this_L
      pL1s[keep_index, , ] <- this_pL1
      mod_pL1s[keep_index, , ] <- this_mod_pL1
      lambdas[keep_index, ] <- this_lambda
      taus_beta[keep_index, , ] <- this_tau_beta
      taus_gamma[keep_index, , ] <- this_tau_gamma
      taus_lambda[keep_index, ] <- this_tau_lambda
      taus_delta[keep_index, ] <- this_tau_delta
      taus_zeta[keep_index, ] <- this_tau_zeta
      betas[keep_index, , ] <- this_beta
      gammas[keep_index, , ] <- this_gamma
      sigmasq_m[keep_index, ] <- this_sigmasq_m
      sigmasq_l[keep_index, ] <- this_sigmasq_l
      sigmasq_pB[keep_index] <- this_sigmasq_pB
      sigmasq_pP[keep_index] <- this_sigmasq_pP
      deltas[keep_index, ] <- this_delta
      zetas[keep_index, ] <- this_zeta
      Us[keep_index, , ] <- this_U
      Vs[keep_index, , ] <- this_V
      vs[keep_index, ] <- this_v
      omegas[keep_index, ] <- this_omega
      zs[keep_index, ] <- this_z
      thetas[keep_index, ] <- this_theta
      pis[keep_index, ] <- this_pi
      pjs[keep_index, ] <- this_pj
      rU[keep_index] <- this_ru
      rV[keep_index] <- this_rv
      
      if (any_X_miss) {
        for (jj in 1 : sum(pB)) {
          Xs[[jj]][keep_index, ] <- this_X[miss_X_ind[[jj]], jj]
        }
      }
      if (any_W_miss) {
        for (jj in 1 : sum(pP)) {
          Ws[[jj]][keep_index, ] <- this_W[miss_W_ind[[jj]], jj]
        }
      }
      
      # Increasing the index by 1.
      keep_index <- keep_index + 1
    }
    
  }
  
  
  # ---------- PART 5 ---------- #
  # Returning the results:
  
  r <- list(Ls = Ls, pL1s = pL1s, mod_pL1s = mod_pL1s,
            lambdas = lambdas, taus_beta = taus_beta,
            taus_gamma = taus_gamma, taus_lambda = taus_lambda,
            taus_delta = taus_delta, taus_zeta = taus_zeta,
            betas = betas, gammas = gammas, sigmasq_m = sigmasq_m,
            sigmasq_l = sigmasq_l, sigmasq_pB = sigmasq_pB,
            sigmasq_pP = sigmasq_pP, deltas = deltas, zetas = zetas,
            Us = Us, Vs = Vs, vs = vs, omegas = omegas, zs = zs,
            thetas = thetas, pis = pis, pjs = pjs, rU = rU, rV = rV,
            pi_accepted = pi_accepted / (Nsims * thin + burn),
            pj_accepted = pj_accepted / (Nsims * thin + burn), 
            ru_accepted = ru_accepted / (Nsims * thin + burn),
            rv_accepted = rv_accepted / (Nsims * thin + burn))
  
  if (any_X_miss) {
    r$Xs <- Xs
  }
  if (any_W_miss) {
    r$Ws <- Ws
  }
  
  return(r)
  
}