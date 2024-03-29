% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UpdExtraVar_function.R
\name{UpdExtraVar}
\alias{UpdExtraVar}
\title{Updating the additional variance component for coefficients.}
\usage{
UpdExtraVar(mod_coef, shr_var, prior_spec)
}
\arguments{
\item{mod_coef}{Values for coefficients whose extra variance will be
updated. These can be the coefficients of the physical traits (matrix format
with rows corresponding to the trait, and columns to the latent factors),
coefficients of the network model (vector with length equal to the number of
latent factors), coefficients of the models for probability of observing a
certain species (vector).}

\item{shr_var}{Value of the variances across the H factor components.}

\item{prior_spec}{Values of parameters in the prior specification.}
}
\description{
Function that updates the extra variance, and can be used for coefficients
across multiple models.
}
