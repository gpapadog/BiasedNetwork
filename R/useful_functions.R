#' Classic expit function.
expit <- function(x) {
  return(exp(x) / (1 + exp(x)))
}

logit <- function(x) {
  return(log(x / (1 - x)))
}