#' Getting omega from v in the increasing shrinkage prior.
#'
#' @param v_val Vector of length equal to the value of H we use including
#' values between 0 and 1.
#' 
OmegaFromV <- function(v_val) {
  
  Hval <- length(v_val)
  
  omega <- v_val
  omega[2 : Hval] <- sapply(2 : Hval, function(x) v_val[x] * prod(1 - v_val[1 : (x - 1)]))

  return(omega)
  
}