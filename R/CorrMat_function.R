#' Specifying a correlation matrix for species.
#' 
#' @param n Dimension of the correlation matrix
#' @param groups List of indices that belong to the same correlation class. The
#' correlation among indices in the same element of the list will be non-zero.
#' @param values Vector of correlation values for the indices belonging to the
#' different elements in the list groups.
#' 
#' @return Correlation matrix.
#' 
#' @export
#' 
CorrMat <- function(n, groups = NULL, values = NULL) {
  
  Cmat <- diag(n)
  for (gg in 1 : length(groups)) {
    Cmat[groups[[gg]], groups[[gg]]] <- values[gg]
  }
  diag(Cmat) <- 1
  return(Cmat)
  
}