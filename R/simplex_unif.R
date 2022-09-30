#' Classical Test of Uniformity
#' 
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param method 
#' 
#' @concept simplex
#' @export
simplex.uniform <- function(X, method=c("LRTsym","LRT")){
  # -------------------------------------------------------------------------
  # Input
  X = check_simplex(X, "simplex.uniform")
}