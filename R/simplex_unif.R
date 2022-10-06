#' Classical Test of Uniformity
#' 
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param method name of the method to be used, including\describe{
#' \item{LRT}{likelihood-ratio test with the Dirichlet distribution.}
#' \item{LRTsym}{likelihood-ratio test using the symmetric Dirichlet distribution (default).}
#' }
#' 
#' @concept simplex
#' @export
simplex.uniform <- function(X, method){
  # -------------------------------------------------------------------------
  # Input
  # data
  X = check_simplex(X, "simplex.uniform")
  
  # choice of a method
  method.all = c("lrt","lrtsym")
  if (missing(method)){
    method.now = "lrtsym"
  } else {
    method.now = match.arg(tolower(method), method.all)
  }
  
  
}