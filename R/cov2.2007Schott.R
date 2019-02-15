#' Two-Sample Test for High-Dimensional Covariances by Schott (2007)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_y\quad vs\quad H_1 : \Sigma_x \neq \Sigma_y}
#' using the procedure by Schott (2007).
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value \eqn{P(H_0|H_1)} under current setting.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(cov2.2007Schott(X, Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov2.2007Schott'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{schott_test_2007}{SHT}
#' 
#' @export
cov2.2007Schott <- function(X, Y){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* cov2.2012LC : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # PARAMETERS AND PRELIMINARY : SRIVASTAVA AND YANAGIHARA (2010)
  N1 = nrow(X); n1 = N1-1
  N2 = nrow(Y); n2 = N2-1
  m  = ncol(X)
  n  = (n1+n2)
  
  V1 = cov(X)*n1; S1 = V1/n1
  V2 = cov(Y)*n2; S2 = V2/n2
  V  = (V1+V2)
  
  a21 = (1.0/(m*(n1-1)*(n1+2)))* (sum(diag(V1%*%V1)) - (1/n1)*(sum(diag(V1))^2))
  a22 = (1.0/(m*(n2-1)*(n2+2)))* (sum(diag(V2%*%V2)) - (1/n2)*(sum(diag(V2))^2))
  a2  = (1/((n-1)*(n+2)*m))* (sum(diag(V%*%V)) - (1/n)*(sum(diag(V))^2))
  
  ##############################################################
  # MAIN COMPONENTS
  thestat = ((n1*n2)/(2*a2*(n1+n2)))*(a21+a22-(2/m)*sum(diag(S1%*%S2)))
  pvalue  = stats::pnorm(thestat, lower.tail = FALSE)
  
  ##############################################################
  # FINALE
  hname   = "Two-Sample Test for High-Dimensional Covariances by Schott (2007)"
  Ha      = "two covariances are not equal."
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
