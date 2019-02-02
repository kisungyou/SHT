#' One-Sample Test for High-Dimensional Mean by Bai and Saranadasa (1996)
#' 
#' Given two multivariate data \eqn{X} and hypothesized mean \eqn{\mu_0}, it tests
#' \deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
#' using the procedure by Bai and Saranadasa (1996).
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param mu0 a length-\eqn{p} mean vector of interest.
#' @param alpha significance level.
#' 
#' @return a (list) object of \code{S3} class \code{hypothesis} containing: \describe{
#' \item{method}{name of the test.}
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under current setting.}
#' \item{significance}{a user-specified significance level.}
#' \item{alternative}{alternative hypothesis.}
#' \item{conclusion}{conclusion by \eqn{p}-value decision rule.}
#' }
#' 
#' @examples 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=25)
#'   counter[i] = ifelse(mean1.1996Bai(X)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean1.1996Bai'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{bai_high_1996}{SHT}
#' 
#' @export
mean1.1996Bai <- function(X, mu0=rep(0,ncol(X)), alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_1d(mu0)      
  check_alpha(alpha)
  if (length(mu0)!=ncol(X)){
    stop("* mean1.1996Bai : mu0 does not have consistent size as data.")
  }
  Xnew = aux_minusvec(X,mu0) # now the problem becomes testing agains (0,0,...,0)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  N      = nrow(Xnew)
  n      = (N-1)
  p      = ncol(Xnew)
  
  xbar   = colMeans(Xnew)
  S      = stats::cov(Xnew)    # sample covariance
  trS    = aux_trace(S)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  term1   = N*sum(xbar*xbar) - trS
  term2   = sqrt((2*n*(n+1)/((n-1)*(n+1)))*(aux_trace(S%*%S) - (trS^2)/n))
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "One-Sample Test for High-Dimensional Mean by Bai and Saranadasa (1996)."
  Ha      = "true mean is different from mu0."
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis."
  } else {
    conclusion = "Not Reject Null Hypothesis."
  }
  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
  
}