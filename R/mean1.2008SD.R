#' One-Sample Test for High-Dimensional Mean by Srivastava and Du (2008)
#' 
#' Given two multivariate data \eqn{X} and hypothesized mean \eqn{\mu_0}, it tests
#' \deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
#' using the procedure by Srivastava and Du (2008).
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param mu0 a length-\eqn{p} mean vector of interest.
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
#'   X = matrix(rnorm(50*5), ncol=5)
#'   counter[i] = ifelse(mean1.2008SD(X)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean1.2008SD'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{srivastava_test_2008}{SHT}
#' 
#' @author Kisung You
#' @export
mean1.2008SD <- function(X, mu0=rep(0,ncol(X))){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_1d(mu0)      
  if (length(mu0)!=ncol(X)){
    stop("* mean1.2008SD : mu0 does not have consistent size as data.")
  }
  Xnew = aux_minusvec(X,mu0) # now the problem becomes testing agains (0,0,...,0)
  
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  N      = nrow(Xnew)
  n      = (N-1)
  p      = ncol(Xnew)
  xbar   = colMeans(Xnew)
  S      = stats::cov(Xnew)    # sample covariance
  R      = stats::cov2cor(S)
  trR2   = aux_trace((R%*%R))
  Dsinv  = diag(1/aux_adjustvec(diag(S)))
  cpn    = (1+(trR2/(p^(3/2))))
  
  term1 = (N*aux_quadform(Dsinv,xbar))-((n*p)/(n-2))
  term2 = sqrt((2*(trR2-((p^2)/n)))*cpn)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "One-Sample Test for High-Dimensional Mean by Srivastava and Du (2008)."
  Ha      = "true mean is different from mu0."
  # if (pvalue < alpha){
  #   conclusion = "Reject Null Hypothesis."
  # } else {
  #   conclusion = "Not Reject Null Hypothesis."
  # }

  DNAME = deparse(substitute(X)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}