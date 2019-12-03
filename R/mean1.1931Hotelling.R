#' One-sample Hotelling's T-squared Test for Multivariate Mean
#' 
#' Given a multivariate sample \eqn{X} and hypothesized mean \eqn{\mu_0}, it tests
#' \deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
#' using the procedure by Hotelling (1931).
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param mu0 a length-\eqn{p} mean vector of interest.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' ## CRAN-purpose small example
#' smallX = matrix(rnorm(10*3),ncol=3)
#' mean1.1931Hotelling(smallX) # run the test
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=5)
#'   counter[i] = ifelse(mean1.1931Hotelling(X)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean1.1931Hotelling'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{hotelling_generalization_1931}{SHT}
#' 
#' @export
mean1.1931Hotelling <- function(X, mu0=rep(0,ncol(X))){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_1d(mu0)      
  if (length(mu0)!=ncol(X)){
    stop("* mean1.1931Hotelling : mu0 does not have consistent size as data.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n      = nrow(X)
  p      = ncol(X)
  xbar   = colMeans(X)
  Sigma  = cov(X)   # sample covariance
  SigBar = Sigma/n  # divided by n
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  vecdiff = as.vector(xbar)-as.vector(mu0)
  t2      = sum(as.vector(solve(SigBar, vecdiff))*vecdiff) # test statistic
  t2adj   = ((n-p)*t2/(p*(n-1)))
  pvalue  = pf(t2adj,p,(n-p),lower.tail = FALSE)
  
  # if (pvalue < alpha){
  #   conclusion = "Reject Null Hypothesis"
  # } else {
  #   conclusion = "Not Reject Null Hypothesis"
  # }

  ##############################################################
  # REPORT
  hname  = "One-sample Hotelling's T-squared Test"
  Ha     = "true mean is different from mu0."
  
  DNAME = deparse(substitute(X)) # borrowed from HDtest
  names(t2) = "statistic"
  res   = list(statistic=t2, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}