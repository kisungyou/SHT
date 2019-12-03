#' One-sample Test for Mean Vector by Dempster (1958, 1960)
#' 
#' Given a multivariate sample \eqn{X} and hypothesized mean \eqn{\mu_0}, it tests
#' \deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
#' using the procedure by Dempster (1958, 1960).
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
#' mean1.1958Dempster(smallX) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=50)
#'   counter[i] = ifelse(mean1.1958Dempster(X)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean1.1958Dempster'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{dempster_high_1958}{SHT}
#' 
#' \insertRef{dempster_significance_1960}{SHT}
#' 
#' @author Kisung You
#' @export
mean1.1958Dempster <- function(X, mu0=rep(0,ncol(X))){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_1d(mu0)      
  if (length(mu0)!=ncol(X)){
    stop("* mean1.1958Dempster : mu0 does not have consistent size as data.")
  }
  Xnew = aux_minusvec(X,mu0) # now the problem becomes testing agains (0,0,...,0)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  N      = nrow(Xnew)
  n      = (N-1)
  p      = ncol(Xnew)
  
  xbar   = as.vector(colMeans(Xnew))
  S      = cov(Xnew)
  trS    = aux_trace(S)
  
  a1 = aux_trace(S)/p
  a2 = ((n^2)/((n-1)*(n+2)*p))*(aux_trace(S%*%S)-((trS^2)/n))
  rr = p*(a1^2)/a2
  
  df1my = floor(rr)
  df2my = floor(n*rr)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = N*sum(xbar*xbar)/trS
  pvalue  = pf(thestat,df1my,df2my,lower.tail = FALSE) # reject if (T2 > thr_alpha)
  
  hname   = "One-sample Test for Mean Vector by Dempster (1958)."
  Ha      = "true mean is different from mu0."
  # if (pvalue < alpha){
  #   conclusion = "Reject Null Hypothesis."
  # } else {
  #   conclusion = "Not Reject Null Hypothesis."
  # }
  # 
  
  DNAME = deparse(substitute(X)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}