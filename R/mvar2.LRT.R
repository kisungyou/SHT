#' Two-sample Simultaneous Likelihood Ratio Test of Mean and Variance
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \mu_x = \mu_y, \sigma_x^2 = \sigma_y^2 \quad vs \quad H_1 : \textrm{ not } H_0}
#' using classical likelihood ratio test.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param y a length-\eqn{m} data vector.
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
#' x = rnorm(10)
#' y = rnorm(10)
#' mvar2.LRT(x, y)
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(100)  # sample x from N(0,1)
#'   y = rnorm(100)  # sample y from N(0,1)
#'   
#'   counter[i] = ifelse(mvar2.LRT(x,y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mvar2.LRT'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @export
mvar2.LRT <- function(x, y){
  ##############################################################
  # Preprocessing & Parameters
  DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
  check_1d(x)        # univariate vector of 1st class
  check_1d(y)        # univariate vector of 2nd class
  
  m = length(x)
  n = length(y)

  ##############################################################
  # Computation
  xbar = base::mean(x)
  ybar = base::mean(y)
  
  sig2 = sum((x-xbar)^2)/m
  tau2 = sum((y-ybar)^2)/n
  mu0  = (m*xbar + n*ybar)/(m+n)
  sig0 = (sum((x-mu0)^2) + sum((y-mu0)^2))/(m+n)
  
  #   statistic & p-value
  loglambda = (m/2)*log(sig2) + (n/2)*log(tau2) - ((m+n)/2)*log(sig0)
  thestat   = -2.0*loglambda
  pvalue    = pchisq(thestat, df=2, lower.tail = FALSE)

  ##############################################################
  # REPORT
  hname  = "Two-sample Simultaneous Likelihood Ratio Test of Mean and Variance."
  Ha     = "true mean and variance of x are different from those of y."
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
