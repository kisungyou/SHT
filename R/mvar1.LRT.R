#' One-sample Simultaneous Likelihood Ratio Test of Mean and Variance
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \mu_x = \mu_0, \sigma_x^2 = \sigma_0^2 \quad vs \quad H_1 : \textrm{ not } H_0}
#' using likelihood ratio test.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param mu0 hypothesized mean \eqn{\mu_0}.
#' @param var0 hypothesized variance \eqn{\sigma_0^2}.
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
#' mvar1.LRT(rnorm(10))
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(100)  # sample x from N(0,1)
#'   
#'   counter[i] = ifelse(mvar1.LRT(x)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mvar1.LRT'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @export
mvar1.LRT <- function(x, mu0=0, var0=1){
  ##############################################################
  # Preprocessing & Parameters
  DNAME = deparse(substitute(x))
  check_1d(x)        # univariate vector of 1st class
  check_number(mu0)  # check : univariate mean
  check_number(var0) # check : univariate variance
  if (var0<=0){
    stop("* mvar1.LRT : var0 should be a nonnegative real number.")
  }
  
  ##############################################################
  # Computation
  n    = length(x)
  xbar = base::mean(x)
  s2   = sum((x-xbar)^2)/n
  
  y1 = n*s2/var0
  y2 = n*((xbar-mu0)^2)/var0
  
  #   statistic & p-value
  loglbd  = (n/2)*log(exp(1)/n) + (n/2)*log(y1) - ((y1+y2)/2)
  thestat = -2*loglbd
  pvalue  = pchisq(thestat, df=2, lower.tail = FALSE)
  
  ##############################################################
  # REPORT
  hname  = "One-sample Simultaneous Likelihood Ratio Test of Mean and Variance."
  Ha     = "true mean and variance of x are different from mu0 and var0."
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}