#' Two-sample Simultaneous Test of Mean and Variance by Pearson and Neyman (1930)
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \mu_x = \mu_y, \sigma_x^2 = \sigma_y^2 \quad vs \quad H_1 : \textrm{ not } H_0}
#' by approximating the null distribution with Beta distribution using the first two moments matching.
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
#' mvar2.1930PN(x, y)
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(100)  # sample x from N(0,1)
#'   y = rnorm(100)  # sample y from N(0,1)
#'   
#'   counter[i] = ifelse(mvar2.1930PN(x,y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mvar2.1930PN'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#'
#' @export
mvar2.1930PN <- function(x, y){
  ##############################################################
  # Preprocessing & Parameters
  DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
  check_1d(x)        # univariate vector of 1st class
  check_1d(y)        # univariate vector of 2nd class
  
  n = length(x)
  m = length(y)
  xbar = base::mean(x)
  ybar = base::mean(y)
  u    = (n*xbar + m*ybar)/(n+m)
  
  ##############################################################
  # Estimator
  term1  = (n/2)*log(sum((x-xbar)^2)/n) + (m/2)*log(sum((y-ybar)^2)/m)
  term2  = ((m+n)/2)*log((sum((x-u)^2) + sum((y-u)^2))/(n+m))
  lambda = exp(term1-term2)
  
  ##############################################################
  # Moment Matching to Beta Distribution
  a1 = ((n+m)/2)*log(n+m) - (n/2)*log(n) - (m/2)*log(m)
  a2 = lgamma((2*n-1)/2) + lgamma((2*m-1)/2) - lgamma((2*(n+m)-1)/2)
  a3 = lgamma((n+m-1)/2) - lgamma((n-1)/2) - lgamma((m-1)/2)
  a  = exp(a1+a2+a3)
  
  b1 = (n+m)*log(n+m) - n*log(n) - m*log(m)
  b2 = lgamma((3*n-1)/2) + lgamma((3*m-1)/2) - lgamma((3*(n+m)-1)/2)
  b3 = lgamma((n+m-1)/2) - lgamma((n-1)/2) - lgamma((m-1)/2)
  b  = exp(b1+b2+b3) - (a^2)

  p = -(a/b)*(a^2 - a + b)
  q = ((a-1)/b)*(a^2 - a + b)
  
  #   statistic & p-value
  thestat = lambda
  pvalue  = stats::pbeta(lambda, shape1=p, shape2=q, lower.tail=FALSE)
  
  
  ##############################################################
  # REPORT
  hname  = "Two-sample Simultaneous Test of Mean and Variance by Muirhead Approximation (1982)."
  Ha     = "true mean and variance of x are different from those of y."
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}