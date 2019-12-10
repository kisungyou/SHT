#' Two-sample Simultaneous Test of Mean and Variance by Perng and Littell (1976)
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \mu_x = \mu_y, \sigma_x^2 = \sigma_y^2 \quad vs \quad H_1 : \textrm{ not } H_0}
#' using Fisher's method of merging two \eqn{p}-values.
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
#' mvar2.1976PL(x, y)
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(100)  # sample x from N(0,1)
#'   y = rnorm(100)  # sample y from N(0,1)
#'   
#'   counter[i] = ifelse(mvar2.1976PL(x,y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mvar2.1976PL'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{perng_test_1976}{SHT}
#' 
#' @export
mvar2.1976PL <- function(x, y){
  ##############################################################
  # Preprocessing & Parameters
  DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
  check_1d(x)        # univariate vector of 1st class
  check_1d(y)        # univariate vector of 2nd class
  
  ##############################################################
  # Merge ! 
  pval1 = SHT::mean2.ttest(x, y, var.equal = TRUE)$p.value
  pval2 = SHT::var2.F(x, y)$p.value
  
  #   statistic & p-value
  thestat = -2*log(pval1) - 2*log(pval2)
  pvalue  = pchisq(thestat, df=4, lower.tail = FALSE)
  
  ##############################################################
  # REPORT
  hname  = "Two-sample Simultaneous Test of Mean and Variance by Perng and Littell (1976)."
  Ha     = "true mean and variance of x are different from those of y."
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}