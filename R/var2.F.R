#' Two-Sample F-Test for Variance
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \sigma_x^2 \left\lbrace =,\geq,\leq \right\rbrace \sigma_y^2\quad vs\quad H_1 : \sigma_x^2 \left\lbrace \neq,<,>\right\rbrace \sigma_y^2}.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param y a length-\eqn{m} data vector.
#' @param alternative specifying the alternative hypothesis.
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
#' var2.F(x, y, alternative="g") ## Ha : var(x) >= var(y)
#' var2.F(x, y, alternative="l") ## Ha : var(x) <= var(y)
#' var2.F(x, y, alternative="t") ## Ha : var(x) =/= var(y)
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(57)  # sample x from N(0,1)
#'   y = rnorm(89)  # sample y from N(0,1)
#'   
#'   counter[i] = ifelse(var2.F(x,y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'var2.F'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{snedecor_statistical_1996}{SHT}
#' 
#' @export
var2.F <- function(x, y, alternative=c("two.sided","less","greater")){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector of 1st class
  check_1d(y)        # univariate vector of 2nd class
  if (missing(alternative)){
    alternative = "two.sided"
  } else {
    if (pracma::strcmp(alternative,"g")){
      alternative = "greater"
    } else if (pracma::strcmp(alternative,"t")){
      alternative = "two.sided"
    } else if (pracma::strcmp(alternative,"l")){
      alternative = "less"
    }
    alternative = match.arg(alternative)
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n = length(x)
  m = length(y)
  varx = stats::var(x)
  vary = stats::var(y)
  thestat = (varx/vary)
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  if (pracma::strcmp(alternative,"two.sided")){
    tmpval = stats::pf(thestat,(n-1),(m-1))
    pvalue = 2*min(tmpval, 1-tmpval)
    Ha     = "two true variances are different."
  } else if (pracma::strcmp(alternative,"less")){
    pvalue = stats::pf(thestat,(n-1),(m-1),lower.tail = TRUE)
    Ha     = "true variance of x is smaller than true variance of y."
  } else if (pracma::strcmp(alternative,"greater")){
    pvalue = stats::pf(thestat,(n-1),(m-1),lower.tail = FALSE)
    Ha     = "true variance of x is greater than true variance of y."
  }

  ##############################################################
  # REPORT
  hname  = "Two-Sample F Test for Variance."
  DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
    