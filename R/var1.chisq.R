#' One-Sample Chi-Square Test for Variance
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : \sigma_x^2 \left\lbrace =,\geq,\leq \right\rbrace \sigma_0^2 \quad vs\quad H_1 : \sigma_x^2 \left\lbrace \neq,<,>\right\rbrace \sigma_0^2}.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param var0 hypothesized variance \eqn{\sigma_0^2}.
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
#' var1.chisq(x, alternative="g") ## Ha : var(x) >= 1
#' var1.chisq(x, alternative="l") ## Ha : var(x) <= 1
#' var1.chisq(x, alternative="t") ## Ha : var(x) =/=1
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(50)  # sample x from N(0,1)
#'   
#'   counter[i] = ifelse(var1.chisq(x,var0=1)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'var1.chisq'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{snedecor_statistical_1996}{SHT}
#' 
#' 
#' @export
var1.chisq <- function(x, var0=1, alternative=c("two.sided","less","greater")){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  check_number(var0) # number to be compared
  if (var0<=0){
    stop("* var1.chisq : var0 should be a nonnegative real number.")
  }
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
  n     = length(x)
  varx  = as.double(stats::var(x))
  xstat = ((n-1)*varx/var0)
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  if (alternative=="two.sided"){
    tmpval = pchisq(xstat, (n-1))
    pvalue = 2*min(tmpval, 1-tmpval)
    Ha     = paste("true variance is different from ",var0,".",sep="")
  } else if (alternative=="less"){
    pvalue = pchisq(xstat, (n-1), lower.tail = TRUE)
    Ha     = paste("true variance is less than ",var0,".",sep="")
  } else if (alternative=="greater"){
    pvalue = pchisq(xstat, (n-1), lower.tail = FALSE)
    Ha     = paste("true variance is greater than ",var0,".",sep="")
  }
  
  ##############################################################
  # REPORT
  hname  = "One-Sample Chi-Square Test for Variance."
  thestat = xstat
  
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}