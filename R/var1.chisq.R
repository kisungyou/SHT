#' One-Sample Chi-Square Test for Variance
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : \sigma_x^2 \left\lbrace =,\geq,\leq \right\rbrace \sigma_0^2 \quad vs\quad H_1 : \sigma_x^2 \left\lbrace \neq,<,>\right\rbrace \sigma_0^2}.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param var0 hypothesized variance \eqn{\sigma_0^2}.
#' @param alternative specifying the alternative hypothesis.
#' @param alpha significance level.
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
#'   x = rnorm(50)  # sample x from N(0,1)
#'   
#'   counter[i] = ifelse(var1.chisq(x,var0=1)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'var1.chisq'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{snedecor_statistical_1996}{SHT}
#' 
#' 
#' @export
var1.chisq <- function(x, var0=1, alternative=c("two.sided","less","greater"), alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  check_number(var0) # number to be compared
  if (var0<=0){
    stop("* var1.chisq : var0 should be a nonnegative real number.")
  }
  check_alpha(alpha) # significance level
  if (missing(alternative)){
    alternative = "two.sided"
  } else {
    if (alternative=="g"){
      alternative = "greater"
    } else if (alternative=="t"){
      alternative = "two.sided"
    } else if (alternative=="l"){
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
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
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