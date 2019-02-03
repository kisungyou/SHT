#' Two-Sample F-Test for Variance
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \sigma_x^2 \left\lbrace =,\geq,\leq \right\rbrace \sigma_y^2\quad vs\quad H_1 : \sigma_x^2 \left\lbrace \neq,<,>\right\rbrace \sigma_y^2}.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param y a length-\eqn{m} data vector.
#' @param alternative specifying the alternative hypothesis.
#' @param alpha significance level.
#' 
#' @return a (list) object of \code{S3} class \code{hypothesis} containing: \describe{
#' \item{method}{name of the test.}
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under current setting.}
#' \item{significance}{a user-specified significance level.}
#' \item{alternative}{alternative hypothesis.}
#' \item{conclusion}{conclusion by \eqn{p}-value decision rule.}
#' }
#' 
#' @examples 
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
#' cat(paste("\n* Example for 'var2.F'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{snedecor_statistical_1996}{SHT}
#' 
#' @export
var2.F <- function(x, y, alternative=c("two.sided","less","greater"), alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector of 1st class
  check_1d(y)        # univariate vector of 2nd class
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
  n = length(x)
  m = length(y)
  varx = aux_var(x)
  vary = aux_var(y)
  thestat = (varx/vary)
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  if (alternative=="two.sided"){
    tmpval = stats::pf(thestat,(n-1),(m-1))
    pvalue = 2*min(tmpval, 1-tmpval)
    Ha     = "two true variances are different."
  } else if (alternative=="less"){
    pvalue = stats::pf(thestat,(n-1),(m-1),lower.tail = TRUE)
    Ha     = "true variance of x is smaller than true variance of y."
  } else if (alternative=="greater"){
    pvalue = stats::pf(thestat,(n-1),(m-1),lower.tail = FALSE)
    Ha     = "true variance of x is greater than true variance of y."
  }
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }
  
  ##############################################################
  # REPORT
  hname  = "Two-Sample F Test for Variance"
  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}
    