#' Two-Sample F Test for Variance
#' 
#' 
#' 
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
}
    