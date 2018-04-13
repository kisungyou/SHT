#' One-Sample Chi-Square Test for Variance
#' 
#' 
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
    alternative = match.arg(alternative)
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n     = length(x)
  varx  = aux_var(x)
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
  hname  = "One-Sample Chi-Square Test for Variance"
  output = hypothesis(hname, t, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}