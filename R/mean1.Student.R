#' One-sample Student's t-test
#' 
#' 
#' 
#' 
#' 
#' 
#' @export
mean1.Student <- function(x, mu0=0, alternative=c("two.sided","less","greater"), alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  check_number(mu0)  # number to be compared
  check_alpha(alpha) # significance level
  if (missing(alternative)){
    alternative = "two.sided"
  } else {
    alternative = match.arg(alternative)
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n    = length(x)
  xbar = mean(x)
  sd   = sd(x)    
  t    = (xbar-mu0)/(sd/sqrt(n))
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  if (alternative=="two.sided"){
    pvalue = 2*pt(abs(t),(n-1),lower.tail = FALSE)
    Ha     = paste("true mean is different from ",mu0,".",sep="")
  } else if (alternative=="less"){
    pvalue = pt(t,(n-1),lower.tail = TRUE)
    Ha     = paste("true mean is less than ",mu0,".",sep="")
  } else if (alternative=="greater"){
    pvalue = pt(t,(n-1),lower.tail = FALSE)
    Ha     = paste("true mean is greater than ",mu0,".",sep="")
  }
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }
  
  ##############################################################
  # REPORT
  hname  = "One-sample Student's t-test"
  output = hypothesis(hname, t, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}