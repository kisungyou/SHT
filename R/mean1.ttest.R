#' One-Sample Student's t-test for Univariate Mean
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
#' using the procedure by Student (1908).
#' 
#' @param x a length-\eqn{n} data vector.
#' @param mu0 hypothesized variance \eqn{\sigma_0^2}.
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
#'   x = rnorm(10)         # sample from N(0,1)
#'   counter[i] = ifelse(mean1.ttest(x)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean1.ttest'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{student_probable_1908}{SHT}
#' 
#' \insertRef{student_probable_1908-1}{SHT}
#' 
#' @export
mean1.ttest <- function(x, mu0=0, alternative=c("two.sided","less","greater"), alpha=0.05){
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
  hname  = "One-Sample Student\'s t-test"
  output = hypothesis(hname, t, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}
