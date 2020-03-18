#' One-sample Student's t-test for Univariate Mean
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : \mu_x = \mu_0\quad vs\quad H_1 : \mu_x \neq \mu_0}
#' using the procedure by Student (1908).
#' 
#' @param x a length-\eqn{n} data vector.
#' @param mu0 hypothesized variance \eqn{\sigma_0^2}.
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
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(10)         # sample from N(0,1)
#'   counter[i] = ifelse(mean1.ttest(x)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean1.ttest'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' 
#' @references 
#' \insertRef{student_probable_1908}{SHT}
#' 
#' \insertRef{student_probable_1908a}{SHT}
#' 
#' @export
mean1.ttest <- function(x, mu0=0, alternative=c("two.sided","less","greater")){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  check_number(mu0)  # number to be compared
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

  ##############################################################
  # REPORT
  hname   = "One-sample Student\'s t-test"
  thestat = t
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  names(thestat) = "t"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
