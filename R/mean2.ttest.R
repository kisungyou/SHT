#' Two-sample Student's t-test for Univariate Means
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \mu_x^2 \left\lbrace =,\geq,\leq \right\rbrace \mu_y^2\quad vs\quad H_1 : \mu_x^2 \left\lbrace \neq,<,>\right\rbrace \mu_y^2}
#' using the procedure by Student (1908) and Welch (1947).
#' 
#' @param x a length-\eqn{n} data vector.
#' @param y a length-\eqn{m} data vector.
#' @param alternative specifying the alternative hypothesis.
#' @param paired a logical; whether consider two samples as paired.
#' @param var.equal a logical; if \code{FALSE}, use Welch's correction.
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
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   x = rnorm(57)  # sample x from N(0,1)
#'   y = rnorm(89)  # sample y from N(0,1)
#'   
#'   counter[i] = ifelse(mean2.ttest(x,y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.ttest'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{student_probable_1908}{SHT}
#' 
#' \insertRef{student_probable_1908a}{SHT}
#' 
#' \insertRef{welch_generalization_1947}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.ttest <- function(x, y, alternative=c("two.sided","less","greater"), paired=FALSE, var.equal=FALSE){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector of 1st class
  check_1d(y)        # univariate vector of 2nd class
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
  # CASE 1 : PAIRED/DEPENDENT T-TEST
  if (paired==TRUE){
    res = mean2.Student.paired(x,y,alternative)
  } else {
  ##############################################################
  # CASE 2 : INDEPENDENT T-TEST
    # 2-1. Common Parameters and Computations
    n1 = length(x)
    n2 = length(y)
    s1 = sd(x)
    s2 = sd(y)
    xbar1 = mean(x)
    xbar2 = mean(y)
    
    # 2-2. branching
    if (var.equal==TRUE){
      hname  = "Student's t-test for Independent Samples with Equal Variance Assumption"  
      
      df = (n1+n2-2)
      sp = sqrt(((n1-1)*(s1^2) + (n2-1)*(s2^2))/(n1+n2-2))
      
      t  = ((xbar1-xbar2)/(sp*sqrt((1/n1)+(1/n2))))
    } else if (var.equal==FALSE){
      hname  = "Student's t-test for Independent Samples with Unequal Variance Assumption"  
      
      df1 = ((((s1^2)/n1) + ((s2^2)/n2))^2)
      df2 = (((((s1^2)/n1)^2)/(n1-1)) + ((((s2^2)/n2)^2)/(n2-1)))
      df  = (df1/df2)
      
      sp  = sqrt(((s1^2)/n1) + ((s2^2)/n2))
      t   = ((xbar1-xbar2)/sp)
    }
    
    # 2-3. hypothesis and determination
    if (alternative=="two.sided"){
      pvalue = 2*pt(abs(t),df,lower.tail = FALSE)
      Ha     = "two true means are different."
    } else if (alternative=="less"){
      pvalue = pt(t,df,lower.tail = TRUE)
      Ha     = "true mean of x is smaller than true mean of y."
    } else if (alternative=="greater"){
      pvalue = pt(t,df,lower.tail = FALSE)
      Ha     = "true mean of x is greater than true mean of y."
    }
    # if (pvalue < alpha){
    #   conclusion = "Reject Null Hypothesis"
    # } else {
    #   conclusion = "Not Reject Null Hypothesis"
    # }
    
    
    thestat = t
    DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
    names(thestat) = "t"
    res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
    class(res) = "htest"
  }
  ##############################################################
  # REPORT
  return(res)
}


# paired case -------------------------------------------------------------
#' @keywords internal
#' @noRd
mean2.Student.paired <- function(x, y, paired.alternative){
  if (length(x)!=length(y)){
    stop("* mean2.Student : for paired t-test, length of x and y should be identical.")
  }
  diff   = (x-y)
  tmpout = mean1.ttest(diff,mu0=0,alternative=paired.alternative)
  
  hname  = "Two-sample Student's t-test for Paired/Dependent Data."
  if (paired.alternative=="two.sided"){
    Ha = "true mean difference is not equal to 0."
  } else if (paired.alternative=="less"){
    Ha = "true mean difference is less than 0."
  } else if (paired.alternative=="greater"){
    Ha = "true mean difference is greater than 0."
  }
  
  thestat = tmpout$statistic
  DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
  names(thestat) = "t"
  res   = list(statistic=thestat, p.value=tmpout$p.value, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
