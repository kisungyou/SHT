#' Two-Sample Student's t-test for Univariate Means
#' 
#' Given two univariate samples \eqn{x} and \eqn{y}, it tests
#' \deqn{H_0 : \mu_x^2 \left\lbrace =,\geq,\leq \right\rbrace \mu_y^2\quad vs\quad H_1 : \mu_x^2 \left\lbrace \neq,<,>\right\rbrace \mu_y^2}
#' using the procedure by Student (1908) and Welch (1947).
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
#'   counter[i] = ifelse(mean2.ttest(x,y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.ttest'\n\n",
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
#' \insertRef{welch_generalization_1947}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.ttest <- function(x, y, alternative=c("two.sided","less","greater"), alpha=0.05, paired=FALSE, var.equal=FALSE){
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
  # CASE 1 : PAIRED/DEPENDENT T-TEST
  if (paired==TRUE){
    output = mean2.Student.paired(x,y,alternative,alpha)
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
    if (pvalue < alpha){
      conclusion = "Reject Null Hypothesis"
    } else {
      conclusion = "Not Reject Null Hypothesis"
    }
    output = hypothesis(hname, t, alpha,
                        pvalue, Ha, 
                        conclusion)
  }
  ##############################################################
  # REPORT
  return(output)
}


# paired case -------------------------------------------------------------
#' @keywords internal
#' @noRd
mean2.Student.paired <- function(x, y, paired.alternative, paired.alpha){
  if (length(x)!=length(y)){
    stop("* mean2.Student : for paired t-test, length of x and y should be identical.")
  }
  diff   = (x-y)
  tmpout = mean1.Student(diff,mu0=0,alternative=paired.alternative,alpha=paired.alpha)
  
  hname  = "Two-sample Student's t-test for Paired/Dependent Data"
  if (paired.alternative=="two.sided"){
    Ha = "true mean difference is not equal to 0."
  } else if (paired.alternative=="less"){
    Ha = "true mean difference is less than 0."
  } else if (paired.alternative=="greater"){
    Ha = "true mean difference is greater than 0."
  }
  
  output = hypothesis(hname, tmpout$statistic, tmpout$significance,
                      tmpout$p.value, Ha, tmpout$conclusion)
  return(output)
}
