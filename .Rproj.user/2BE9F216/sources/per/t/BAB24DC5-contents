#' Two-sample Student's t-test
#' 
#' 
#' 
#' @param x a vector of length \eqn{m}.
#' @param y a vector of length \eqn{n}.
#' @param alternative
#' @param alpha significance level, default set as 0.05.
#' @param paired a logical indicating to use a paired \eqn{t}-test.
#' @param var.equal a logical indicating to use Welch's variant.
#' 
#' 
#' 
#' @export
mean2.Student <- function(x, y, alternative=c("two.sided","less","greater"), alpha=0.05, paired=FALSE, var.equal=FALSE){
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
      Ha     = "two means are different."
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