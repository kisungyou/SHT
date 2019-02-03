#' Analysis of Variance for Homogeneity of Means
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \mu_1^2 = \cdots \mu_k^2\quad vs\quad H_1 : \textrm{at least one equality does not hold.}}
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample vector.
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
#' ## test when k=5 (samples)
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   mylist = list()
#'   for (j in 1:5){
#'      mylist[[j]] = rnorm(50)   
#'   }
#'   
#'   counter[i] = ifelse(meank.anova(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'meank.anova'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @export
meank.anova <- function(dlist, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 
  check_alpha(alpha)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY FOR USING ANOVA
  K = length(dlist)
  if (K < 2){
    stop("* meank.anova : we need at least 2 groups of data.")
  }
  labellist = list()
  for (i in 1:K){
    labellist[[i]] = rep(i,length(dlist[[i]]))
  }
  
  data  = unlist(dlist)
  group = as.factor(unlist(labellist))
  
  ##############################################################
  # COMPUTATION : USE AOV INTERFACE
  aovout = unlist(summary(aov(data~group)))
  
  ##############################################################
  # REPORT
  hname   = "Analysis of Variance for Homogeneity of Means"
  Ha      = "at least one of equalities does not hold."
  thestat = as.double(aovout[7])
  pvalue  = as.double(aovout[9])
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }

  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}

