#' Analysis of Variance for Equality of Means
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \mu_1^2 = \cdots \mu_k^2\quad vs\quad H_1 : \textrm{at least one equality does not hold.}}
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample vector.
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
#' cat(paste("\n* Example for 'meank.anova'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @export
meank.anova <- function(dlist){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 

  ##############################################################
  # COMPUTATION : PRELIMINARY FOR USING ANOVA
  K = length(dlist)
  if (K < 2){
    stop("* meank.anova : we need at least 2 sets of data.")
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
  hname   = "Analysis of Variance for Equality of Means"
  Ha      = "at least one of equalities does not hold."
  thestat = as.double(aovout[7])
  pvalue  = as.double(aovout[9])
  # if (pvalue < alpha){
  #   conclusion = "Reject Null Hypothesis"
  # } else {
  #   conclusion = "Not Reject Null Hypothesis"
  # }

  DNAME = deparse(substitute(dlist))
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}

