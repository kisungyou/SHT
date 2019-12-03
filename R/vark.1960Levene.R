#' Levene's Test for Homogeneity of Variance
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \sigma_1^2 = \cdots \sigma_k^2\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Levene (1960).
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
#' ## CRAN-purpose small example
#' small1d = list()
#' for (i in 1:5){ # k=5 sample
#'   small1d[[i]] = rnorm(20)
#' }
#' vark.1960Levene(small1d) # run the test
#' 
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
#'   counter[i] = ifelse(vark.1960Levene(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'vark.1960Levene'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{levene_robust_1960}{SHT}
#' 
#' @export
vark.1960Levene <- function(dlist){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 

  ##############################################################
  # COMPUTATION : PRELIMINARY
  k        = length(dlist)
  ni       = unlist(lapply(dlist, length))
  N        = sum(ni)
  
  ##############################################################
  # COMPUTATION : LEVENE'S
  zlist = list()
  for (i in 1:k){
    tgt  = dlist[[i]]
    ybar = mean(tgt)      # this must be changed for Brown-Forsyth
    zlist[[i]] = abs(tgt-ybar)
  }
  z_meanvec = unlist(lapply(zlist, mean))
  z_meanall = mean(unlist(zlist))
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  term1 = 0
  for (i in 1:k){
    term1 = term1 + (ni[i])*((z_meanvec[i]-z_meanall)^2)
  }
  term2 = 0
  for (i in 1:k){
    vecdiff = (zlist[[i]]-z_meanvec[i])
    term2   = term2 + sum(vecdiff*vecdiff)
  }
  thestat = (((N-k)*term1)/((k-1)*term2))
  pvalue  = stats::pf(thestat,(k-1),(N-k),lower.tail = FALSE)
  Ha      = "at least one of equalities does not hold."
  
  # if (pvalue < alpha){
  #   conclusion = "Reject Null Hypothesis"
  # } else {
  #   conclusion = "Not Reject Null Hypothesis"
  # }
  # 
  ##############################################################
  # REPORT
  hname  = "Levene\'s Test for Homogeneity of Variance."
  DNAME = deparse(substitute(dlist)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}