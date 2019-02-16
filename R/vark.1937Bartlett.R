#' Bartlett's Test for Homogeneity of Variance
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \sigma_1^2 = \cdots \sigma_k^2\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Bartlett (1937).
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample vector.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value \eqn{P(H_0|H_1)} under current setting.}
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
#'   counter[i] = ifelse(vark.1937Bartlett(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'vark.1937Bartlett'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{bartlett_properties_1937}{SHT}
#' 
#' @export
vark.1937Bartlett <- function(dlist){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  k        = length(dlist)
  vec_n    = unlist(lapply(dlist, length))
  vec_Si2  = unlist(lapply(dlist, stats::var))
  
  N        = sum(vec_n)
  Sp2      = (sum((vec_n-1)*vec_Si2))/(N-k)
  
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  term1 = (((N-k)*log(Sp2)) - (sum((vec_n-1)*log(vec_Si2))))
  term2 = (((sum(1/(vec_n-1)))-(1/(N-k)))/(3*(k-1))) + 1
  
  thestat = (term1/term2)
  pvalue  = stats::pchisq(thestat,(k-1),lower.tail = FALSE)
  Ha      = "at least one of equalities does not hold."
  
  ##############################################################
  # REPORT
  hname  = "Bartlett\'s Test for Homogeneity of Variance."
  DNAME = deparse(substitute(dlist)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}