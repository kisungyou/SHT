#' Test for Equality of Means by Schott (2007)
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \mu_1 = \cdots \mu_k\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Schott (2007). It can be considered as a generalization 
#' of two-sample testing procedure proposed by \code{\link[SHT:mean2.1996BS]{Bai and Saranadasa (1996)}}.
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample matrix of same dimension.
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
#' tinylist = list()
#' for (i in 1:3){ # consider 3-sample case
#'   tinylist[[i]] = matrix(rnorm(10*3),ncol=3)
#' }
#' meank.2007Schott(tinylist)
#' 
#' \donttest{
#' ## test when k=5 samples with (n,p) = (10,50)
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   mylist = list()
#'   for (j in 1:5){
#'      mylist[[j]] = matrix(rnorm(10*5),ncol=5)
#'   }
#'   
#'   counter[i] = ifelse(meank.2007Schott(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'meank.2007Schott'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' 
#' @references 
#' \insertRef{schott_highdimensional_2007}{SHT}
#' 
#' @export
meank.2007Schott <- function(dlist){
  ##############################################################
  # PREPROCESSING AND PARAMETERS
  check_dlistnd(dlist) 
  vec.ni = unlist(lapply(dlist, nrow)) # obs : per-class
  n = sum(vec.ni)                      # obs : total number
  p = ncol(dlist[[1]])
  g = length(dlist)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  vec.E = lapply(dlist, cov)       # per-class covariance
  E = array(0,c(p,p))
  for (i in 1:g){
    E = E + (vec.E[[i]]*vec.ni[i]) # sum of scatter
  }
  vec.Xbar = lapply(dlist, colMeans)
  Xbar = rep(0,p)
  for (i in 1:g){
    Xbar = Xbar + (vec.ni[i]*vec.Xbar[[i]])
  }
  Xbar = Xbar/n
  H = array(0,c(p,p))
  for (i in 1:g){
    bdf = as.vector(vec.Xbar[[i]]-Xbar) # bars' difference
    H = H + outer(bdf,bdf)*vec.ni[i]
  }
  e = (n-g)
  h = (g-1)
  
  ##############################################################
  # COMPUTATION : MAIN PART
  # 1. the statistic
  Tnp = (1/sqrt(n-1))*((aux_trace(H)/h) - (aux_trace(E)/e))
  # 2. variance term
  a = (1/((e+2)*(e-1)))*(aux_trace(E%*%E) - (1/e)*(aux_trace(E)^2))
  Tnp.var = (2/h)*(a/e)
  # 3. adjusted statistic and p-value
  thestat = Tnp
  thestat.adj = Tnp/sqrt(Tnp.var)
  pvalue  = stats::pnorm(thestat.adj, lower.tail = FALSE)
  
  
  ##############################################################
  # REPORT
  hname   = "Test for Equality of Means by Schott (2007)"
  DNAME = deparse(substitute(dlist))
  Ha    = "one of equalities does not hold."
  names(thestat) = "Tnp"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}