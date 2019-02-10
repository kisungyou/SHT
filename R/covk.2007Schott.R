#' Test for Homogeneity of High-Dimensional Covariances by Schott (2007).
#' 
#' Given multivariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \Sigma_1 = \cdots \Sigma_k\quad vs\quad H_1 : \textrm{ not } H_0}
#' by a procedure by Schott (2007).
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample matrix of dimension \eqn{p}.
#' @param alpha significance level.
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
#'      mylist[[j]] = matrix(rnorm(32*8), ncol=8)
#'   }
#'   
#'   counter[i] = ifelse(covk.2007Schott(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'covk.2007Schott'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{schott_test_2007}{SHT}
#' 
#' @export
covk.2007Schott <- function(dlist, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_dlistnd(dlist) 
  check_alpha(alpha)
  
  ##############################################################
  # PARAMETERS
  g = length(dlist) # number of groups
  if (g < 2){
    stop("* covk.2007Schott : we need at least 2 sets of data.")
  }
  p     = ncol(dlist[[1]])    # dimension
  n.vec = unlist(dlist, nrow) # per-class observations
  n     = sum(n.vec)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  #   1. sample covariance matrices as 3d array : use implicit knowledge
  S3d = array(0,c(p,p,g))
  for (i in 1:g){
    # tgt = as.matrix(scale(dlist[[i]], center=TRUE, scale=FALSE))
    # S3d[,,i] = stats::cov(dlist[[i]])*((n.vec[i]-1)/n.vec[i])
    S3d[,,i] = cov(dlist[[i]])*(n.vec[i]-1)/n.vec[i]
  }
  S   = array(0,c(p,p)) # pooled covariance
  for (i in 1:g){
    S = S + (n.vec[i]*S3d[,,i]/n)
  }
  #   2. eta
  eta.vec = (n.vec+2)*(n.vec-1)
  #   3. trace of Si^2 and Si
  tr2.vec = rep(0,g)
  tr1.vec = rep(0,g)
  for (i in 1:g){
    Si = S3d[,,i]
    tr2.vec[i] = sum(diag(Si%*%Si)) 
    tr1.vec[i] = sum(diag(Si))
  }
  
  ##############################################################
  # COMPUTATION : MAIN COMPUTATION
  #   1. tnm
  tnm = 0
  for (j in 2:g){
    Sj = S3d[,,j]
    for (i in 1:(j-1)){
      Si = S3d[,,i]
      
      term1 = (1.0-((n.vec[i]-2.0)/eta.vec[i]))*tr2.vec[i]
      term2 = (1.0-((n.vec[j]-2.0)/eta.vec[j]))*tr2.vec[j]
      term3 = 2*sum(diag(Si%*%Sj))
      term4 = (n.vec[i]/eta.vec[i])*((tr1.vec[i])^2)
      term5 = (n.vec[j]/eta.vec[j])*((tr1.vec[j])^2)
      
      tnm   = tnm + (term1+term2-term3-term4-term5)
    }
  }
  #   2. that2
  a = ((n^2)/((n+2)*(n-1)))*(sum(diag(S%*%S)) - (1/n)*(sum(diag(S))^2))
  
  nsums = 0
  for (j in 2:g){
    nj = n.vec[j]
    for (i in 1:(j-1)){
      ni = n.vec[i]
      nsums = nsums + ((1/ni + 1/nj)^2)
    }
  }
  that2 = 4*(nsums + (g-1)*(g-2)*sum(1/(n.vec^2)) )*(a^2)
  
  
  ##############################################################
  # FINALE
  hname   = "Test for Homogeneity of High-Dimensional Covariances by Schott (2007)"
  Ha      = "at least one of equalities does not hold."
  thestat = tnm/sqrt(that2)
  pvalue  = pnorm(thestat, lower.tail = FALSE) # right-hand side
  
  DNAME = deparse(substitute(dlist))
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}