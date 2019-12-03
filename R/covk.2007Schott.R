#' Test for Homogeneity of Covariances by Schott (2007)
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \Sigma_1 = \cdots \Sigma_k\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Schott (2007).
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
#' covk.2007Schott(tinylist) # run the test
#' 
#' \donttest{
#' ## test when k=4 samples with (n,p) = (100,20)
#' ## empirical Type 1 error 
#' niter   = 1234
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   mylist = list()
#'   for (j in 1:4){
#'      mylist[[j]] = matrix(rnorm(100*20),ncol=20)
#'   }
#'   
#'   counter[i] = ifelse(covk.2007Schott(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'covk.2007Schott'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{schott_test_2007}{SHT}
#' 
#' @export
covk.2007Schott <- function(dlist){
  ##############################################################
  # PREPROCESSING
  check_dlistnd(dlist) 
  
  ##############################################################
  # PREPARATION
  g     = length(dlist)  # g-sample case
  p     = ncol(dlist[[1]])
  vec.n = unlist(lapply(dlist, nrow)) - 1 ## nS ~ W(\Sigma, n) means nS = scatter, S = (1/n)*scatter ## ERRATA it gotta be.
  vec.S = array(0,c(p,p,g))
  for (i in 1:g){
    vec.S[,,i] = stats::cov(dlist[[i]])
  }
  
  n = sum(vec.n)
  S = array(0,c(p,p))
  for (i in 1:g){
    S = S + ((vec.n[i]/n)*vec.S[,,i])
  }
  
  # tr(Si) and tr(Si^2)
  vec.trS  = rep(0,g)
  vec.trS2 = rep(0,g)
  for (i in 1:g){
    Si = vec.S[,,i]
    vec.trS[i]  = sum(diag(Si))
    vec.trS2[i] = sum(diag(Si%*%Si))
  }
  
  ##############################################################
  # COMPUTATION 1 : tnm
  tnm = 0
  for (i in 1:g){
    ni = vec.n[i]
    ei = (ni+2)*(ni-1)
    Si = vec.S[,,i]
    for (j in 1:g){
      nj = vec.n[j]
      ej = (nj+2)*(nj-1)
      Sj = vec.S[,,j]
      
      if (i<j){
        add1 = (1 - (ni-2)/ei)*vec.trS2[i]
        add2 = (1 - (nj-2)/ej)*vec.trS2[j]
        min1 = 2*sum(diag(Si%*%Sj))
        min2 = (ni/ei)*((vec.trS[i])^2)
        min3 = (nj/ej)*((vec.trS[j])^2)
        
        tnm = tnm + (add1+add2) - (min1+min2+min3)
      }
    }
  }
  
  ##############################################################
  # COMPUTATION 2 : variance of tnm
  a = ((n^2)/((n+2)*(n-1)))*(sum(diag(S%*%S)) - (1/n)*(sum(diag(S))^2))
  
  inner1 = 0
  for (i in 1:g){
    ni = vec.n[i]
    for (j in 1:g){
      nj = vec.n[j]
      if (i<j){
        inner1 = inner1 + (((ni+nj)/(ni*nj))^2)
      }
    }
  }
  inner2 = (g-1)*(g-2)*sum(1/(vec.n^2))
  theta  = 2.0*sqrt(inner1+inner2)*a
  
  thestat = tnm/theta
  pvalue  = pnorm(tnm/theta, lower.tail = FALSE) 

  ##############################################################
  # FINALE
  hname   = "Test for Homogeneity of Covariances by Schott (2007)"
  Ha      = "at least one of equalities does not hold."
  
  DNAME = deparse(substitute(dlist)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}