#' Test for Homogeneity of Covariances by Schott (2001)
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \Sigma_1 = \cdots \Sigma_k\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Schott (2001) using Wald statistics. In the original paper, it provides 4 
#' different test statistics for general elliptical distribution cases. However, we only deliver 
#' the first one with an assumption of multivariate normal population.
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
#' covk.2001Schott(tinylist) # run the test
#' 
#' \dontrun{
#' ## test when k=5 samples with (n,p) = (100,20)
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   mylist = list()
#'   for (j in 1:5){
#'      mylist[[j]] = matrix(rnorm(100*20),ncol=20)
#'   }
#'   
#'   counter[i] = ifelse(covk.2001Schott(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'covk.2001Schott'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{schott_tests_2001}{SHT}
#' 
#' @export
covk.2001Schott <- function(dlist){
  ##############################################################
  # PREPROCESSING
  check_dlistnd(dlist) 
  
  ##############################################################
  # PREPARATION
  g     = length(dlist)  # g-sample case
  p     = ncol(dlist[[1]])
  vec.n = unlist(lapply(dlist, nrow))-1
  vec.S = array(0,c(p,p,g)) 
  for (i in 1:g){
    vec.S[,,i] = stats::cov(dlist[[i]])
  }
  
  n = sum(vec.n)
  S = array(0,c(p,p))
  for (i in 1:g){
    S = S + ((vec.n[i]/n)*vec.S[,,i])
  }
  vec.gamma = vec.n/n
  
  # Sinv = pracma::pinv(S)  # use pinv for Sinv
  Sinv = tryCatch({solve(S)},error = function(e){pracma::pinv(S)})
  vec.Sinv2 = rep(0,g) # tr((SSinv)^2)
  for (i in 1:g){
    SSinv = vec.S[,,i]%*%Sinv
    vec.Sinv2[i] = sum(diag(SSinv%*%SSinv))
  }
  
  ##############################################################
  # MAIN COMPUTATION
  # 1. first term
  term1 = 0
  for (i in 1:g){
    term1 = term1 + vec.gamma[i]*vec.Sinv2[i]
  }
  
  # 2. second term
  term2 = 0
  for (i in 1:g){
    gi = vec.gamma[i]
    Si = vec.S[,,i]
    for (j in 1:g){
      gj = vec.gamma[j]
      Sj = vec.S[,,j]
      
      term2 = term2 + (gi*gj)*sum(diag(Si%*%Sinv%*%Sj%*%Sinv))
    }
  }
  
  # wrap up
  thestat = (term1-term2)*n/2
  thedf   = (g-1)*p*(p+1)/2
  pvalue  = stats::pchisq(thestat, df=thedf, lower.tail = FALSE)
  
  ##############################################################
  # FINALE
  hname   = "Test for Homogeneity of Covariances by Schott (2001)"
  Ha      = "at least one of equalities does not hold."
  
  DNAME = deparse(substitute(dlist)) # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}