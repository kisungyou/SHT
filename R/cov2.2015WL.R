#' Two-sample Test for Covariance Matrices by Wu and Li (2015)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_y\quad vs\quad H_1 : \Sigma_x \neq \Sigma_y}
#' using the procedure by Wu and Li (2015).
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param m the number of random projections to be applied.
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
#' smallX = matrix(rnorm(10*3),ncol=3)
#' smallY = matrix(rnorm(10*3),ncol=3)
#' cov2.2015WL(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(cov2.2015WL(X, Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov2.2015WL'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{wu_tests_2015}{SHT}
#' 
#' @export
cov2.2015WL <- function(X, Y, m=50){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* cov2.2015WL : two samples X and Y should be of same dimension.")
  }
  m = as.integer(m)
  
  ##############################################################
  # PARAMETERS and CENTERING
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  
  Xnew = as.matrix(scale(X, center=TRUE, scale=FALSE))
  Ynew = as.matrix(scale(Y, center=TRUE, scale=FALSE))
  
  ##############################################################
  # LET'S RUN MULTIPLE ITERATIONS 
  rec.stat = rep(0,m)
  for (i in 1:m){
    projvec = rnorm(p)
    projvec = projvec/sqrt(sum(projvec*projvec))
    
    Xproj = as.vector(Xnew%*%projvec) # projection onto 1-dimensional space
    Yproj = as.vector(Ynew%*%projvec)
    
    s1 = sum(Xproj^2)/n1
    s2 = sum(Yproj^2)/n2
    rec.stat[i] = (((2/n1)+(2/n2))^(-1/2))*log(s1/s2)
  }
  thestat = max(rec.stat)
  pvalue  = 1-(pnorm(thestat, lower.tail=TRUE)^m)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  hname   = "Two-sample Test for Covariance Matrices by Wu and Li (2015)"
  Ha      = "two covariances are not equal."
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "T2m"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}

