#' One-sample Covariance Test with Maximum Pairwise Bayes Factor
#' 
#' Not Written Here - No Reference Yet.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param a0 shape parameter for inverse-gamma prior.
#' @param b0 scale parameter for inverse-gamma prior.
#' @param gamma non-negative variance scaling parameter.
#' @param nthreads number of threads for parallel execution via OpenMP.
#'
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{maximum of pairwise Bayes factor.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' \item{log.BF.mat}{matrix of pairwise Bayes factors in natural log.}
#' }
#'
#' @examples
#' \dontrun{
#' ## generate from multivariate normal with identity covariance.
#' p = 10 # dimensionality
#' data = matrix(rnorm(100*p), ncol=p)
#'
#' ## run test with different parameters
#' out1 = cov1.mxPBF(data)
#' out2 = cov1.mxPBF(data, a0=5.0, b0=5.0) # change some params
#'
#' ## visualize two Bayes Factor matrices
#' par(mfrow=c(1,2), pty="s")
#' image(exp(out1$log.BF.mat)[,p:1], main="a0=b0=1.0")
#' image(exp(out2$log.BF.mat)[,p:1], main="a0=b0=5.0")
#' }
#' 
#' @keywords internal
#' @noRd
cov1.mxPBF <- function(X, Sigma0=diag(ncol(X)), a0=2.0, b0=2.0, gamma=1.0, nthreads=1){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  n = nrow(X)
  p = ncol(X)
  if ((nrow(X)==1)||(p<2)){
    stop("* cov1.mxPBF : input 'X' is invalid. Provide a multivariate sample of multiple observations.")
  }
  if ((length(a0)!=1)||(a0<=0)){
    stop("* cov1.mxPBF : 'a0' should be a nonnegative number.")
  }
  if ((length(b0)!=1)||(b0<=0)){
    stop("* cov1.mxPBF : 'b0' should be a nonnegative number.")
  }
  if ((length(gamma)!=1)||(gamma<=0)){
    stop("* cov1.mxPBF : 'gamma' should be a nonnegative number.")
  }
  nCores = as.integer(nthreads)
  if (nCores < 1){
    stop("* cov1.mxPBF : 'nthreads' should be a positive integer.")
  }
  
  ##############################################################
  # PREPROCESSING : ADJUST THE DATA FOR TESTING
  scaler = aux_getinvroot(Sigma0)
  X.centered = scale(X, center=TRUE, scale=FALSE)
  X.adjusted = (matrix(X.centered,nrow=n) %*% scaler)
  
  ##############################################################
  # MAIN COMPUTATION
  if (nCores==1){
    log.BF.mat = cpp_cov1_mxPBF_single(X.adjusted, a0, b0, gamma)
  } else {
    log.BF.mat = cpp_cov1_mxPBF_multiple(X.adjusted , a0, b0, gamma, nCores)
  }
  diag(log.BF.mat) = -Inf
  
  ##############################################################
  # FINALE
  hname   = "One-sample Covariance Test with Maximum Pairwise Bayes Factor"
  Ha      = "true covariance is different from Sigma0."
  
  thestat = max(exp(log.BF.mat))
  DNAME = deparse(substitute(X)) # borrowed from HDtest
  names(thestat) = "maximum BF"
  res   = list(statistic=thestat, alternative = Ha, method=hname, data.name = DNAME, log.BF.mat = log.BF.mat)
  class(res) = "htest"
  return(res)
}


# DEVELOPMENTAL TEST CODE
# p = 50
# X = matrix(rnorm(100*p), ncol=p)
# aa = mxPBF::testcov1(X)
# bb = cov1.mxPBF(X)
# c(max(aa$log.BF.mat), max(bb$log.BF.mat))
# 
# b1 = cov1.mxPBF(X,nthreads=1)
# b2 = cov1.mxPBF(X,nthreads=2)
# b4 = cov1.mxPBF(X,nthreads=4)
# b8 = cov1.mxPBF(X,nthreads=8)
# c(max(aa$log.BF.mat), max(b1$log.BF.mat), max(b2$log.BF.mat), max(b4$log.BF.mat), max(b8$log.BF.mat))
# 
# library(microbenchmark)
# timeout = microbenchmark(list=alist(R=mxPBF::testcov1(X),
#                                     C1=cov1.mxPBF(X,nthreads=1),
#                                     C2=cov1.mxPBF(X,nthreads=2),
#                                     C4=cov1.mxPBF(X,nthreads=4),
#                                     C8=cov1.mxPBF(X,nthreads=8), times=10))