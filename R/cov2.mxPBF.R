#' Two-sample Covariance Test with Maximum Pairwise Bayes Factor
#' 
#' Not Written Here - No Reference Yet.
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
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
#' ## empirical Type 1 error with BF threshold = 20
#' niter   = 12345
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(cov2.mxPBF(X,Y)$statistic > 20, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov2.mxPBF'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @export
cov2.mxPBF <- function(X, Y, a0=2.0, b0=2.0, gamma=1.0, nthreads=1){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* cov2.mxPBF : two samples X and Y should be of same dimension.")
  }
  p = ncol(X)
  if ((nrow(X)==1)||(nrow(Y)==1)||(p<2)){
    stop("* cov2.mxPBF : inputs are invalid. Provide multivariate samples with multiple observations.")
  }
  if ((length(a0)!=1)||(a0<=0)){
    stop("* cov2.mxPBF : 'a0' should be a nonnegative number.")
  }
  if ((length(b0)!=1)||(b0<=0)){
    stop("* cov2.mxPBF : 'b0' should be a nonnegative number.")
  }
  if ((length(gamma)!=1)||(gamma<=0)){
    stop("* cov2.mxPBF : 'gamma' should be a nonnegative number.")
  }
  nCores = as.integer(nthreads)
  if (nCores < 1){
    stop("* cov2.mxPBF : 'nthreads' should be a positive integer.")
  }
  
  ##############################################################
  # PRELIMINARY
  Xnew = as.matrix(scale(X, center=TRUE, scale=FALSE))
  Ynew = as.matrix(scale(Y, center=TRUE, scale=FALSE))

  ##############################################################
  # MAIN COMPUTATION
  if (nCores==1){
  log.BF.mat = cpp_cov2_mxPBF_single(Xnew, Ynew, a0, b0, gamma)
  } else {
    log.BF.mat = cpp_cov2_mxPBF_multiple(Xnew, Ynew, a0, b0, gamma, nCores)
  }
  diag(log.BF.mat) = -Inf
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Covariance Test with Maximum Pairwise Bayes Factor"
  Ha      = "two covariances are not equal."
  
  thestat = max(exp(log.BF.mat))
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "maximum BF"
  res   = list(statistic=thestat, alternative = Ha, method=hname, data.name = DNAME, log.BF.mat = log.BF.mat)
  class(res) = "htest"
  return(res)
}

# 
# count = 0
# for (i in 1:1000){
#   X = matrix(rnorm(50*100), ncol=100)
#   Y = matrix(rnorm(30*100), ncol=100)
  # if (max(exp(cov2.mxPBF(X,Y)$log.BF.mat))>10){
#     count = count + 1
#   }
# }
# # 
# p = 50
# X = matrix(rnorm(100*p), ncol=p)
# Y = matrix(rnorm(80*p), ncol=p)
# aa = mxPBF::testcov2(X, Y)
# b1 = cov2.mxPBF(X,Y,nthreads=1)
# b2 = cov2.mxPBF(X,Y,nthreads=2)
# b4 = cov2.mxPBF(X,Y,nthreads=4)
# b8 = cov2.mxPBF(X,Y,nthreads=8)
# c(max(aa$log.BF.mat), max(b1$log.BF.mat), max(b2$log.BF.mat), max(b4$log.BF.mat), max(b8$log.BF.mat))
# 
# library(microbenchmark)
# timeout = microbenchmark(list=alist(R=mxPBF::testcov2(X,Y),
#                                     C1=cov2.mxPBF(X,Y,nCores=1),
#                                     C2=cov2.mxPBF(X,Y,nCores=2),
#                                     C4=cov2.mxPBF(X,Y,nCores=4),
#                                     C8=cov2.mxPBF(X,Y,nCores=8), times=10))