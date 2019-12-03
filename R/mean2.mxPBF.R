#' Two-sample Mean Test with Maximum Pairwise Bayes Factor
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
#' \item{log.BF.vec}{vector of pairwise Bayes factors in natural log.}
#' }
#' 
#' @examples 
#' \dontrun{
#' ## empirical Type 1 error with BF threshold = 10
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(100*10), ncol=10)
#'   Y = matrix(rnorm(200*10), ncol=10)
#'   
#'   counter[i] = ifelse(mean2.mxPBF(X,Y)$statistic > 10, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.mxPBF'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#'
#' 
#' @export
mean2.mxPBF <- function(X, Y, a0=2.0, b0=2.0, gamma=1.0, nthreads=1){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.mxPBF : two samples X and Y should be of same dimension.")
  }
  p = ncol(X)
  if ((nrow(X)<2)||(nrow(Y)<2)||(p<2)){
    stop("* mean2.mxPBF : inputs are invalid. Provide multivariate samples with multiple observations.")
  }
  if ((length(a0)>1)||(a0<=0)){
    stop("* mean2.mxPBF : 'a0' should be a nonnegative number.")
  }
  if ((length(b0)>1)||(b0<=0)){
    stop("* mean2.mxPBF : 'b0' should be a nonnegative number.")
  }
  if ((length(gamma)>1)||(gamma<=0)){
    stop("* mean2.mxPBF : 'gamma' should be a nonnegative number.")
  }
  nCores = as.integer(nthreads)
  if (nCores < 1){
    stop("* mean2.mxPBF : 'nthreads' should be a positive integer.")
  }
  
  ##############################################################
  # MAIN COMPUTATION
  if (nCores==1){
    log.BF.vec = as.vector(cpp_mean2_mxPBF_single(X, Y, a0, b0, gamma))  
  } else {
    log.BF.vec = as.vector(cpp_mean2_mxPBF_multiple(X, Y, a0, b0, gamma, nCores))
  }
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Mean Test with Maximum Pairwise Bayes Factor"
  Ha      = "two means are not equal."
  
  thestat = max(exp(log.BF.vec))
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "maximum BF"
  res   = list(statistic=thestat, alternative = Ha, method=hname, data.name = DNAME, log.BF.vec = log.BF.vec)
  class(res) = "htest"
  return(res)
}
# 
# p = 100
# X = matrix(rnorm(100*p), ncol=p)
# Y = matrix(rnorm(80*p), ncol=p)
# aa = mxPBF::testmean2(X, Y)
# bb = mean2.mxPBF(X,Y)
# 
# b1 = mean2.mxPBF(X,Y,nthreads=1)
# b2 = mean2.mxPBF(X,Y,nthreads=2)
# b4 = mean2.mxPBF(X,Y,nthreads=4)
# b8 = mean2.mxPBF(X,Y,nthreads=8)
# c(max(aa$log.BF.vec), max(b1$log.BF.vec), max(b2$log.BF.vec), max(b4$log.BF.vec), max(b8$log.BF.vec))
# 
# library(microbenchmark)
# timeout = microbenchmark(list=alist(R=mxPBF::testmean2(X,Y),
#                                     C1=mean2.mxPBF(X,Y,nthreads=1),
#                                     C2=mean2.mxPBF(X,Y,nthreads=2),
#                                     C4=mean2.mxPBF(X,Y,nthreads=4),
#                                     C8=mean2.mxPBF(X,Y,nthreads=8), times=10))