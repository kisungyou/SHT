#' Two-Sample Covariance Test with Maximum Pairwise Bayes Factor
#' 
#' 
#' @export
cov2.mxPBF <- function(X, Y, a0=2.0, b0=2.0, gamma=1.0, nCores=2){
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
  
  ##############################################################
  # PRELIMINARY
  Xnew = as.matrix(scale(X, center=TRUE, scale=FALSE))
  Ynew = as.matrix(scale(Y, center=TRUE, scale=FALSE))
  Znew = rbind(Xnew, Ynew)
  
  ##############################################################
  # MAIN COMPUTATION
  log.BF.mat = cpp_cov2_mxPBF(Xnew, Ynew, Znew, a0, b0, gamma, nCores)
  diag(log.BF.mat) = -Inf
  
  ##############################################################
  # FINALE
  hname   = "Two-Sample Covariance Test with Maximum Pairwise Bayes Factor"
  Ha      = "two covariances are not equal."
  
  thestat = max(exp(log.BF.mat))
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "BF"
  res   = list(statistic=thestat, alternative = Ha, method=hname, data.name = DNAME, log.BF.mat = log.BF.mat)
  class(res) = "htest"
  return(res)
}

# 
# count = 0
# for (i in 1:1000){
#   X = matrix(rnorm(50*100), ncol=100)
#   Y = matrix(rnorm(30*100), ncol=100)
#   if (max(exp(cov2.mxPBF(X,Y)$log.BF.mat))>10){
#     count = count + 1
#   }
# }
# 
# X = matrix(rnorm(50*100), ncol=100)
# Y = matrix(rnorm(30*100), ncol=100)
# aa = mxPBF::testcov2(X, Y)
# bb = cov2.mxPBF(X,Y)