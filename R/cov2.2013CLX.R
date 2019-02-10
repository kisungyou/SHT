#' Two-Sample Test for High-Dimensional Covariances by Cai, Liu, and Xia (2013)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_y\quad vs\quad H_1 : \Sigma_x \neq \Sigma_y}
#' using the procedure by Cai, Liu, and Xia (2013).
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
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
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(cov2.2013CLX(X, Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov2.2013CLX'\n\n",
#' sprintf("* number of rejections   : %d\n",sum(counter)),
#' sprintf("* total number of trials : %d\n",niter),
#' sprintf("* empirical Type 1 error : %.4f\n", sum(counter/niter)),sep=""))
#' }
#' 
#' @references 
#' \insertRef{cai_two-sample_2013}{SHT}
#' 
#' 
#' @export
cov2.2013CLX <- function(X, Y, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  check_alpha(alpha)
  if (ncol(X)!=ncol(Y)){
    stop("* cov2.2013CLX : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # BORROWED FROM JAY
  # parameter setting
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  # elementary computation : sample covariance with new df
  Sigma1.hat = cov(X)*(n1-1)/n1
  Sigma2.hat = cov(Y)*(n2-1)/n2
  # mean estimation
  bar.X = colMeans(X)
  theta1.hat = matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      theta1.hat[i,j] = sum( ((X[,i]-bar.X[i])*(X[,j]-bar.X[j]) - Sigma1.hat[i,j])^2 )/n1
    }
  }
  bar.Y = colMeans(Y)
  theta2.hat = matrix(0, nrow=p, ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      theta2.hat[i,j] = sum( ((Y[,i]-bar.Y[i])*(Y[,j]-bar.Y[j]) - Sigma2.hat[i,j])^2 )/n2
    }
  }
  # statistic and results
  M.mat = (Sigma1.hat - Sigma2.hat)^2 / (theta1.hat/n1 + theta2.hat/n2)
  Mn = max(M.mat) # test statistic
  pval.num = Mn - 4*log(p) + log(log(p))
  pvalue  = 1-exp(-(1/sqrt(8*pi))*exp(-pval.num/2))
  
  
  
  ##############################################################
  # FINALE
  hname   = "Two-Sample Test for High-Dimensional Covariances by Cai, Liu, and Xia (2013)"
  Ha      = "two covariances are not equal."
  thestat = Mn
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
