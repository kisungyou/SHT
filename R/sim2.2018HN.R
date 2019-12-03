#' Two-sample Simultaneous Test of Means and Covariances by Hyodo and Nishiyama (2018)
#' 
#' Given a multivariate sample \eqn{X}, hypothesized mean \eqn{\mu_0} and covariance \eqn{\Sigma_0}, it tests
#' \deqn{H_0 : \mu_x = \mu_y \textrm{ and } \Sigma_x = \Sigma_y \quad vs\quad H_1 : \textrm{ not } H_0}
#' using the procedure by Hyodo and Nishiyama (2018) in a similar fashion to that of Liu et al. (2017) for 
#' one-sample test.
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
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
#' sim2.2018HN(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(121*10), ncol=10)
#'   Y = matrix(rnorm(169*10), ncol=10)
#'   counter[i] = ifelse(sim2.2018HN(X,Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'sim2.2018HN'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{hyodo_simultaneous_2018}{SHT}
#' 
#' @export
sim2.2018HN <- function(X, Y){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* sim2.2018HN : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # PARAMETERS 
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  
  ##############################################################
  # COMPUTATION 1 : PRELIMINARY
  S1 = cov(X) # Sg
  S2 = cov(Y)
  
  xbar = as.vector(colMeans(X))
  ybar = as.vector(colMeans(Y))
  
  K1 = 0
  K2 = 0
  for (i in 1:n1){
    vecdiff = as.vector(X[i,])-xbar
    K1 = K1 + (sum(vecdiff^2)^2)
  }
  for (j in 1:n2){
    vecdiff = as.vector(Y[j,])-ybar
    K2 = K2 + (sum(vecdiff^2)^2)
  }
  K1 = K1/(n1-1)
  K2 = K2/(n2-1)
  
  ##############################################################
  # COMPUTATION 2 : UNBIASED ESTIMATORS
  trS12  = aux_trace(S1%*%S2)
  Sigsq1 = ((n1-1)/(n1*(n1-2)*(n1-3)))*((n1-1)*(n1-2)*aux_trace(S1%*%S1) + (aux_trace(S1)^2) - (n1*K1))
  Sigsq2 = ((n2-1)/(n2*(n2-2)*(n2-3)))*((n2-1)*(n2-2)*aux_trace(S2%*%S2) + (aux_trace(S2)^2) - (n2*K2))
  
  est_delsq = sum((xbar-ybar)^2) - aux_trace(S1)/n1 - aux_trace(S2)/n2 # different of means
  est_Delsq = Sigsq1 + Sigsq2 - 2*trS12
  
  ##############################################################
  # COMPUTATION 3 : LEADING VARIANCES
  sig10sq = 2*Sigsq1/(n1^2) + 2*Sigsq2/(n2^2) + 4*trS12/(n1*n2)
  sig20sq = 4*(Sigsq1^2)/(n1^2) + 4*(Sigsq2^2)/(n2^2) + 8*(trS12^2)/(n1*n2)
  
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = (est_delsq/sqrt(sig10sq)) + (est_Delsq/sqrt(sig20sq))
  pvalue  = pnorm(thestat/sqrt(2),lower.tail=FALSE) # reject if (T/sqrt(2) > thr_alpha)
  
  hname   = "Two-sample Simultaneous Test of Means and Covariances by Hyodo and Nishiyama (2018)"
  Ha      = "both means and covariances are not equal."
  
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "T"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
  
}