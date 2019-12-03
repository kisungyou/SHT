#' One-sample Test for Covariance Matrix by Fisher (2012)
#' 
#' Given a multivariate sample \eqn{X} and hypothesized covariance matrix \eqn{\Sigma_0}, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_0\quad vs\quad H_1 : \Sigma_x \neq \Sigma_0}
#' using the procedure by Fisher (2012). This method utilizes the generalized form of the inequality
#' \deqn{\frac{1}{p} \sum_{i=1}^p (\lambda_i^r - 1)^{2s} \ge 0} and offers two types of 
#' test statistics \eqn{T_1} and \eqn{T_2} corresponding to the case \eqn{(r,s)=(1,2)} and \eqn{(2,1)} respectively.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
#' @param type \code{1} or \code{2} for corresponding statistic from the paper. 
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
#' cov1.2012Fisher(smallX) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter1 = rep(0,niter)  # p-values of the type 1
#' counter2 = rep(0,niter)  # p-values of the type 2
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=50) # (n,p) = (5,50)
#'   counter1[i] = ifelse(cov1.2012Fisher(X, type=1)$p.value < 0.05, 1, 0)
#'   counter2[i] = ifelse(cov1.2012Fisher(X, type=2)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov1.2012Fisher' \n","*\n",
#' "* empirical error with statistic 1 : ", round(sum(counter1/niter),5),"\n",
#' "* empirical error with statistic 2 : ", round(sum(counter2/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{fisher_testing_2012}{SHT}
#' 
#' @export
cov1.2012Fisher <- function(X, Sigma0=diag(ncol(X)), type){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  n = nrow(X)-1 # be careful about this notation !
  p = ncol(X)
  c = (p/n)
  S = cov(X)
  if (missing(type)){
    type = 1
  }
  mytype = as.character(as.integer(type))
  
  ##############################################################
  # CENTER AND SCALE PROPERLY
  X.centered  = as.matrix(scale(X, center=TRUE, scale=FALSE))
  scaler      = aux_getinvroot(Sigma0)
  X.processed = X.centered%*%scaler
  
  ##############################################################
  # MAIN COMPUTATION
  vec.a = cov1.2012Fisher.estimators(n,p,S)
  a1 = vec.a[1]
  a2 = vec.a[2]
  a3 = vec.a[3]
  a4 = vec.a[4]
  
  if (pracma::strcmp(mytype,"1")){
    thestat = (n/(c*sqrt(8)))*(a4-(4*a3)+(6*a2)-(4*a1)+1)
  } else if (pracma::strcmp(mytype,"2")){
    thestat = (n/sqrt(8*((c^2)+(12*c)+8)))*(a4-(2*a2)+1)
  } else {
    stop("* cov1.2012Fisher : 'type' should be either an integer 1 or 2.")
  }
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "One-sample Test for Covariance Matrix by Fisher (2012)."
  Ha      = "true covariance is different from Sigma0."
  DNAME = deparse(substitute(X)) # borrowed from HDtest
  names(thestat) = "T"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}



# auxiliary variables -----------------------------------------------------
#' @keywords internal
#' @noRd
cov1.2012Fisher.estimators <- function(n,p,S){
  S2  = S%*%S
  S3  = S2%*%S
  S4  = S3%*%S
  trS  = aux_trace(S)
  trS2 = aux_trace(S2)
  trS3 = aux_trace(S3)
  trS4 = aux_trace(S4)

  gamma.top = (n^5)*((n^2)+n+2)
  gamma.bot = (n+1)*(n+2)*(n+4)*(n+6)*(n-1)*(n-2)*(n-3)
  gamma = (gamma.top/gamma.bot)
  
  tau = ((n^4)/((n-1)*(n-2)*(n+2)*(n+4)))
    
  a1 = (1/p)*trS
  a2 = ((n^2)/((n-1)*(n+2)*p))*(trS2 - (1/n)*(trS^2))
  a3 = (tau/p)*(trS3 - (3/n)*trS2*trS + (2/(n^2))*(trS^3))
  a4 = (gamma/p)*(trS4 - (4/n)*trS3*trS - ((2*(n^2)+3*n-6)/(n*(n^2 +n+2)))*(trS2^2) + 2*((5*n)+6)*trS2*(trS^2)/(n*((n^2)+n+2)) - ((5*n)+6)*(trS^4)/(n^4+n^3+(2*(n^2))))

  return(c(a1,a2,a3,a4))
}