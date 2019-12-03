#' One-sample Test for Covariance Matrix by Wu and Li (2015)
#' 
#' Given a multivariate sample \eqn{X} and hypothesized covariance matrix \eqn{\Sigma_0}, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_0\quad vs\quad H_1 : \Sigma_x \neq \Sigma_0}
#' using the procedure by Wu and Li (2015). They proposed to use \eqn{m} number of multiple random projections 
#' since only a single operation might attenuate the efficacy of the test.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
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
#' cov1.2015WL(smallX) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' ##   compare effects of m=5, 10, 50
#' niter = 1000
#' rec1  = rep(0,niter) # for m=5
#' rec2  = rep(0,niter) #     m=10
#' rec3  = rep(0,niter) #     m=50
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*10), ncol=50) # (n,p) = (10,50)
#'   rec1[i] = ifelse(cov1.2015WL(X, m=5)$p.value < 0.05, 1, 0)
#'   rec2[i] = ifelse(cov1.2015WL(X, m=10)$p.value < 0.05, 1, 0)
#'   rec3[i] = ifelse(cov1.2015WL(X, m=50)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov1.2015WL'\n","*\n",
#' "* Type 1 error with m=5   : ",round(sum(rec1/niter),5),"\n",
#' "* Type 1 error with m=10  : ",round(sum(rec2/niter),5),"\n",
#' "* Type 1 error with m=50  : ",round(sum(rec3/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{wu_tests_2015}{SHT}
#' 
#' @export
cov1.2015WL <- function(X, Sigma0=diag(ncol(X)), m=25){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  n = nrow(X)
  p = ncol(X)
  m = as.integer(m)
  
  ##############################################################
  # CENTER AND SCALE PROPERLY
  X.centered  = as.matrix(scale(X, center=TRUE, scale=FALSE))
  scaler      = aux_getinvroot(Sigma0)
  X.processed = X.centered%*%scaler
  
  ##############################################################
  # LET'S RUN MULTIPLE ITERATIONS 
  rec.stat = rep(0,m)
  for (i in 1:m){
    projvec = rnorm(p)
    projvec = projvec/sqrt(sum(projvec*projvec))
    Y       = as.vector(X.processed%*%projvec)
    rec.stat[i] = sqrt(2*sum(Y^2)) - sqrt((2*n)-1)
  }
  thestat = max(rec.stat)
  pvalue  = 1-(pnorm(thestat, lower.tail=TRUE)^m)

  ##############################################################
  # COMPUTATION : DETERMINATION
  hname   = "One-sample Test for Covariance Matrix by Wu and Li (2015)."
  Ha      = "true covariance is different from Sigma0."
  DNAME = deparse(substitute(X)) # borrowed from HDtest
  names(thestat) = "T1m"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}