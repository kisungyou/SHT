#' One-sample Simultaneous Likelihood Ratio Test of Mean and Covariance
#' 
#' Given a multivariate sample \eqn{X}, hypothesized mean \eqn{\mu_0} and covariance \eqn{\Sigma_0}, it tests
#' \deqn{H_0 : \mu_x = \mu_0 \textrm{ and } \Sigma_x = \Sigma_0 \quad vs\quad H_1 : \textrm{ not } H_0}
#' using the standard likelihood-ratio test procedure.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param mu0 a length-\eqn{p} mean vector of interest.
#' @param Sigma0 a \eqn{(p\times p)} given covariance matrix.
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
#' sim1.LRT(smallX) # run the test
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(100*10), ncol=10)
#'   counter[i] = ifelse(sim1.LRT(X)$p.value < 0.05, 1, 0)
#'   print(paste("* iteration ",i,"/1000 complete..."))
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'sim1.LRT'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @concept simtest
#' @export
sim1.LRT <- function(X, mu0=rep(0,ncol(X)),  Sigma0=diag(ncol(X))){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  n = nrow(X)
  p = ncol(X)

  ##############################################################
  # CENTER AND SCALE PROPERLY
  minusX = matrix(rep(mu0,n), ncol=p, byrow=TRUE)
  scaler = aux_getinvroot(Sigma0)
  Xnew   = (X-minusX)%*%scaler
  
  ##############################################################
  # COMPUTATION
  xbar = as.vector(colMeans(Xnew))
  S    = stats::cov(Xnew)*(n-1)/n
  LR   = aux_trace(S) - base::log(base::det(S)) - p + base::sum(xbar*xbar)
  
  thestat = n*LR
  mydf    = p*(1+p)/2 + p
  pvalue  = stats::pchisq(thestat, df=mydf, lower.tail = FALSE)
  
  ##############################################################
  # REPORT
  hname   = "One-sample Simultaneous Likelihood Ratio Test of Mean and Covariance"
  Ha      = "both mean and covariance are not equal to mu0 and Sigma0."
  
  DNAME = deparse(substitute(X))
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}

