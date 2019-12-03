#' Multivariate Test of Uniformity based on Normal Quantiles by Yang and Modarres (2017) 
#' 
#' Given a multivariate sample \eqn{X}, it tests
#' \deqn{H_0 : \Sigma_x = \textrm{ uniform on } \otimes_{i=1}^p [a_i,b_i] \quad vs\quad H_1 : \textrm{ not } H_0}
#' using the procedure by Yang and Modarres (2017). Originally, it tests the goodness of fit 
#' on the unit hypercube \eqn{[0,1]^p} and modified for arbitrary rectangular domain. Since 
#' this method depends on quantile information, every observation should strictly reside within 
#' the boundary so that it becomes valid after transformation.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param lower length-\eqn{p} vector of lower bounds of the test domain.
#' @param upper length-\eqn{p} vector of upper bounds of the test domain.
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
#' smallX = matrix(runif(10*3),ncol=3)
#' unif.2017YMq(smallX) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1234
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(runif(50*5), ncol=25)
#'   counter[i] = ifelse(unif.2017YMq(X)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'unif.2017YMq'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{yang_multivariate_2017}{SHT}
#' 
#' @export
unif.2017YMq <- function(X, lower=rep(0,ncol(X)), upper=rep(1,ncol(X))){
  ##############################################################
  # PREPROCESSING
  check_nd(X)    # univariate vector
  n = nrow(X)
  d = ncol(X)
  if (d<2){
    stop("* unif.2017YMq : use univariate functions for univariate data.")
  }
  if ((any(lower!=rep(0,d)))||(any(upper!=rep(1,d)))){
    X = aux_adjustcube(X,lower,upper)  
  }
  if ((any(X<0))||(any(X>1))){
    stop("* unif.2017YMq : since it utilizes quantile function, elements in X should be in [0,1] even after transformation.")
  }
  if (length(lower)!=d){
    stop("* unif.2017YMq : 'lower' should be a vector of length 'ncol(X)'.")
  }
  if (length(upper)!=d){
    stop("* unif.2017YMq : 'upper' should be a vector of length 'ncol(X)'.")
  }
  
  ##############################################################
  # COMPUTATION
  Z    = qnorm(X)
  Zbar = colMeans(Z)
  Cn   = n*sum(Zbar^2)
  
  thestat = Cn
  pvalue  = pchisq(thestat, d, lower.tail = FALSE)
  names(thestat) = "Cn"
  
  ##############################################################
  # REPORT
  hname = "Multivariate Test of Uniformity based on Normal Quantiles by Yang and Modarres (2017)"
  DNAME = deparse(substitute(X))
  Ha    = paste("Sample ", DNAME, " does not follow uniform distribution.",sep="")
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
