#' Two-sample Test for High-Dimensional Means by Srivastava and Du (2008)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \mu_x = \mu_y\quad vs\quad H_1 : \mu_x \neq \mu_y}
#' using the procedure by Srivastava and Du (2008).
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
#' mean2.2008SD(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(mean2.2008SD(X,Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.2008SD'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{srivastava_test_2008}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.2008SD <- function(X, Y){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.2008SD : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n1 = nrow(X)
  n2 = nrow(Y)
  n  = (n1+n2-2) # from highmean package
  p  = ncol(X)
  
  x1 = colMeans(X)
  x2 = colMeans(Y)
  Sn = (((n1-1)*cov(X) + (n2-1)*cov(Y))/n)
  
  dsvec  = aux_adjustvec(diag(Sn))
  Dsinv  = diag(1/dsvec)
  Dshalf = diag(1/sqrt(dsvec))
  
  R      = stats::cov2cor(Sn)
  R2     = (R%*%R)
  trR2   = aux_trace(R2)
  cpn    = (1+(trR2/(p^(3/2))))
  
  xdiff = as.vector(x1-x2)
  term1 = ((((n1*n2)/(n1+n2))*aux_quadform(Dsinv,xdiff))-((n*p)/(n-2)))
  term2 = sqrt((2*(trR2-((p^2)/n)))*cpn)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "Two-sample Test for High-Dimensional Means by Srivastava and Du (2008)."
  Ha      = "true means are different."
  # if (pvalue < alpha){
  #   conclusion = "Reject Null Hypothesis."
  # } else {
  #   conclusion = "Not Reject Null Hypothesis."
  # }
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}