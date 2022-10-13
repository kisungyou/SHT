#' Two-sample Test for High-Dimensional Covariances by Li and Chen (2012)
#'
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_y\quad vs\quad H_1 : \Sigma_x \neq \Sigma_y}
#' using the procedure by Li and Chen (2012). 
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param use.unbiased a logical; \code{TRUE} to use up to 4th-order U-statistics as proposed in the paper, \code{FALSE} for faster run under an assumption that \eqn{\mu_h = 0} (default: \code{TRUE}).
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
#' smallX = matrix(rnorm(10*4),ncol=5)
#' smallY = matrix(rnorm(10*4),ncol=5)
#' cov2.2012LC(smallX, smallY) # run the test
#' 
#' \dontrun{
#' ## empirical Type 1 error : use 'biased' version for faster computation
#' niter   = 1000
#' counter = rep(0,niter)
#' for (i in 1:niter){
#'   X = matrix(rnorm(500*25), ncol=10)
#'   Y = matrix(rnorm(500*25), ncol=10)
#'   
#'   counter[i] = ifelse(cov2.2012LC(X,Y,use.unbiased=FALSE)$p.value  < 0.05,1,0)
#'   print(paste0("iteration ",i,"/1000 complete.."))
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'cov2.2012LC'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#'
#' @references 
#' \insertRef{li_two_2012}{SHT}
#' 
#' @concept covariance
#' @export
cov2.2012LC <- function(X, Y, use.unbiased=TRUE){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* cov2.2012LC : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # COMPUTATION
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  
  if (use.unbiased){ # unbiased / slower / full
    A1  = cov2_2012LC_A(X)
    A2  = cov2_2012LC_A(Y)
    C12 = cov2_2012LC_C(X, Y)
  } else {
    X1 = as.matrix(scale(X, center = TRUE, scale = FALSE))
    X2 = as.matrix(scale(Y, center = TRUE, scale = FALSE))
    
    A1  = cov2_2012LC_A_biased(X1) # elements for test statistics
    A2  = cov2_2012LC_A_biased(X2)
    C12 = cov2_2012LC_C_biased(X1, X2)
  }
  
  Tn1n2 = (A1 + A2 - 2*C12)  # test statistic
  shat  = 2*(A1/n2 + A2/n1)  # variance term
  pvalue = stats::pnorm((Tn1n2/sqrt(shat)), lower.tail = FALSE)
  
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Test for High-Dimensional Covariances by Li and Chen (2012)"
  Ha      = "two covariances are not equal."
  thestat = Tn1n2/sqrt(shat)
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}