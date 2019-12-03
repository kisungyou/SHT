#' Two-sample Test for High-Dimensional Covariances by Li and Chen (2012)
#'
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \Sigma_x = \Sigma_y\quad vs\quad H_1 : \Sigma_x \neq \Sigma_y}
#' using the procedure by Li and Chen (2012). In accordance with a proposal 
#' by authors, we offer an option to use biased estimator instead for faster computation.
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param unbiased a logical; \code{FALSE} to use biased estimator with faster speed, \code{TRUE} otherwise.
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
#' cov2.2012LC(smallX, smallY) # run the test
#' 
#' \dontrun{
#' ## comparison of biased and unbiased estimator
#' ## empirical Type 1 error 
#' niter   = 100
#' vec.slow = rep(0,niter)  # record p-values
#' vec.fast = rep(0,niter)
#' for (i in 1:niter){
#'   X = matrix(rnorm(500*25), ncol=10)
#'   Y = matrix(rnorm(500*25), ncol=10)
#'   
#'   vec.slow[i] = ifelse(cov2.2012LC(X,Y,unbiased=TRUE)$p.value  < 0.05,1,0)
#'   vec.fast[i] = ifelse(cov2.2012LC(X,Y,unbiased=FALSE)$p.value < 0.05,1,0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* EMPIRICAL TYPE 1 ERROR COMPARISON \n","*\n",
#' "* Biased   case : ", round(sum(vec.fast/niter),5),"\n",
#' "* Unbiased case : ", round(sum(vec.slow/niter),5),"\n",sep=""))
#' }
#'
#' @references 
#' \insertRef{li_two_2012}{SHT}
#' 
#' @export
cov2.2012LC <- function(X, Y, unbiased=FALSE){
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
  
  X1 = as.matrix(scale(X, center = TRUE, scale = FALSE))
  X2 = as.matrix(scale(Y, center = TRUE, scale = FALSE))
  
  if (unbiased==TRUE){ # unbiased / slower / full
    A1  = cpp_cov2_2012LC_computeA(X1)
    A2  = cpp_cov2_2012LC_computeA(X2)
    C12 = cpp_cov2_2012LC_computeC(X1, X2)
  } else {
    A1  = cpp_cov2_2012LC_biased_computeA(X1) # elements for test statistics
    A2  = cpp_cov2_2012LC_biased_computeA(X2)
    C12 = cpp_cov2_2012LC_biased_computeC(X1, X2)
  }
  
  Tn1n2 = (A1 + A2 - 2*C12)  # test statistic
  shat  = 2*(A1/n2 + A2/n1)  # variance term
  pvalue = pnorm((Tn1n2/shat), lower.tail = FALSE)
  
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Test for High-Dimensional Covariances by Li and Chen (2012)"
  Ha      = "two covariances are not equal."
  thestat = Tn1n2
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "statistic"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}