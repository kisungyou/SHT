#' Two-sample Test for Multivariate Means by Johansen (1980)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \mu_x = \mu_y\quad vs\quad H_1 : \mu_x \neq \mu_y}
#' using the procedure by Johansen (1980) by adapting Welch-James approximation 
#' of the degree of freedom for Hotelling's \eqn{T^2} test.
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
#' mean2.1980Johansen(smallX, smallY) # run the test
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(mean2.1980Johansen(X,Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.1980Johansen'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{johansen_welchjames_1980}{SHT}
#' 
#' @export
mean2.1980Johansen <- function(X, Y){
  # First two parts are commonly available for 
  #   mean2.1965Yao
  #   mean2.1980Johansen
  #   mean2.1986NVM
  #   mean2.2004KY
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.1980Johansen : two samples X and Y should be of same dimension.")
  }
  p = ncol(X)
  
  ##############################################################
  # PARAMETERS AND PRELIMINARY COMPUTATIONS
  N1 = nrow(X); n1 = N1-1
  N2 = nrow(Y); n2 = N2-1
  
  x1bar = as.vector(colMeans(X)) # means
  x2bar = as.vector(colMeans(Y))
  xbardiff = (x1bar-x2bar)
  
  S1 = cov(X)/N1 # sample tilde' covariances
  S2 = cov(Y)/N2
  SS  = (S1+S2)
  
  S1inv = pracma::pinv(S1) # inverse of covariances
  S2inv = pracma::pinv(S2)
  SSinv = pracma::pinv(SS)
  
  T2 = aux_quadform(SSinv, xbardiff) # Hotelling's T statistic
  
  ##############################################################
  # SPECIFICS
  #   some numbers
  A1 = diag(p) - solve(S1inv+S2inv, S1inv)
  A2 = diag(p) - solve(S1inv+S2inv, S2inv)
  D  = ((1/n1)*(aux_trace(A1%*%A1) + (aux_trace(A1)^2)) + (1/n2)*(aux_trace(A2%*%A2) + (aux_trace(A2)^2)))/2
  v  = p*(p+2)/(3*D)
  q  = p+(2*D)-((6*D)/(p*(p-1)+2))
  
  #   2. adjust statistic and compute p-value
  thestat = T2
  T2adj   = T2/q
  pvalue  = pf(T2adj, p, v, lower.tail = FALSE)
  
  
  ##############################################################
  # REPORT
  hname   = "Two-sample Test for Multivariate Means by Johansen (1980)"
  Ha      = "true means are different."
  
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "T2"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}