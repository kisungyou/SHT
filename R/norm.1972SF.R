#' Univariate Test of Normality by Shapiro and Francia (1972)
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : x\textrm{ is from normal distribution} \quad vs\quad H_1 : \textrm{ not } H_0}
#' using a test procedure by Shapiro and Francia (1972), which is an approximation to Shapiro and Wilk (1965).
#' 
#' @param x a length-\eqn{n} data vector.
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
#' x = rnorm(10)
#' norm.1972SF(x) # run the test
#' 
#' \donttest{
#' ## generate samples from several distributions
#' x = stats::runif(496)            # uniform
#' y = stats::rgamma(496, shape=2)  # gamma
#' z = stats::rlnorm(496)           # log-normal
#' 
#' ## test above samples
#' test.x = norm.1972SF(x) # uniform
#' test.y = norm.1972SF(y) # gamma
#' test.z = norm.1972SF(z) # log-normal
#' }
#' 
#' @references 
#' \insertRef{shapiro_approximate_1972}{SHT}
#' 
#' @export
norm.1972SF <- function(x){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  DNAME <- deparse(substitute(x))
  
  ##############################################################
  # COMPUTATION
  x = sort(x)   # in an increasing order
  n = length(x)
  n <- length(x)
  if ((n < 5 || n > 5000)){
    stop("* norm.1972SF : we only take care of sample size between (5,5000).")
  }
  y = qnorm(ppoints(n, a = 3/8))
  W =  cor(x, y)^2
  u = log(n)
  v = log(u)
  mu  = -1.2725 + 1.0521 * (v - u)
  sig = 1.0308 - 0.26758 * (v + 2/u)
  z   = (log(1 - W) - mu)/sig
  
  ##############################################################
  # REPORT
  thestat = W
  hname   = "Univariate Test of Normality by Shapiro and Francia (1972)"
  Ha    = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  names(thestat) = "W"
  pvalue = stats::pnorm(z, lower.tail = FALSE)
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}