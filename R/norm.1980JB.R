#' Univariate Test of Normality by Jarque and Bera (1980)
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : x\textrm{ is from normal distribution} \quad vs\quad H_1 : \textrm{ not } H_0}
#' using a test procedure by Jarque and Bera (1980).
#' 
#' @param x a length-\eqn{n} data vector.
#' @param method method to compute \eqn{p}-value. Using initials is possible, \code{"a"} for asymptotic for example. Case insensitive.
#' @param nreps the number of Monte Carlo simulations to be run when \code{method="MC"}.
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
#' ## generate samples from uniform distribution
#' x = runif(28)
#' 
#' ## test with both methods of attaining p-values
#' test1 = norm.1980JB(x, method="a") # Asymptotics
#' test2 = norm.1980JB(x, method="m") # Monte Carlo 
#' 
#' @references 
#' \insertRef{jarque_efficient_1980}{SHT}
#' 
#' \insertRef{jarque_test_1987}{SHT}
#' 
#' @export
norm.1980JB <- function(x, method=c("asymptotic","MC"), nreps=2000){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  myrule = tolower(method)
  if (pracma::strcmp(myrule,"a")){
    myrule = "asymptotic"
  } else if (pracma::strcmp(myrule,"m")){
    myrule = "mc"
  }
  finrule = match.arg(myrule, c("asymptotic","mc"))
  nreps   = as.integer(nreps)
  
  ##############################################################
  # BRANCHING
  if (finrule=="asymptotic"){
    thestat = norm_1987JB_single(x)
    pvalue  = pchisq(thestat, df=2, lower.tail = FALSE)
  } else {
    tmpout  = norm_1987JB_mcarlo(x, nreps)
    thestat = tmpout$statistic
    pvalue  = tmpout$counts/nreps
  }
  
  ##############################################################
  # REPORT
  hname   = "Univariate Test of Normality by Jarque and Bera (1980)"
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  Ha    = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  names(thestat) = "JB"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
