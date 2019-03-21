#' Robust Jarque-Bera Test of Univariate Normality by Gel and Gastwirth (2008)
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : x\textrm{ is from normal distribution} \quad vs\quad H_1 : \textrm{ not } H_0}
#' using a test procedure by Gel and Gastwirth (2008), which is a robustified version Jarque-Bera test.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param C1 a control constant. Authors proposed \eqn{C1=6} for nominal level of \eqn{\alpha=0.05}.
#' @param C2 a control constant. Authors proposed \eqn{C2=24} for nominal level of \eqn{\alpha=0.05}.
#' @param method method to compute \eqn{p}-value. Using initials is possible, \code{"a"} for asymptotic for example.
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
#' test1 = norm.2008RJB(x, method="a") # Asymptotics
#' test2 = norm.2008RJB(x, method="m") # Monte Carlo 
#' 
#' @references 
#' \insertRef{gel_robust_2008}{SHT}
#' 
#' @export
norm.2008RJB <- function(x, C1=6, C2=24, method=c("asymptotic","MC"), nreps=2000){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  myrule = tolower(method)
  if (myrule=="a"){
    myrule = "asymptotic"
  } else if (myrule=="m"){
    myrule = "mc"
  }
  finrule = match.arg(myrule, c("asymptotic","mc"))
  nreps   = as.integer(nreps)
  if ((C1<=0)||(length(C1)>1)){
    stop("* norm.2008RJB : 'C1' must be a nonnegative constant number.")
  }
  if ((C2<=0)||(length(C2)>1)){
    stop("* norm.2008RJB : 'C2' must be a nonnegative constant number.")
  }
  
  ##############################################################
  # BRANCHING
  if (finrule=="asymptotic"){
    thestat = norm_2008RJB_single(x, C1, C2)
    pvalue  = pchisq(thestat, df=2, lower.tail = FALSE)
  } else {
    tmpout  = norm_2008RJB_mcarlo(x, nreps, C1, C2)
    thestat = tmpout$statistic
    pvalue  = tmpout$counts/nreps
  }
  
  ##############################################################
  # REPORT
  hname   = "Robust Jarque-Bera Test of Univariate Normality by Gel and Gastwirth (2008)"
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  Ha    = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  names(thestat) = "RJB"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}