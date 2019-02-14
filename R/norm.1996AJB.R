#' Adjusted Jarque-Bera Test of Univariate Normality by Urzua (1996)
#' 
#' 
#' 
#' 
#' 
#' 
#' @examples 
#' \donttest{
#' ## generate samples from uniform distribution
#' x = runif(496)
#' 
#' ## test with both methods of attaining p-values
#' test1 = norm.1996AJB(x, method="a") # Asymptotics
#' test2 = norm.1996AJB(x, method="m") # Monte Carlo 
#' }
#' 
#' @export
norm.1996AJB <- function(x, method=c("asymptotic","MC"), nreps=2000){
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
  
  ##############################################################
  # BRANCHING
  if (finrule=="asymptotic"){
    thestat = norm_1996AJB_single(x)
    pvalue  = pchisq(thestat, df=2, lower.tail = FALSE)
  } else {
    tmpout  = norm_1996AJB_mcarlo(x, nreps)
    thestat = tmpout$statistic
    pvalue  = tmpout$counts/nreps
  }
  
  ##############################################################
  # REPORT
  hname   = "Adjusted Jarque-Bera Test of Univariate Normality by Urzua (1996)"
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  Ha    = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  names(thestat) = "AJB"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
