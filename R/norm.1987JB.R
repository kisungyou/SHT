#' Univariate Test of Normality by Jarque and Bera (1987)
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
#' test1 = norm.1987JB(x, method="a") # Asymptotics
#' test2 = norm.1987JB(x, method="m") # Monte Carlo 
#' }
#' 
#' @export
norm.1987JB <- function(x, method=c("asymptotic","MC"), nreps=2000){
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
    thestat = norm_1987JB_single(x)
    pvalue  = pchisq(thestat, df=2, lower.tail = FALSE)
  } else {
    tmpout  = norm_1987JB_mcarlo(x, nreps)
    thestat = tmpout$statistic
    pvalue  = tmpout$counts/nreps
  }
  
  ##############################################################
  # REPORT
  hname   = "Univariate Test of Normality by Jarque and Bera (1987)"
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  Ha    = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  names(thestat) = "JB"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
