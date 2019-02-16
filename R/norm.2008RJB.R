#' Robust Jarque-Bera Test of Univariate Normality by Gel and Gastwirth (2008)
#' 
#' 
#' @examples 
#' \donttest{
#' ## generate samples from uniform distribution
#' x = runif(496)
#' 
#' ## test with both methods of attaining p-values
#' test1 = norm.2008RJB(x, method="a") # Asymptotics
#' test2 = norm.2008RJB(x, method="m") # Monte Carlo 
#' }
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