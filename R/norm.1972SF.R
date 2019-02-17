#' Univariate Test of Normality by Shapiro and Francia (1972)
#' 
#' 
#' @examples 
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
  DNAME = deparse(substitute(x)) # borrowed from HDtest
  Ha    = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  names(thestat) = "W"
  pvalue = stats::pnorm(z, lower.tail = FALSE)
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}