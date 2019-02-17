#' Univariate Test of Normality by Shapiro and Wilk (1965)
#' 
#' 
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
#' test.x = norm.1965SW(x) # uniform
#' test.y = norm.1965SW(y) # gamma
#' test.z = norm.1965SW(z) # log-normal
#' }
#' 
#' @export
norm.1965SW <- function(x){
  ##############################################################
  # PREPROCESSING
  check_1d(x)        # univariate vector
  
  ##############################################################
  # MAIN CALL OF 'SHAPIRO.WILK'
  DNAME = deparse(substitute(x))
  tmp = stats::shapiro.test(x)
  Ha  = paste("Sample ", DNAME, " does not follow normal distribution.",sep="")
  tmp$method = "Univariate Test of Normality by Shapiro and Wilk (1965)"
  tmp$alternative = Ha
  tmp$data.name   = DNAME
  
  
  
  return(tmp)
}