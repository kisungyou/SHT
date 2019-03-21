#' Univariate Test of Normality by Shapiro and Wilk (1965)
#' 
#' Given an univariate sample \eqn{x}, it tests
#' \deqn{H_0 : x\textrm{ is from normal distribution} \quad vs\quad H_1 : \textrm{ not } H_0}
#' using a test procedure by Shapiro and Wilk (1965). Actual computation of \eqn{p}-value 
#' is done via an approximation scheme by Royston (1992).
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
#' ## generate samples from several distributions
#' x = stats::runif(28)            # uniform
#' y = stats::rgamma(28, shape=2)  # gamma
#' z = stats::rlnorm(28)           # log-normal
#' 
#' ## test above samples
#' test.x = norm.1965SW(x) # uniform
#' test.y = norm.1965SW(y) # gamma
#' test.z = norm.1965SW(z) # log-normal
#' 
#' @references 
#' \insertRef{shapiro_analysis_1965}{SHT}
#' 
#' \insertRef{royston_approximating_1992}{SHT}
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