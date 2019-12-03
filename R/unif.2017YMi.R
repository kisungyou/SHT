#' Multivariate Test of Uniformity based on Interpoint Distances by Yang and Modarres (2017) 
#' 
#' Given a multivariate sample \eqn{X}, it tests
#' \deqn{H_0 : \Sigma_x = \textrm{ uniform on } \otimes_{i=1}^p [a_i,b_i]  \quad vs\quad H_1 : \textrm{ not } H_0}
#' using the procedure by Yang and Modarres (2017). Originally, it tests the goodness of fit 
#' on the unit hypercube \eqn{[0,1]^p} and modified for arbitrary rectangular domain.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param type type of statistic to be used, one of \code{"Q1"},\code{"Q2"}, and \code{"Q3"}.
#' @param lower length-\eqn{p} vector of lower bounds of the test domain.
#' @param upper length-\eqn{p} vector of upper bounds of the test domain.
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
#' unif.2017YMi(smallX) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' ##   compare performances of three methods 
#' niter = 1234
#' rec1  = rep(0,niter) # for Q1
#' rec2  = rep(0,niter) #     Q2
#' rec3  = rep(0,niter) #     Q3
#' for (i in 1:niter){
#'   X = matrix(runif(50*10), ncol=50) # (n,p) = (10,50)
#'   rec1[i] = ifelse(unif.2017YMi(X, type="Q1")$p.value < 0.05, 1, 0)
#'   rec2[i] = ifelse(unif.2017YMi(X, type="Q2")$p.value < 0.05, 1, 0)
#'   rec3[i] = ifelse(unif.2017YMi(X, type="Q3")$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'unif.2017YMi'\n","*\n",
#' "* Type 1 error with Q1 : ", round(sum(rec1/niter),5),"\n",
#' "*                   Q2 : ", round(sum(rec2/niter),5),"\n",
#' "*                   Q3 : ", round(sum(rec3/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{yang_multivariate_2017}{SHT}
#' 
#' @export
unif.2017YMi <- function(X, type=c("Q1","Q2","Q3"), lower=rep(0,ncol(X)), upper=rep(1,ncol(X))){
  ##############################################################
  # PREPROCESSING
  check_nd(X)    # univariate vector
  n = nrow(X)
  d = ncol(X)
  
  if (d<2){
    stop("* unif.2017YMi : use univariate functions for univariate data.")
  }
  if ((any(lower!=rep(0,d)))||(any(upper!=rep(1,d)))){
    X = aux_adjustcube(X,lower,upper)  
  }
  ttype = match.arg(type)
  
  if (length(lower)!=d){
    stop("* unif.2017YMi : 'lower' should be a vector of length 'ncol(X)'.")
  }
  if (length(upper)!=d){
    stop("* unif.2017YMi : 'upper' should be a vector of length 'ncol(X)'.")
  }
  
  ##############################################################
  # COMPUTATION
  #   1. pairwise squared distance matrix
  D2 = as.matrix(dist(X))^2
  #   2. Q1 (mean distance)
  Q1.top = (2/(n*(n-1)))*sum(D2[upper.tri(D2)])
  Q1.bot = d*(2*n+3)/(90*n*(n-1))
  Q1     = ((Q1.top-(d/6))^2)/Q1.bot
  #   3. Q2 (variance distance)
  Q2.top = sum((D2[upper.tri(D2)]-(d/6))^2)*(2/(n*(n-1)))
  if (d==2){
    Q2.bot = ((989+202*(n-2))/56700)*(2/(n*(n-1)))
  } else if (d==3){
    Q2.bot = ((37+6*(n-2))/1050)*(2/(n*(n-1)))
  } else {
    Q2.bot = ((49*(d^2)/16200) + (101*d/37800) + 2*(n-2)*(((d^2)/16200)+(29*d/37800)))*(2/(n*(n-1)))
  }
  Q2 = ((Q2.top - (7*d/180))^2)/Q2.bot

  #   5. case by case reporting
  if (pracma::strcmp(ttype,"Q1")){
    thestat = Q1
    names(thestat) = "Q1"
    pvalue  = pchisq(Q1, 1, lower.tail = FALSE)
  } else if (pracma::strcmp(ttype, "Q2")){
    thestat = Q2
    names(thestat) = "Q2"
    pvalue  = pchisq(Q2, 1, lower.tail = FALSE)
  } else if (pracma::strcmp(ttype, "Q3")){
    thestat = Q1+Q2
    names(thestat) = "Q3"
    pvalue  = pchisq(thestat, 2, lower.tail = FALSE)
  } else {
    stop("* something wrong.")
  }
  
  ##############################################################
  # REPORT
  hname = "Multivariate Test of Uniformity based on Interpoint Distances by Yang and Modarres (2017)"
  DNAME = deparse(substitute(X))
  Ha    = paste("Sample ", DNAME, " does not follow uniform distribution.",sep="")
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}