#' Test for Equality of Means by Zhang and Xu (2009)
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \mu_1 = \cdots \mu_k\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Zhang and Xu (2009) by applying multivariate extension of Scheffe's method 
#' of transformation.
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample matrix of same dimension.
#' @param method a method to be applied for the transformed problem. \code{"L"} for \eqn{L^2}-norm based method, and 
#' \code{"T"} for Hotelling's test, which might fail due to dimensionality. Case insensitive.
#' 
#' @examples 
#' ## CRAN-purpose small example
#' tinylist = list()
#' for (i in 1:3){ # consider 3-sample case
#'   tinylist[[i]] = matrix(rnorm(10*3),ncol=3)
#' }
#' meank.2009ZX(tinylist) # run the test
#' 
#' \donttest{
#' ## test when k=5 samples with (n,p) = (100,20)
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   mylist = list()
#'   for (j in 1:5){
#'      mylist[[j]] = matrix(rnorm(100*10),ncol=10)
#'   }
#'   
#'   counter[i] = ifelse(meank.2009ZX(mylist, method="L")$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'meank.2009ZX'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @references 
#' \insertRef{zhang_ksample_2009}{SHT}
#' 
#' @export
meank.2009ZX <- function(dlist, method=c("L","T")){
  ##############################################################
  # PREPROCESSING
  check_dlistnd(dlist) 
  mymethod = tolower(method)
  method   = match.arg(mymethod, c("l","t"))

  ## parameters
  k = length(dlist)
  p = ncol(dlist[[1]])
  
  ##############################################################
  # DATA MANIPULATION
  # 1. reordering by increasing size
  newdat = dlist[order(unlist(lapply(dlist, nrow)))]
  n1     = nrow(newdat[[1]])
  # 2. take the process
  Y3d = array(0,c(n1,p,(k-1)))
  for (l in 2:k){
    Y3d[,,(l-1)] = meank_2009ZX_pairwise(newdat[[1]], newdat[[l]])
  }
  # 3. Z : this part can be a bit tricky
  q = as.integer((k-1)*p)
  Z = array(0,c(n1,q))
  for (i in 1:n1){
    Z[i,] = as.vector(t(Y3d[i,,]) ) # or transposed
  }
  
  ##############################################################
  # APPLY EITHER METHODS
  if (pracma::strcmp(method,"t")){ # Hotelling's
    tmpout = SHT::mean1.1931Hotelling(Z)
  } else if (pracma::strcmp(method, "l")) {
    tmpout = SHT::mean1.1996BS(Z)
  } else {
    stop("* meank.2009ZX : 'method' should be one of two options.")
  }
  
  ##############################################################
  # MODIFY AND RETURN
  hname = "Test for Equality of Means by Zhang and Xu (2009)"
  DNAME = deparse(substitute(dlist))
  Ha    = "one of equalities does not hold."
  
  tmpout$method      = hname
  tmpout$data.name   = DNAME
  tmpout$alternative = Ha
  return(tmpout)
}


# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
meank_2009ZX_pairwise <- function(X1,Xl){
  n1 = nrow(X1)
  nl = nrow(Xl)
  
  Xlbar = as.vector(colMeans(Xl))
  Xlbar_part = as.vector(colMeans(Xl[1:n1,]))
  Yl = array(0,c(n1,ncol(X1)))
  for (j in 1:n1){
    X1j = as.vector(X1[j,])
    Xlj = as.vector(Xl[j,])
    Yl[j,] = (X1j-Xlbar) + sqrt(n1/nl)*(Xlj-Xlbar_part)
  }
  return(Yl)
}
