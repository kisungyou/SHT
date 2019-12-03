#' Test for Equality of Two Distributions by Biswas and Ghosh (2014)
#' 
#' Given two samples (either univariate or multivariate) \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : F_X = F_Y\quad vs\quad H_1 : F_X \neq F_Y}
#' using the procedure by Biswas and Ghosh (2014) in a nonparametric way based on 
#' pairwise distance measures. Both asymptotic and permutation-based determination of 
#' \eqn{p}-values are supported.
#' 
#' @param X a vector/matrix of 1st sample.
#' @param Y a vector/matrix of 2nd sample.
#' @param method method to compute \eqn{p}-value. Using initials is possible, \code{"p"} for permutation tests. Case insensitive.
#' @param nreps the number of permutations to be run when \code{method="permutation"}.
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
#' smallY = matrix(rnorm(10*3),ncol=3)
#' eqdist.2014BG(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## compare asymptotic and permutation-based powers
#' set.seed(777)
#' ntest  = 1000
#' pval.a = rep(0,ntest)
#' pval.p = rep(0,ntest)
#' 
#' for (i in 1:ntest){
#'   x = matrix(rnorm(100), nrow=5)
#'   y = matrix(rnorm(100), nrow=5)
#'   
#'   pval.a[i] = ifelse(eqdist.2014BG(x,y,method="a")$p.value<0.05,1,0)
#'   pval.p[i] = ifelse(eqdist.2014BG(x,y,method="p",nreps=100)$p.value <0.05,1,0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* EMPIRICAL TYPE 1 ERROR COMPARISON \n","*\n",
#' "* Asymptotics : ", round(sum(pval.a/ntest),5),"\n",
#' "* Permutation : ", round(sum(pval.p/ntest),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{biswas_nonparametric_2014}{SHT}
#' 
#' @export
eqdist.2014BG <- function(X, Y, method=c("asymptotic","permutation"), nreps=2000){
  ##############################################################
  # PREPROCESSING : accept both vector and matrix-valued data
  if (is.vector(X)){
    check_1d(X)
    check_1d(Y)
    X = matrix(X) 
    Y = matrix(Y)
  } else {
    check_nd(X)
    check_nd(Y)
    if (ncol(X)!=ncol(Y)){
      stop("* eqdist.2014BG : two input matrices should have same number of columns.")
    }
  }
  mymethod = tolower(method)
  if (pracma::strcmp(mymethod,"a")){
    mymethod = "asymptotic"
  } else if (pracma::strcmp(mymethod,"p")){
    mymethod = "permutation"
  }
  mymethod = match.arg(mymethod, c("asymptotic","permutation"))
  nreps = as.integer(nreps)
  
  ##############################################################
  # PARAMETERS AND COMPUTING TOTAL DISTANCE AND STATISTIC
  m = nrow(X)
  n = nrow(Y)
  
  XY  = rbind(X,Y) # smart way of doing this is to concatenate all data
  DXY = as.matrix(dist(XY))
  DX0 = DXY[1:m,1:m]                   # under null
  DY0 = DXY[(m+1):(m+n),(m+1):(m+n)]
  DZ0 = DXY[1:m,(m+1):(m+n)]
  Tmn = R_eqdist_2014BG_statistic(DX0,DY0,DZ0)
      
  ##############################################################
  if (pracma::strcmp(mymethod,"permutation")){
    Tvec = rep(0,nreps)
    for (i in 1:nreps){
      idx = sample(1:(m+n), m, replace=FALSE)
      idy = setdiff(1:(m+n), idx)
      
      DX1 = DXY[idx,idx]
      DY1 = DXY[idy,idy]
      DZ1 = DXY[idx,idy]
      Tvec[i] = R_eqdist_2014BG_statistic(DX1,DY1,DZ1)
    }
    pvalue = sum(Tvec>=Tmn)/nreps
  } else if (pracma::strcmp(mymethod,"asymptotic")){
    lbd = (m/(m+n))
    S1  = cpp_eqdist_2014BG_computeS(DX0)
    S2  = cpp_eqdist_2014BG_computeS(DY0)
    
    sig2 = (m*S1 + n*S2)/(m+n)
    Tmnstar = ((m+n)*lbd*(1.0-lbd)/(2*sig2))*Tmn
    pvalue  = pchisq(Tmnstar, 1, lower.tail = FALSE)
  }
  
  
  
  ##############################################################
  # REPORT
  thestat = Tmn
  hname   = "Test for Equality of Two Distributions by Biswas and Ghosh (2014)"
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") 
  Ha    = "two distributions are not equal"
  names(thestat) = "Tmn"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}

#' @keywords internal
#' @noRd
R_eqdist_2014BG_statistic <- function(DX,DY,DXY){
  m = nrow(DXY)
  n = ncol(DXY)

  muff = sum(DX[upper.tri(DX)])/(m*(m-1)/2)
  mufg = sum(DXY)/(m*n)
  mugg = sum(DY[upper.tri(DY)])/(n*(n-1)/2)
  
  vec1 = c(muff,mufg)
  vec2 = c(mufg,mugg)
  output = sum((vec1-vec2)^2)
  return(output)
}



# for exportation given distance matrix -----------------------------------
#' @keywords internal
#' @noRd
eqdist.2014BG.givenD <- function(D, size1, size2, method=c("asymptotic","permutation"), nreps=2000){
  ##############################################################
  ## conversion setting
  if (is.vector(X)){
    check_1d(X)
    check_1d(Y)
    X = matrix(X) 
    Y = matrix(Y)
  } else {
    check_nd(X)
    check_nd(Y)
    if (ncol(X)!=ncol(Y)){
      stop("* eqdist.2014BG : two input matrices should have same number of columns.")
    }
  }
  mymethod = tolower(method)
  if (pracma::strcmp(mymethod,"a")){
    mymethod = "asymptotic"
  } else if (pracma::strcmp(mymethod,"p")){
    mymethod = "permutation"
  }
  nreps = as.integer(nreps)
  
  m = as.integer(size1)
  n = as.integer(size2)
  DXY = D
  DX = DXY[1:m,1:m]                   # under null
  DY = DXY[(m+1):(m+n),(m+1):(m+n)]
  DZ = DXY[1:m,(m+1):(m+n)]
  Tmn = R_eqdist_2014BG_statistic(DX,DY,DZ)
  
  ##############################################################
  if (pracma::strcmp(mymethod,"permutation")){
    Tvec = rep(0,nreps)
    for (i in 1:nreps){
      idx = sample(1:(m+n), m, replace=FALSE)
      idy = setdiff(1:(m+n), idx)
      
      DX1 = DXY[idx,idx]
      DY1 = DXY[idy,idy]
      DZ1 = DXY[idx,idy]
      Tvec[i] = R_eqdist_2014BG_statistic(DX1,DY1,DZ1)
    }
    pvalue = sum(Tvec>=Tmn)/nreps
  } else if (pracma::strcmp(mymethod,"asymptotic")){
    lbd = (m/(m+n))
    S1  = cpp_eqdist_2014BG_computeS(DX)
    S2  = cpp_eqdist_2014BG_computeS(DY)
    
    sig2 = (m*S1 + n*S2)/(m+n)
    Tmnstar = ((m+n)*lbd*(1.0-lbd)/(2*sig2))*Tmn
    pvalue  = pchisq(Tmnstar, 1, lower.tail = FALSE)
  }
  
  
  
  ##############################################################
  # REPORT
  thestat = Tmn
  hname   = "Test for Equality of Two Distributions by Biswas and Ghosh (2014)"
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") 
  Ha    = "two distributions are not equal"
  names(thestat) = "Tmn"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
  
}
  