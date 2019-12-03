#' Two-sample Test for Multivariate Means by Lopes, Jacob, and Wainwright (2011)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \mu_x = \mu_y\quad vs\quad H_1 : \mu_x \neq \mu_y}
#' using the procedure by Lopes, Jacob, and Wainwright (2011) using random projection. 
#' Due to solving system of linear equations, we suggest you to opt for asymptotic-based 
#' \eqn{p}-value computation unless truly necessary for random permutation tests.
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param method method to compute \eqn{p}-value. \code{"asymptotic"} for using approximating null distribution, 
#' and \code{"MC"} for random permutation tests. Using initials is possible, \code{"a"} for asymptotic for example.
#' @param nreps the number of permutation iterations to be run when \code{method="MC"}.
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
#' smallX = matrix(rnorm(10*3),ncol=10)
#' smallY = matrix(rnorm(10*3),ncol=10)
#' mean2.2011LJW(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(10*20), ncol=20)
#'   Y = matrix(rnorm(10*20), ncol=20)
#'   
#'   counter[i] = ifelse(mean2.2011LJW(X,Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.2011LJW'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{lopes_more_2011}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.2011LJW <- function(X, Y, method=c("asymptotic","MC"), nreps=1000){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.2011LJW : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # BRANCHING
  myrule = tolower(method)
  if (strcmp(myrule,"a")){
    myrule = "asymptotic"
  } else if (strcmp(myrule,"b")){
    myrule = "mc"
  }
  finrule = match.arg(myrule, c("asymptotic","mc"))
  nreps   = as.integer(nreps)
  
  ##############################################################
  # BRANCHING
  p  = ncol(X)
  n1 = nrow(X)
  n2 = nrow(Y)
  n  = (n1+n2-2)
  if (p < (n/2)){
    stop("* mean2.2011LJW : the spirit of this method requires '2*p >= nrow(X)+nrow(Y)-2'.")
  }
  k  = floor((n1+n2-2)/2)  
  
  if (strcmp(finrule,"asymptotic")){
    Tk2 = mean2_2011LJW_statistic(X,Y)
    Tk2adj = ((n-k+1)/(k*n))*Tk2
    pvalue  = pf(Tk2adj, k, (n-k+1), lower.tail = FALSE)
  } else {
    Tk2 = mean2_2011LJW_statistic(X,Y)
    Tk2rec = rep(0,nreps)
    Z       = rbind(X,Y) # let's stack for ease
    for (i in 1:nreps){
      reode = sample(1:(n1+n2), (n1+n2), replace = FALSE) # random permutation
      id1   = reode[1:n1]
      id2   = reode[(n1+1):(n1+n2)]
      Tk2rec[i] = mean2_2011LJW_statistic(Z[id1,], Z[id2,])
    }
    pvalue = sum(Tk2rec>=Tk2)/nreps
  }
  
  thestat = Tk2
  ##############################################################
  # FINALE
  hname   = "Two-sample Test for Multivariate Means by Lopes, Jacob, and Wainwright (2011)"
  Ha      = "true means are different."
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "T2"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}




# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
mean2_2011LJW_statistic <- function(X, Y){
  p  = ncol(X)
  n1 = nrow(X)
  n2 = nrow(Y)
  n  = (n1+n2-2)
  Shat = ((n1-1)*cov(X) + (n2-1)*cov(Y))/(n1+n2-2)
  
  xbar = as.vector(colMeans(X))
  ybar = as.vector(colMeans(Y))
  
  k    = floor((n1+n2-2)/2)  
  Pk   = matrix(rnorm(p*k),nrow=p)
  Pkxy = as.vector(t(Pk)%*%(xbar-ybar)) # now it's a matrix
  
  
  solveLHS = t(Pk)%*%Shat%*%Pk
  solveRHS = Pkxy
  solved = tryCatch({
    as.vector(solve(solveLHS, solveRHS))
  }, error = function(e){
    as.vector(pracma::pinv(solveLHS)%*%solveRHS)
  }
  )
  Tk2     = ((n1*n2)/(n1+n2))*sum(Pkxy*solved)
  return(Tk2)
}
