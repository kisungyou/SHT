#' Two-sample Test for Multivariate Means by Thulin (2014)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \mu_x = \mu_y\quad vs\quad H_1 : \mu_x \neq \mu_y}
#' using the procedure by Thulin (2014) using random subspace methods. 
#' We did not enable parallel computing schemes for this in that it might incur 
#' huge computational burden since it entirely depends on random permutation scheme.
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param B the number of selected subsets for averaging. \eqn{B\geq 100} is recommended.
#' @param nreps the number of permutation iterations to be run.
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
#' mean2.2014Thulin(smallX, smallY, B=10, nreps=10) # run the test
#' 
#' \donttest{
#' ## Compare with 'mean2.2011LJW' 
#' ## which is based on random projection.
#' n = 33    # number of observations for each sample
#' p = 100   # dimensionality
#' 
#' X = matrix(rnorm(n*p), ncol=p)
#' Y = matrix(rnorm(n*p), ncol=p)
#' 
#' ## run both methods with 100 permutations
#' mean2.2011LJW(X,Y,nreps=100,method="m")  # 2011LJW requires 'm' to be set.
#' mean2.2014Thulin(X,Y,nreps=100)
#' }
#' 
#' @references 
#' \insertRef{thulin_highdimensional_2014}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.2014Thulin <- function(X, Y, B=100, nreps=1000){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.2014Thulin : two samples X and Y should be of same dimension.")
  }
  n1 = nrow(X)
  n2 = nrow(Y)
  nn = (n1+n2)
  seqnn = 1:nn
  
  B     = as.integer(B)
  nreps = as.integer(nreps)
  k     = floor((n1+n2-2)/2)
  Z     = rbind(X,Y)

  ##############################################################
  # MAIN COMPUTATION
  T2     = mean2_2014_Thulin(X,Y,k,B) # our main statistic
  T2perm = rep(0,nreps)
  for (i in 1:nreps){
    id1   = base::sample(seqnn, n1, replace=FALSE) # random permutation
    id2   = setdiff(seqnn, id1)
    T2perm[i] = mean2_2014_Thulin(Z[id1,], Z[id2,], k, B)
  }
  
  thestat = T2
  pvalue  = sum(T2perm>=T2)/nreps
  
  ##############################################################
  # FINALE
  hname   = "Two-sample Test for Multivariate Means by Thulin (2014)"
  Ha      = "true means are different."
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "T2"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}



# auxiliary ---------------------------------------------------------------
#' @keywords internal
#' @noRd
mean2_2014_Thulin <- function(X,Y,k,B){
  # we select k-indices
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  
  # compute
  T2rec = rep(0,B)
  for (i in 1:B){
    colk = sample(1:p, k, replace = FALSE)
    pX   = X[,colk]
    pY   = Y[,colk]
    
    xbar = as.vector(colMeans(pX)) #### computing standard hotelling's T^2
    ybar = as.vector(colMeans(pY))
    Shat = ((n1-1)*cov(pX) + (n2-1)*cov(pY))/(n1+n2-2)
    
    T2rec[i] = sum(as.vector(solve(Shat,(xbar-ybar)))*(xbar-ybar))*((n1*n2)/(n1+n2))
  }
  
  # output
  T2fin = base::mean(T2rec)
  return(T2fin)
}


xx = list()
xx[[1]] = rnorm(10)
xx[[2]] = rnorm(5)
xx[[3]] = rnorm(7)
xx[[4]] = rnorm(21)

xx[order(unlist(lapply(xx, length)))]
