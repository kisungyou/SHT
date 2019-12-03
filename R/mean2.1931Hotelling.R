#' Two-sample Hotelling's T-squared Test for Multivariate Means
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \mu_x = \mu_y\quad vs\quad H_1 : \mu_x \neq \mu_y}
#' using the procedure by Hotelling (1931).
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param paired a logical; whether you want a paired Hotelling's test.
#' @param var.equal a logical; whether to treat the two covariances as being equal.
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
#' mean2.1931Hotelling(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## generate two samples from standard normal distributions.
#' X = matrix(rnorm(50*5), ncol=5)
#' Y = matrix(rnorm(77*5), ncol=5)
#' 
#' ## run single test
#' print(mean2.1931Hotelling(X,Y))
#' 
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=5)
#'   Y = matrix(rnorm(77*5), ncol=5)
#'   
#'   counter[i] = ifelse(mean2.1931Hotelling(X,Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.1931Hotelling'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{hotelling_generalization_1931}{SHT}
#' 
#' @export
mean2.1931Hotelling <- function(X, Y, paired=FALSE, var.equal=TRUE){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  alpha = 0.05
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.1931Hotelling : two samples X and Y should be of same dimension.")
  }
  if (!is.logical(paired)){   stop("* mean2.1931Hotelling : 'paired' should be a logical.")  }
  if (!is.logical(var.equal)){stop("* mean2.1931Hotelling : 'var.equal' should be a logical.")}
  
  ##############################################################
  # BRANCHING 
  if (paired==TRUE){
    if (nrow(X)!=nrow(Y)){
      stop("* mean2.1931Hotelling : for paired different test, number of observations for X and Y should be equal.")
    } 
    
    diff = X-Y
    mu0  = rep(0,ncol(X))
    tmpout = mean1.1931Hotelling(diff, mu0=mu0)
      
    hname  = "Two-sample Hotelling's T-squared Test for Paired/Dependent Data."
    Ha     = "true means are different."
 
    thestat = tmpout$statistic
    DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
    names(thestat) = "T2"
    res   = list(statistic=thestat, p.value=tmpout$p.value, alternative = Ha, method=hname, data.name = DNAME)
    class(res) = "htest"
  } else {
    nx = nrow(X)
    ny = nrow(Y) 
    Sx = cov(X)
    Sy = cov(Y)
    xbar = as.vector(colMeans(X))
    ybar = as.vector(colMeans(Y))
    p    = ncol(X)
    Ha   = "true means are different."
    
    if (var.equal==TRUE){ # pooled variance-covariance is used
      vecdiff = (xbar-ybar)
      Spool   = ((nx-1)*Sx + (ny-1)*Sy)/(nx+ny-2)
      t2      = (sum(as.vector(solve(Spool, vecdiff))*vecdiff))*(nx*ny/(nx+ny))
      t2adj   = ((nx+ny-p-1)/(p*(nx+ny-2)))*t2
      pvalue  = pf(t2adj,p,(nx+ny-1-p),lower.tail = FALSE)
      hname   = "Hotelling's T-squared Test for Independent Samples with Equal Covariance Assumption."
    } else {
      hname   = "Hotelling's T-squared Test for Independent Samples with Unequal Covariance Assumption."
      S1tilde = Sx/nx
      S2tilde = Sy/ny
      Stilde  = (S1tilde+S2tilde)
      vecdiff = (xbar-ybar)
      Sinv    = aux_pinv(Stilde)
      Z1      = (S1tilde%*%Sinv)
      Z2      = (S2tilde%*%Sinv)
      Z1sq    = (Z1%*%Z1)
      Z2sq    = (Z2%*%Z2)
      
      v1      = (p+(p^2))
      v2      = (((aux_trace(Z1sq)+((aux_trace(Z1))^2))/nx)+((aux_trace(Z2sq)+((aux_trace(Z2))^2))/ny))
      v       = (v1/v2)
      
      t2      = (sum(as.vector(solve(Stilde, vecdiff))*vecdiff))
      t2adj   = (t2*(v-p+1)/(v*p))
      pvalue  = pf(t2adj,p,(v-p+1),lower.tail = FALSE)
    }
    # if (pvalue < alpha){
    #   conclusion = "Reject Null Hypothesis."
    # } else {
    #   conclusion = "Not Reject Null Hypothesis."
    # }
    
    thestat = t2
    DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
    names(thestat) = "T2"
    res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
    class(res) = "htest"
  }

  ##############################################################
  # REPORT
  return(res)
}
