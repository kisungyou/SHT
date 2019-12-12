#' Temporary : Random Projection for One-Sample Simultaneous Testing
#' 
#' 
#' 
#' @examples 
#' ## CRAN-purpose small example
#' smallX = matrix(rnorm(10*3),ncol=3)
#' sim1.2020YLL(smallX) # run the test
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 10000
#' pvalues = rep(0,niter)
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*10), ncol=25)
#'   pvalues[i] = sim1.2020YLL(X)$p.value
#'   counter[i] = ifelse(pvalues[i] < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'sim1.2020YLL'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' 
#' ## visualize
#' vec.alpha = seq(from=0.01, to=0.99, length.out=100)
#' vec.error = rep(0,100)
#' for (i in 1:100){
#'    vec.error[i] = sum((pvalues <= vec.alpha[i]))/length(pvalues)
#' }
#' opar <- par(pty="s")
#' plot(vec.alpha, vec.error, type="b", main="Type 1 Error")
#' abline(a=0, b=1, lwd=2, col="red")
#' par(opar)
#' }
#' 
#' @export
sim1.2020YLL <- function(X, mu0=rep(0,ncol(X)), Sigma0=diag(ncol(X)), m=50){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  n = nrow(X)
  p = ncol(X)

  ##############################################################
  # CENTER AND SCALE PROPERLY
  minusX = matrix(rep(mu0,n), ncol=p, byrow=TRUE)
  scaler = aux_getinvroot(Sigma0)
  Xnew   = (X-minusX)%*%scaler
  
  ##############################################################
  # LET'S RUN MULTIPLE ITERATIONS 
  rec.stat = rep(0,m)
  for (i in 1:m){
    projvec = rnorm(p)
    projvec = projvec/sqrt(sum(projvec*projvec))
    
    Xproj = as.vector(Xnew%*%projvec) # projection onto 1-dimensional space

    #------- this part is for defining statistics
    # try 1. LRT and 
    tmpval = mvar1.LRT(Xproj)$statistic
    rec.stat[i] = (((tmpval/2)^(1/3)) - 8/9)/sqrt(1/9)
  }
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = max(rec.stat)
  pvalue  = 1-(pnorm(thestat, lower.tail=TRUE)^m)
  
  hname   = "One-Sample Simultaneous Test of Mean and Covariance by You, Lin, and Lee (2020)"
  Ha      = "both mean and covariance are not equal to mu0 and Sigma0."
  
  DNAME = deparse(substitute(X))
  names(thestat) = "T"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}
