#' Temporary : Random Projection for Two-Sample Simultaneous Testing
#' 
#' @examples 
#' ## CRAN-purpose small example
#' smallX = matrix(rnorm(10*3),ncol=3)
#' smallY = matrix(rnorm(10*3),ncol=3)
#' sim2.2020YLL(smallX, smallY) # run the test
#' 
#' \donttest{
#' ## empirical Type 1 error 
#' niter   = 1000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(100*10), ncol=50)
#'   Y = matrix(rnorm(100*10), ncol=50)
#'   counter[i] = ifelse(sim2.2020YLL(X,Y,m=100)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'sim2.2020YLL'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @keywords internal
#' @noRd
sim2.2020YLL <- function(X, Y, m=50){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* sim2.2020YLL : two samples X and Y should be of same dimension.")
  }
  p = ncol(X)
  
  ##############################################################
  # LET'S RUN MULTIPLE ITERATIONS 
  rec.stat = rep(0,m)
  for (i in 1:m){
    projvec = rnorm(p)
    projvec = projvec/sqrt(sum(projvec*projvec))
    
    Xproj = as.vector(X%*%projvec) # projection onto 1-dimensional space
    Yproj = as.vector(Y%*%projvec)
    
    #------- this part is for defining statistics
    # try 1. LRT and 
    tmpval = mvar2.LRT(Xproj, Yproj)$statistic # this is -2*loglbd
    rec.stat[i] = (((tmpval/2)^(1/3)) - 8/9)/sqrt(1/9)
  }
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = max(rec.stat)
  pvalue  = 1-(pnorm(thestat, lower.tail=TRUE)^m)
  
  hname   = "Two-sample Simultaneous Test of Means and Covariances by You, Lin, and Lee (2020)"
  Ha      = "both means and covariances are not equal."
  
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  names(thestat) = "T"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}