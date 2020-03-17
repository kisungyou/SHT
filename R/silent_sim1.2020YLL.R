#' Temporary : Random Projection for One-Sample Simultaneous Testing
#' 
#' @keywords internal
#' @noRd
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
