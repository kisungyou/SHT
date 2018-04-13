#' One-Sample Hotelling's T-squared Test for Multivariate Mean
#' 
#' 
#' 
#' @export
mean1.Hotelling <- function(X, mu0=rep(0,ncol(X)), alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_1d(mu0)      
  check_alpha(alpha)
  if (length(mu0)!=ncol(X)){
    stop("* mean1.Hotelling : mu0 does not have consistent size as data.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n      = nrow(X)
  p      = ncol(X)
  xbar   = colMeans(X)
  Sigma  = cov(X)   # sample covariance
  SigBar = Sigma/n  # divided by n
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  vecdiff = as.vector(xbar)-as.vector(mu0)
  t2      = sum(as.vector(Rlinsolve::lsolve.bicgstab(SigBar, vecdiff, verbose=FALSE)$x)*vecdiff) # test statistic
  t2adj   = ((n-p)*t2/(p*(n-1)))
  pvalue  = pf(t2adj,p,(n-p),lower.tail = FALSE)
  
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }

  ##############################################################
  # REPORT
  hname  = "One-Sample Hotelling's T-squared Test"
  Ha     = "true mean is different from mu0."
  output = hypothesis(hname, t2, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}