#' One-Sample Test for High-Dimensional Mean by Srivastava and Du (2008)
#' 
#' 
#' 
#' @references 
#' \insertRef{srivastava_test_2008}{SHT}
#' 
#' @author Kisung You
#' @export
mean1.Srivastava2008 <- function(X, mu0=rep(0,ncol(X)), alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_1d(mu0)      
  check_alpha(alpha)
  if (length(mu0)!=ncol(X)){
    stop("* mean1.Srivastava2008 : mu0 does not have consistent size as data.")
  }
  Xnew = aux_minusvec(X,mu0) # now the problem becomes testing agains (0,0,...,0)
  
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  N      = nrow(Xnew)
  n      = (N-1)
  p      = ncol(Xnew)
  xbar   = colMeans(Xnew)
  S      = stats::cov(Xnew)    # sample covariance
  R      = stats::cov2cor(S)
  trR2   = aux_trace((R%*%R))
  Dsinv  = diag(1/aux_adjustvec(diag(S)))
  cpn    = (1+(trR2/(p^(3/2))))
  
  term1 = (N*aux_quadform(Dsinv,xbar))-((n*p)/(n-2))
  term2 = sqrt((2*(trR2-((p^2)/n)))*cpn)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "One-Sample Test for High-Dimensional Mean by Srivastava and Du (2008)."
  Ha      = "true mean is different from mu0."
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis."
  } else {
    conclusion = "Not Reject Null Hypothesis."
  }
  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}