#' Two-Sample Test for High-Dimensional Means by Srivastava and Du (2008)
#' 
#' 
#' 
#' @references 
#' \insertRef{srivastava_test_2008}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.Srivastava2008 <- function(X, Y, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  check_alpha(alpha)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.Srivastava2008 : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n1 = nrow(X)
  n2 = nrow(Y)
  n  = (n1+n2-2) # from highmean package
  p  = ncol(X)
  
  x1 = colMeans(X)
  x2 = colMeans(Y)
  Sn = (((n1-1)*cov(X) + (n2-1)*cov(Y))/n)
  
  dsvec  = aux_adjustvec(diag(Sn))
  Dsinv  = diag(1/dsvec)
  Dshalf = diag(1/sqrt(dsvec))
  
  R      = stats::cov2cor(Sn)
  R2     = (R%*%R)
  trR2   = aux_trace(R2)
  cpn    = (1+(trR2/(p^(3/2))))
  
  xdiff = as.vector(x1-x2)
  term1 = ((((n1*n2)/(n1+n2))*aux_quadform(Dsinv,xdiff))-((n*p)/(n-2)))
  term2 = sqrt((2*(trR2-((p^2)/n)))*cpn)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "Two-Sample Test for High-Dimensional Means by Srivastava and Du (2008)."
  Ha      = "true means are different."
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