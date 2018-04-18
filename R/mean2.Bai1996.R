#' Two-Sample Test for High-Dimensional Means by Bai and Saranadasa (1996)
#' 
#' 
#' 
#' @references 
#' \insertRef{bai_high_1996}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.Bai1996 <- function(X, Y, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  check_alpha(alpha)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.Bai1996 : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n1 = nrow(X)
  n2 = nrow(Y)
  n  = (n1+n2-2) # from highmean package
  
  x1 = colMeans(X)
  x2 = colMeans(Y)
  Sn = (((n1-1)*cov(X) + (n2-1)*cov(Y))/n)
  
  xdiff = as.vector(x1-x2)
  Bn2   = ((n^2)/((n+2)*(n-1)))*(aux_trace((Sn%*%Sn))-(1/n)*((aux_trace(Sn))^2))
  
  term1 = ((sum(xdiff*xdiff))*((n1*n2)/(n1+n2)) - aux_trace(Sn))
  term2 = (sqrt(Bn2)*sqrt(2*(n+1)/n))
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail=FALSE) # reject if (Z > thr_alpha)
  
  hname   = "Two-Sample Test for High-Dimensional Means by Bai and Saranadasa (1996)."
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
