#' Two-Sample Test for High-Dimensional Means by Cai et al (2014)
#' 
#' 
#' @references 
#' \insertRef{cai_two-sample_2014}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.Cai2014 <- function(X, Y, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  check_alpha(alpha)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.Cai2014 : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)

  x1 = colMeans(X)
  x2 = colMeans(Y)
  Sn = (((n1-1)*cov(X) + (n2-1)*cov(Y))/(n1+n2-2))
  
  xdiff2 = (as.vector(x1-x2)^2)
  diagSn = diag(Sn)
  diagSn[diagSn <= (100*.Machine$double.eps)]=(100*.Machine$double.eps)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat  = ((n1*n2)/(n1+n2))*max(xdiff2/diagSn)
  adjstat  = (-(2*log(p))+log(log(p)))
  stanstat = (thestat+adjstat)
  pvalue  = 1 - exp(-exp(-stanstat/2)/sqrt(3.1415926535897932384626433))
    
  hname   = "Two-Sample Test for High-Dimensional Means by Cai et al (2014)."
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