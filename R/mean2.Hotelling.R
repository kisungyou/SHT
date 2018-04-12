#' Two-sample Hotelling's T-squared test
#' 
#' 
#' 
#' @export
mean2.Hotelling <- function(X, Y, alpha=0.05, paired=FALSE, var.equal=TRUE){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  check_alpha(alpha)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.Hotelling : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # BRANCHING 
  if (paired==TRUE){
    if (nrow(X)!=nrow(Y)){
      stop("* mean2.Hotelling : for paired different test, number of observations for X and Y should be equal.")
    } 
    
    diff = X-Y
    mu0  = rep(0,ncol(X))
    tmpout = mean1.Hotelling(diff, mu0=mu0, alpha=alpha)
      
    hname  = "Two-sample Hotelling's T-squared test for Paired/Dependent Data"
    Ha     = "true mean are different."
    output = hypothesis(hname, tmpout$statistic, tmpout$alpha,
                        tmpout$p.value, Ha, tmpout$conclusion)
  } else {
    n1 = nrow(X)
    n2 = nrow(Y) ##############################################################################
    
  }

  ##############################################################
  # REPORT
  return(output)
}
