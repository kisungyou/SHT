#' Two-Sample Hotelling's T-squared Test for Multivariate Mean
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
  if (!is.logical(paired)){   stop("* mean2.Hotelling : 'paired' should be a logical.")  }
  if (!is.logical(var.equal)){stop("* mean2.Hotelling : 'var.equal' should be a logical.")}
  
  ##############################################################
  # BRANCHING 
  if (paired==TRUE){
    if (nrow(X)!=nrow(Y)){
      stop("* mean2.Hotelling : for paired different test, number of observations for X and Y should be equal.")
    } 
    
    diff = X-Y
    mu0  = rep(0,ncol(X))
    tmpout = mean1.Hotelling(diff, mu0=mu0, alpha=alpha)
      
    hname  = "Two-Sample Hotelling's T-squared Test for Paired/Dependent Data."
    Ha     = "true mean are different."
    output = hypothesis(hname, tmpout$statistic, tmpout$alpha,
                        tmpout$p.value, Ha, tmpout$conclusion)
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
      t2      = (sum(as.vector(Rlinsolve::lsolve.bicgstab(Spool, vecdiff, verbose=FALSE)$x )*vecdiff))*(nx*ny/(nx+ny))
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
      
      t2      = (sum(as.vector(Rlinsolve::lsolve.bicgstab(Stilde, vecdiff, verbose=FALSE)$x)*vecdiff))
      t2adj   = (t2*(v-p+1)/(v*p))
      pvalue  = pf(t2adj,p,(v-p+1),lower.tail = FALSE)
    }
    if (pvalue < alpha){
      conclusion = "Reject Null Hypothesis."
    } else {
      conclusion = "Not Reject Null Hypothesis."
    }
    output = hypothesis(hname, t2, alpha,
                        pvalue, Ha, 
                        conclusion)
  }

  ##############################################################
  # REPORT
  return(output)
}
