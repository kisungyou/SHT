#' Brown-Forsythe Test for Homogeneity of Variance
#' 
#' 
#' @export
vark.BrownForsythe <- function(dlist, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 
  check_alpha(alpha)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  k        = length(dlist)
  ni       = unlist(lapply(dlist, length))
  N        = sum(ni)
  
  ##############################################################
  # COMPUTATION : BROWN-FORSYTHE
  zlist = list()
  for (i in 1:k){
    tgt  = dlist[[i]]
    ybar = stats::median(tgt)
    zlist[[i]] = abs(tgt-ybar)
  }
  z_meanvec = unlist(lapply(zlist, mean))
  z_meanall = mean(unlist(zlist))
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  term1 = 0
  for (i in 1:k){
    term1 = term1 + (ni[i])*((z_meanvec[i]-z_meanall)^2)
  }
  term2 = 0
  for (i in 1:k){
    vecdiff = (zlist[[i]]-z_meanvec[i])
    term2   = term2 + sum(vecdiff*vecdiff)
  }
  thestat = (((N-k)*term1)/((k-1)*term2))
  pvalue  = stats::pf(thestat,(k-1),(N-k),lower.tail = FALSE)
  Ha      = "at least one of equalities does not hold."
  
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }
  
  ##############################################################
  # REPORT
  hname  = "Brown-Forsythe Test for Homogeneity of Variance"
  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}

#for (i in 1:1000){data=list();data[[1]]=rnorm(100);data[[2]]=rnorm(300);data[[3]]=rnorm(500);if (vark.Levene(data)$p.value<0.05){output[i]=1}else{output[i]=0}}
