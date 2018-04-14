#' Analysis of Variance for Equality of Means
#' 
#' 
#' 
#' @export
meank.anova <- function(dlist, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 
  check_alpha(alpha)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY FOR USING ANOVA
  K = length(dlist)
  if (K < 2){
    stop("* meank.anova : we need at least 2 groups of data.")
  }
  labellist = list()
  for (i in 1:K){
    labellist[[i]] = rep(i,length(dlist[[i]]))
  }
  
  data  = unlist(dlist)
  group = as.factor(unlist(labellist))
  
  ##############################################################
  # COMPUTATION : USE AOV INTERFACE
  aovout = unlist(summary(aov(data~group)))
  
  ##############################################################
  # REPORT
  hname   = "Analysis of Variance for Equality of Means"
  Ha      = "at least one of equalities does not hold."
  thestat = as.double(aovout[7])
  pvalue  = as.double(aovout[9])
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }

  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}

