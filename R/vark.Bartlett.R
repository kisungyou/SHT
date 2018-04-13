#' Bartlett's Test for Homogeneity of Variance
#' 
#' 
#' @export
vark.Bartlett <- function(dlist, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_dlist1d(dlist) 
  check_alpha(alpha)
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  k        = length(dlist)
  vec_n    = unlist(lapply(dlist, length))
  vec_Si2  = unlist(lapply(dlist, aux_var))
  
  N        = sum(vec_n)
  Sp2      = (sum((vec_n-1)*vec_Si2))/(N-k)
  
  
  
  ##############################################################
  # COMPUTATION : HYPOTHESIS and DETERMINATION
  term1 = (((N-k)*log(Sp2)) - (sum((vec_n-1)*log(vec_Si2))))
  term2 = (((sum(1/(vec_n-1)))-(1/(N-k)))/(3*(k-1))) + 1
  
  thestat = (term1/term2)
  pvalue  = stats::pchisq(thestat,(k-1),lower.tail = FALSE)
  Ha      = "at least one of equalities does not hold."
  
  if (pvalue < alpha){
    conclusion = "Reject Null Hypothesis"
  } else {
    conclusion = "Not Reject Null Hypothesis"
  }
  
  ##############################################################
  # REPORT
  hname  = "Bartlett\'s Test for Homogeneity of Variance"
  output = hypothesis(hname, thestat, alpha,
                      pvalue, Ha, 
                      conclusion)
  return(output)
}