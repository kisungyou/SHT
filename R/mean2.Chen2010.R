#' Two-Sample Test for High-Dimensional Means by Chen and Qin (2014)
#' 
#' 
#' @references 
#' \insertRef{chen_two-sample_2010}{SHT}
#' 
#' @author Kisung You
#' @export
mean2.Chen2010 <- function(X, Y, alpha=0.05){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  check_alpha(alpha)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.Chen2010 : two samples X and Y should be of same dimension.")
  }
  
  ##############################################################
  # COMPUTATION : PRELIMINARY
  n1 = nrow(X)
  n2 = nrow(Y)
  p  = ncol(X)
  n  = (n1+n2-2)
  
  Sn  = ((((n1-1)*cov(X))+((n2-1)*cov(Y)))/n)
  trS = sum(diag(Sn))
  trcov2 = ((n^2)/((n+2)*(n-1)))*(sum(Sn^2)-((trS^2)/n))
  T1 = (X%*%t(X))
  T2 = (Y%*%t(Y))
  
  part1 = (sum(T1) - sum(diag(T1)))/(n1*(n1 - 1))
  part2 = (sum(T2) - sum(diag(T2)))/(n2*(n2 - 1))
  part3 = 2*sum(X %*% t(Y))/(n1*n2)
  
  term1 = (part1+part2-part3)
  term2 = sqrt((2/(n1*(n1 - 1)) + 2/(n2*(n2 - 1)) + 4/(n1*n2))*trcov2)
  
  ##############################################################
  # COMPUTATION : DETERMINATION
  thestat = (term1/term2)
  pvalue  = pnorm(thestat,lower.tail = FALSE)
  
  hname   = "Two-Sample Test for High-Dimensional Means by Chen and Qin (2010)."
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



apval_Chen2010_samecov <- function(sam1, sam2){
  n1 <- dim(sam1)[1]
  n2 <- dim(sam2)[1]
  p <- dim(sam1)[2]
  n <- n1 + n2 - 2
  sam.cov <- ((n1 - 1)*cov(sam1) + (n2 - 1)*cov(sam2))/n
  trS <- sum(diag(sam.cov))
  tr.cov2 <- n^2/((n + 2)*(n - 1))*(sum(sam.cov^2) - trS^2/n)
  T1 <- sam1 %*% t(sam1)
  T2 <- sam2 %*% t(sam2)
  P1 <- (sum(T1) - sum(diag(T1)))/(n1*(n1 - 1))
  P2 <- (sum(T2) - sum(diag(T2)))/(n2*(n2 - 1))
  P3 <- -2*sum(sam1 %*% t(sam2))/(n1*n2)
  T <- P1 + P2 + P3
  test.stat <- T/sqrt((2/(n1*(n1 - 1)) + 2/(n2*(n2 - 1)) + 4/(n1*n2))*tr.cov2)
  test.stat <- as.numeric(test.stat)
  pval <- 1 - pnorm(test.stat)
  names(pval) <- "Chen2010"
  out <- NULL
  out$sam.info <- c("n1" = n1, "n2" = n2, "p" = p)
  out$cov.assumption <- "the two groups have same covariance"
  out$method <- "asymptotic distribution"
  out$pval <- pval
  return(out)
}