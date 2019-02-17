#' Test for Homogeneity of Variance by Schott (2007)
#' 
#' 
#' @examples 
#' \donttest{
#' ### generate k-sample from multivariate normal distribution
#' k = 2           # set 5-sample case
#' n = 100          # 100 observations in R^5
#' p = 20
#' mydata = list()
#' for (i in 1:k){
#'   mydata[[i]] = matrix(rnorm(n*p), ncol=p)
#' }
#' 
#' X = mydata[[1]]
#' Y = mydata[[2]]
#' 
#' covk.2007Schott(mydata)
#' cov2.2007Schott(X,Y)
#' 
#' }
#' @export
covk.2007Schott <- function(dlist){
  ##############################################################
  # PREPROCESSING
  check_dlistnd(dlist) 
  
  ##############################################################
  # PREPARATION
  g     = length(dlist)  # g-sample case
  p     = ncol(dlist[[1]])
  vec.n = unlist(lapply(dlist, nrow))-1
  vec.S = array(0,c(p,p,g))
  for (i in 1:g){
    vec.S[,,i] = stats::cov(dlist[[i]])
  }
  
  n = sum(vec.n)
  S = array(0,c(p,p))
  for (i in 1:g){
    S = S + ((vec.n[i]/n)*vec.S[,,i])
  }
  
  # tr(Si) and tr(Si^2)
  vec.trS  = rep(0,g)
  vec.trS2 = rep(0,g)
  for (i in 1:g){
    Si = vec.S[,,i]
    vec.trS[i]  = sum(diag(Si))
    vec.trS2[i] = sum(diag(Si%*%Si))
  }
  
  ##############################################################
  # COMPUTATION 1 : tnm
  tnm = 0
  for (i in 1:g){
    ni = vec.n[i]
    ei = (ni+2)*(ni-1)
    Si = vec.S[,,i]
    for (j in 1:g){
      nj = vec.n[j]
      ej = (nj+2)*(nj-1)
      Sj = vec.S[,,j]
      
      if (i<j){
        add1 = (1 - (ni-2)/ei)*vec.trS2[i]
        add2 = (1 - (nj-2)/ej)*vec.trS2[j]
        min1 = 2*sum(diag(Si%*%Sj))
        min2 = (ni/ei)*((vec.trS[i])^2)
        min3 = (nj/ej)*((vec.trS[j])^2)
        
        tnm = tnm + (add1+add2) - (min1+min2+min3)
      }
    }
  }
  
  ##############################################################
  # COMPUTATION 2 : variance of tnm
  a = ((n^2)/((n+2)*(n-1)))*(sum(diag(S%*%S)) - (1/n)*(sum(diag(S))^2))
  
  inner1 = 0
  for (i in 1:g){
    ni = vec.n[i]
    for (j in 1:g){
      nj = vec.n[j]
      if (i<j){
        inner1 = inner1 + (((ni+nj)/(ni*nj))^2)
      }
    }
  }
  inner2 = (g-1)*(g-2)*sum(1/(vec.n^2))
  theta  = 2.0*sqrt(inner1+inner2)*a
  
  pvalue = pnorm(tnm/theta, lower.tail = FALSE)
  return(pvalue)
}