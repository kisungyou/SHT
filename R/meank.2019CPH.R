#' Test for Equality of Means by Cao, Park, and He (2019)
#' 
#' Given univariate samples \eqn{X_1~,\ldots,~X_k}, it tests
#' \deqn{H_0 : \mu_1 = \cdots \mu_k\quad vs\quad H_1 : \textrm{at least one equality does not hold}}
#' using the procedure by Cao, Park, and He (2019). 
#' 
#' @param dlist a list of length \eqn{k} where each element is a sample matrix of same dimension.
#' @param method a method to be applied to estimate variance parameter. \code{"original"} for the estimator 
#' proposed in the paper, and \code{"Hu"} for the one used in 2017 paper by Hu et al. Case insensitive and initials can be used as well.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples 
#' ## CRAN-purpose small example
#' tinylist = list()
#' for (i in 1:3){ # consider 3-sample case
#'   tinylist[[i]] = matrix(rnorm(10*3),ncol=3)
#' }
#' meank.2019CPH(tinylist, method="o") # newly-proposed variance estimator
#' meank.2019CPH(tinylist, method="h") # adopt one from 2017Hu
#' 
#' \dontrun{
#' ## test when k=5 samples with (n,p) = (10,50)
#' ## empirical Type 1 error 
#' niter   = 10000
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   mylist = list()
#'   for (j in 1:5){
#'      mylist[[j]] = matrix(rnorm(10*50),ncol=50)
#'   }
#'   
#'   counter[i] = ifelse(meank.2019CPH(mylist)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'meank.2019CPH'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{cao_test_2019}{SHT}
#' 
#' @export
meank.2019CPH <- function(dlist, method=c("original","Hu")){
  ##############################################################
  # PREPROCESSING
  check_dlistnd(dlist) 
  mymethod = tolower(method)
  if (strcmp(mymethod,"o")){
    mymethod = "original"
  } else if (strcmp(mymethod,"h")){
    mymethod = "hu"
  }
  method   = match.arg(mymethod, c("original","hu"))
  
  ## parameters
  k = length(dlist)
  p = ncol(dlist[[1]])
  vec.nl = unlist(lapply(dlist, nrow)) # observations per class
  n      = sum(vec.nl)                 # total number of observations
  
  ##############################################################
  # PRELIMINARY COMPUTATION
  vec.xbar = (lapply(dlist, colMeans)) # sample means in a list
  vec.S    = (lapply(dlist, cov))      # sample covariances
  for (l in 1:k){
    nl = vec.nl[l]
    vec.S[[l]] = vec.S[[l]]*(nl-1)/nl
  }
  
  ##############################################################
  # COMPUTE 1 : Statistic T
  term1 = 0
  for (i in 1:k){
    ni    = vec.nl[i]
    tgtd  = dlist[[i]]     # target data
    gmat  = tgtd%*%t(tgtd) # gram matrix
    diag(gmat) = 0
    term1 = term1 + (sum(gmat)/2)*((n-ni)/(n*(ni-1)))
  }
  term2 = 0
  for (l in 1:(k-1)){
    nl   = vec.nl[l]
    barl = vec.xbar[[l]]
    for (s in (l+1):k){
      ns   = vec.nl[s]
      bars = vec.xbar[[s]]
      term2 = term2 + (2*nl*ns/n)*sum(barl*bars)
    }
  }
  thestat = term1-term2
  
  ##############################################################
  # COMPUTE 2 : RATIO-CONSISTENT VARIANCE
  if (strcmp(method,"original")){
    term1 = 0
    for (l in 1:k){
      nl = vec.nl[l]
      tgtdata = dlist[[i]]
      nlhalf = floor(nl/2)
      Sn1 = cov(tgtdata[1:nlhalf,])*((nlhalf-1)/nlhalf)
      Sn2 = cov(tgtdata[(nlhalf+1):nl,])*((nlhalf-1)/nlhalf)
      term1 = term1 + (nl*((n-nl)^2)/(nl-1))*sum(diag(Sn1%*%Sn2))
    }
    term2 = 0
    for (l in 1:(k-1)){
      Sl = vec.S[[l]]
      nl = vec.nl[l]
      for (s in (l+1):k){
        Ss = vec.S[[s]]
        ns = vec.nl[s]
        term2 = term2 + 2*(nl*ns*sum(diag(Sl%*%Ss)))
      }
    }
    sig2 = (2/(n^2))*(term1+term2)
  } else if (strcmp(method,"hu")){
    term1 = 0
    for (l in 1:k){
      nl = vec.nl[l]
      Sl = vec.S[[l]]
      term1 = term1 + ((nl*((n-nl)^2)*(nl-1))/((nl+1)*(nl-2)))*(sum(diag(Sl%*%Sl)) - (sum(diag(Sl))^2)/(nl-1))
    }
    term2 = 0
    for (l in 1:(k-1)){
      Sl = vec.S[[l]]
      nl = vec.nl[l]
      for (s in (l+1):k){
        Ss = vec.S[[s]]
        ns = vec.nl[s]
        term2 = term2 + 2*(nl*ns*sum(diag(Sl%*%Ss)))
      }
    }
    sig2 = (2/(n^2))*(term1+term2)
  }
  adjstat = thestat/sqrt(sig2)
  pvalue  = pnorm(adjstat, lower.tail = FALSE)
  
  ##############################################################
  # REPORT
  hname   = "Test for Equality of Means by Cao, Park, and He (2019)"
  DNAME = deparse(substitute(dlist))
  Ha    = "one of equalities does not hold."
  names(thestat) = "T"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname, data.name = DNAME)
  class(res) = "htest"
  return(res)
}