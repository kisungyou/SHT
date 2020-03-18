#' Two-sample Test for High-Dimensional Means by Cai, Liu, and Xia (2014)
#' 
#' Given two multivariate data \eqn{X} and \eqn{Y} of same dimension, it tests
#' \deqn{H_0 : \mu_x = \mu_y\quad vs\quad H_1 : \mu_x \neq \mu_y}
#' using the procedure by Cai, Liu, and Xia (2014). Since the test is equivalent to testing 
#' \deqn{H_0 : \Omega(\mu_x - \mu_y)=0},
#' if \eqn{\Omega}, an inverse covariance (or precision), is not provided, it finds 
#' a sparse precision matrix via CLIME estimator. Also, if two samples 
#' are assumed to have different covariance structure, it uses weighting scheme for adjustment.
#' 
#' @param X an \eqn{(n_x \times p)} data matrix of 1st sample.
#' @param Y an \eqn{(n_y \times p)} data matrix of 2nd sample.
#' @param Omega precision matrix; if \code{NULL}, it applies CLIME estimation. Otherwise, 
#' a \eqn{(p\times p)} inverse covariance should be provided.
#' @param cov.equal a logical to determine homogeneous covariance assumption.
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
#' smallX = matrix(rnorm(10*3),ncol=3)
#' smallY = matrix(rnorm(10*3),ncol=3)
#' mean2.2014CLX(smallX, smallY) # run the test
#' 
#' \dontrun{
#' ## empirical Type 1 error 
#' niter   = 100
#' counter = rep(0,niter)  # record p-values
#' for (i in 1:niter){
#'   X = matrix(rnorm(50*5), ncol=10)
#'   Y = matrix(rnorm(50*5), ncol=10)
#'   
#'   counter[i] = ifelse(mean2.2014CLX(X, Y)$p.value < 0.05, 1, 0)
#' }
#' 
#' ## print the result
#' cat(paste("\n* Example for 'mean2.2014CLX'\n","*\n",
#' "* number of rejections   : ", sum(counter),"\n",
#' "* total number of trials : ", niter,"\n",
#' "* empirical Type 1 error : ",round(sum(counter/niter),5),"\n",sep=""))
#' }
#' 
#' @references 
#' \insertRef{cai_twosample_2014}{SHT}
#' 
#' @export
mean2.2014CLX <- function(X, Y, Omega=NULL, cov.equal = TRUE){
  ##############################################################
  # PREPROCESSING
  check_nd(X)
  check_nd(Y)
  if (ncol(X)!=ncol(Y)){
    stop("* mean2.2014CLX : two samples X and Y should be of same dimension.")
  }
  if ((length(Omega)==0)&&is.null(Omega)){
    clime.flag = TRUE
  } else {
    p = ncol(X)
    if ((is.matrix(Omega))&&(nrow(Omega)==p)&&(ncol(Omega)==p)){
      clime.flag = FALSE
    } else {
      warning(" mean2.2014CLX : provided 'Omega' is not a valid matrix. We use data-driven approach with 'CLIME' estimator.")
      clime.flag = TRUE
    }
  }

  # Parameters and Functions
  n1   = nrow(X)
  n2   = nrow(Y)
  p    = ncol(X)

  # extra for future report
  Ha    = "two means are different."
  DNAME = paste(deparse(substitute(X))," and ",deparse(substitute(Y)),sep="") # borrowed from HDtest
  
  ##############################################################
  # CASE 1 : Omega is 'known'
  if (clime.flag==FALSE){
    xbar = base::colMeans(X)
    ybar = base::colMeans(Y)
    zbar = as.vector(Omega%*%(xbar-ybar))
    
    wvec = diag(Omega)                 # correction
    eeps = (100*.Machine$double.eps)
    wvec[abs(wvec)<eeps] = eeps
    
    # compute
    thestat   = ((n1*n2)/(n1+n2))*max((zbar^2)/wvec)    # compute and report
    # threshold = 2*log(p) - log(log(p)) - log(pi) - 2*log(log(1/(1-alpha))) # just in case for future use
    pvalue    = 1-extreme_type1(thestat-2*log(p)+log(log(p)))
    
    # report
    
    hname   = "Two-sample Mean Test with Known Precision by Cai, Liu, and Xia (2014)."
    names(thestat) = "statistic"
    res     = list(statistic=thestat, p.value=pvalue, alternative=Ha, method=hname, data.name=DNAME)
  } else {
  # CASE 2 : Omega is 'unknown' : need to estimate using clime 
    if (cov.equal==TRUE){ # 2-1.
      S.pooled = ((n1-1)*cov(X) + (n2-1)*cov(Y))/(n1+n2)
      sink(tempfile())    ########################################################################### no CV yet
      fastout = fastclime(S.pooled, lambda.min = 0.01)
      sink()

      Omega.hat = fastout$icovlist[[length(fastout$icovlist)]]
      # Omega.hat = cv_fastclime_equalcov(X, Y)
      
      X.omega = X%*%Omega.hat
      Y.omega = Y%*%Omega.hat
      W0      = ((n1-1)*cov(X.omega) + (n2-1)*cov(Y.omega))/(n1+n2)
      
      eeps    = 100*.Machine$double.eps
      wvec    = diag(W0)
      # wvec[abs(wvec)<eeps] = eeps
      zbar    = Omega.hat%*% as.vector(base::colMeans(X) - base::colMeans(Y))
      
      thestat = ((n1*n2)/(n1+n2))*max(zbar^2/wvec)
      # threshold = 2*log(p) - log(log(p)) - log(pi) - 2*log(log(1/(1-alpha))) # just in case for future use
      pvalue    = 1-extreme_type1(thestat-2*log(p)+log(log(p)))
      
      # report
      
      hname   = "Two-sample Mean Test with Equal Precisions with Estimation by Cai, Liu, and Xia (2014)."
      names(thestat) = "statistic"
      res     = list(statistic=thestat, p.value=pvalue, alternative=Ha, method=hname, data.name=DNAME)
      
    } else {              # 2-2. unequal covariance
      S.pooled = ((n1-1)*cov(X)/n1 + ((n1/n2)*(n2-1)*cov(Y))/n2)
      sink(tempfile())    ########################################################################### no CV yet
      fastout = fastclime(S.pooled, lambda.min = 0.01)
      sink()
      
      Omega.hat = fastout$icovlist[[length(fastout$icovlist)]]
      
      X.omega = X%*%Omega.hat
      Y.omega = Y%*%Omega.hat
      W0      = (n1-1)*cov(X.omega)/n1 + (n2-1)*cov(Y.omega)/n2*(n1/n2)
      
      eeps    = 100*.Machine$double.eps
      wvec    = diag(W0)
      # wvec[abs(wvec)<eeps] = eeps
      zbar    = Omega.hat%*% as.vector(base::colMeans(X) - base::colMeans(Y))
      
      
      thestat   = n1*max(zbar^2/wvec)
      # threshold = 2*log(p) - log(log(p)) - log(pi) - 2*log(log(1/(1-alpha))) # just in case for future use
      pvalue    = 1-extreme_type1(thestat-2*log(p)+log(log(p)))
      
      # report
      
      hname   = "Two-sample Mean Test with Unequal Precisions with Estimation by Cai, Liu, and Xia (2014)."
      names(thestat) = "statistic"
      res     = list(statistic=thestat, p.value=pvalue, alternative=Ha, method=hname, data.name=DNAME)
    }
  }
  
  ##############################################################
  # set up the output
  class(res) = "htest"
  return(res)
}


#' @keywords internal
#' @noRd
extreme_type1 <- function(t){
  return(exp(-(1/sqrt(pi))*exp(-t/2)))
}

# #' @keywords internal
# #' @noRd
#' cv_fastclime_equalcov <- function(X, Y){
#'   K     = 5
#'   npath = 10
#'   
#'   p  = ncol(X)
#'   n1 = nrow(X); cvid_X = aux_CVsplit(n1, K)
#'   n2 = nrow(Y); cvid_Y = aux_CVsplit(n2, K)
#'   
#'   # parallel in the inner loop
#'   nCore = max(round(detectCores()/2), 1)
#'   
#'   # iterate to go over from lambda = 1 to 0.01 with length 50
#'   all.scores = array(0, c(K,npath)) # (K cross validation , 50 lambdas)
#'   all.lambda = seq(from=0.999, to=0.01, length.out=npath)
#'   for (i in 1:K){
#'     Xtrain = X[cvid_X$large[[i]],]
#'     Ytrain = Y[cvid_Y$large[[i]],]
#'     Xtest = X[cvid_X$small[[i]],]
#'     Ytest = Y[cvid_Y$small[[i]],]
#'     
#'     Spool.train = ((nrow(Xtrain)-1)*cov(Xtrain) + (nrow(Ytrain)-1)*cov(Ytrain))/(nrow(Xtrain) + nrow(Ytrain))
#'     Spool.test  = ((nrow(Xtest)-1)*cov(Xtest) + (nrow(Ytest)-1)*cov(Ytest))/(nrow(Xtest) + nrow(Ytest))
#'     sink(tempfile())    
#'     fastout = fastclime(Spool.train, lambda.min = 0.01)
#'     sink()
#'     
#'     
#'     cl = makeCluster(nCore)
#'     registerDoParallel(cl)
#'     itforeach=NULL
#'     Rs = foreach (itforeach=1:npath, .combine=cbind) %dopar% {
#'       tgticov = fastclime::fastclime.selector(fastout$lambdamtx, fastout$icovlist, all.lambda[itforeach])
#'       max(abs(tgticov$icov%*%Spool.test - diag(p)))
#'     }
#'     stopCluster(cl)
#'     all.scores[i,] = as.vector(Rs)
#'   }
#'   # which lambda has the least ? Not Squared
#'   fin.score = base::colSums(all.scores)
#'   lbd.final = all.lambda[which.min(fin.score)[1]]
#'   print(lbd.final)
#'   
#'   # finalize
#'   fin.S     = ((n1-1)*cov(X) + (n2-1)*cov(Y))/(n1+n2)
#'   sink(tempfile())
#'   fin.fast  = fastclime(fin.S, lambda.min=0.01)
#'   sink()
#'   fin.out   = fastclime::fastclime.selector(fin.fast$lambdamtx, fin.fast$icovlist, lbd.final)
#'   fin.Omega = as.matrix(fin.out$icov)
#'   return(fin.Omega)
#' }