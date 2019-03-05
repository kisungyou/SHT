# 01. aux_pinv       : PseudoInverse using SVD and NumPy Scheme
# 02. aux_trace      : trace of a matrix
# 03. aux_var        : sample variance for univariate data
# 04. aux_adjustvec  : adjust small values (in magnitude) by some number
#     aux_adjustmat  : adjust small values (in magnitude) by some number
# 05. aux_quadform   : compute x'Ax
# 06. aux_minusvec   : for one-sample test
# 07. aux_CVsplit    : generate CV splits as list
# 08. aux_scatter    : sample covariance without scaling
# 09. aux_trace      : trace of a matrix
# 10. aux_getinvroot : scaling for one-sample covariance test
# 11. aux_adjustcube : given multivariate data, it turns them to fit in the cube

# 01. PseudoInverse using SVD and NumPy Scheme ----------------------------
# https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Singular_value_decomposition_(SVD)
#' @keywords internal
#' @noRd
aux_pinv <- function(A){
  svdA      = base::svd(A)
  tolerance = (.Machine$double.eps)*max(c(nrow(A),ncol(A)))*as.double(max(svdA$d))
  
  idxcut    = which(svdA$d <= tolerance)
  invDvec   = (1/svdA$d)
  invDvec[idxcut] = 0
  
  output = (svdA$v%*%diag(invDvec)%*%t(svdA$u))
  return(output)
}

# 02. trace of a matrix ---------------------------------------------------
#' @keywords internal
#' @noRd
aux_trace <- function(A){
  return(sum(diag(A)))
}

# 03. sample variance (univariate) ----------------------------------------

# 04. aux_adjustvec -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_adjustvec <- function(x, val=1e-10){
 n = length(x)
 y = rep(0,n)
 
 idx_small = which(abs(x)<=val)
 if (length(idx_small)==0){
   return(x)
 } else {
   idx_keep  = setdiff(1:n, idx_small)
   
   y[idx_small] = sign(x[idx_small])*val
   y[idx_keep]  = x[idx_keep]
   return(y)
 }
}
#' @keywords internal
#' @noRd
aux_adjustmat <- function(A, val=1e-10){
  n = nrow(A)
  p = ncol(A)
  
  idx_small = which(abs(A)<=val)
  if (length(idx_small)==0){
    return(A)
  } else {
    B = A
    B[idx_small] = sign(A[idx_small])*val
    return(B)
  }
}

# 05. aux_quadform --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_quadform <- function(A,x){
  return(sum((as.vector(A%*%x))*x))
}

# 06. aux_minusvec --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_minusvec <- function(X,vec){
  n = nrow(X)
  p = ncol(X)
  Y = array(0,c(n,p))
  for (i in 1:n){
    Y[i,] = as.vector(X[i,])-vec
  }
  return(Y)
}


# 07. aux_CVsplit   : generate CV splits as list --------------------------
#' @keywords internal
#' @noRd
aux_CVsplit <- function(n, K){
  allseq= 1:n
  folds = cut(seq(1,n),breaks=K,labels=FALSE)
  idxx  = list()
  for (i in 1:K){
    idxx[[i]] = which(folds==i, arr.ind=TRUE)
  }
  idyy = list()
  for (i in 1:K){
    idyy[[i]] = setdiff(allseq, idxx[[i]])
  }
  
  output = list()
  output$large = idyy
  output$small = idxx
  return(output)
}


# 08. aux_scatter ---------------------------------------------------------
#' @keywords internal
#' @noRd
aux_scatter <- function(X){
  nX = as.matrix(scale(X, center=TRUE, scale=FALSE))
  return(t(nX)%*%nX)
}


# 09. aux_trace -----------------------------------------------------------
#' @keywords internal
#' @noRd
aux_trace <- function(A){
  return(sum(diag(A)))
}


# 10. aux_getinvroot ------------------------------------------------------
#' @keywords internal
#' @noRd
aux_getinvroot <- function(X){
  eigs = eigen(X)
  if (any(eigs$values < .Machine$double.eps*10)){
    stop("** The desired covariance 'Sigma0' is invalid.")
  }
  out = eigs$vectors %*% diag((eigs$values)^(-0.5)) %*% t(eigs$vectors)
  return(out)
}


# 11. aux_adjustcube  -----------------------------------------------------
#' @keywords internal
#' @noRd
aux_adjustcube <- function(X,lower,upper){
  n = nrow(X)
  p = ncol(X)
  Y = array(0,c(n,p))
  for (i in 1:p){
    a = lower[i]
    b = upper[i]
    
    tgt   = as.vector(X[,i])
    Y[,i] = tgt/(b-a) - a/(b-a)
  }
  return(Y)
}
