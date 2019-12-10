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
# 12. aux_plrt       : mvar2.2012ZXC

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


# 12. aux_plrt ------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_plrt <- function(x, y, k = 500){
  # x is the vector containg the x sample values
  # y is the vector containg the y sample values
  # k is the number of Legendre quadrature points
  
  # The output of plrt includes the value lam of the likelihood ratio test statistic 
  # for testing the equality of two independent normal distributions, and the P-value
  # of the test P(lambda_{n,m} <= lam)
  
  # Generate the Legendre quadrature points
  # r contains the nodes, w contains the weights
  Legendre <- function(p) {
    x <- matrix(0, p, p)
    b <- c(1:(p - 1))^2
    b <- sqrt(b/(4 * b - 1))
    for(i in 2:p) {
      x[i, i - 1] <- b[i - 1]
    }
    x <- x + t(x)
    a <- eigen(x)
    u <- a$values
    v <- a$vectors
    w <- 2 * v[1,  ]^2
    cbind(u[p:1], w[p:1])
  }
  z <- Legendre(k)
  r <- z[, 1]
  w <- z[, 2]
  
  # Compute the likelihood ratio test statistic
  n <- length(x)
  m <- length(y)
  xb <- mean(x)
  yb <- mean(y)
  u <- mean(c(x, y))
  num <- (sum((x - xb)^2)/n)^(n/2) * (sum((y - yb)^2)/m)^(m/2)
  den <- ((sum((x - u)^2) + sum((y - u)^2))/(n + m))^((n + m)/2)
  lam <- num/den
  
  # Locate the intervals containing the two roots a and b, and find a and b
  u <- seq(0,1,0.005)
  v <- (exp(log(lam)+n/2*log(n)+m/2*log(m)-(n+m)/2*log(n+m)-n/2*log(u)))^(2/m)
  u <- u[1-u>v]
  fn <- function(w, ...)
  {
    1-w - (exp(log(lam)+n/2*log(n)+m/2*log(m)-(n+m)/2*log(n+m)-n/2*log(w)))^(2/m)
  }
  a <- uniroot(fn, lower = 1e-10, upper = u[1], tol=1e-8, lam = lam, n = n, m = m)$root
  b <- uniroot(fn, lower = u[1], upper = 1-1e-10, tol=1e-8, lam = lam, n = n, m = m)$root
  
  # Compute the double integral
  h1 <- (b - a)/2
  h2 <- (b + a)/2
  w1 <- h1 * r + h2
  d1 <- 1 - w1
  c1 <- (exp(log(lam)+n/2*log(n)+m/2*log(m)-(n+m)/2*log(n+m)-n/2*log(w1)))^(2/m)
  k1 <- (d1 - c1)/2
  k2 <- (d1 + c1)/2
  s <- rep(0,k)
  for(i in 1:k) {
    w2 <- k1[i] * r + k2[i]
    s[i] <- k1[i] * sum((w1[i]^((n - 1)/2 - 1) * w2^((m - 1)/2 - 1))/sqrt(1 - w1[i] - w2) * w)
  }                            
  v <- h1 * sum(s*w)
  pv <- 1 - exp(lgamma((n + m - 1)/2) - lgamma((n - 1)/2) - lgamma((m - 1)/2) - lgamma(0.5))*v
  
  return(list(lrt=lam, pvalue=pv))
}
