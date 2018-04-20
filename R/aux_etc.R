# 01. aux_pinv  : PseudoInverse using SVD and NumPy Scheme
# 02. aux_trace : trace of a matrix
# 03. aux_var   : sample variance for univariate data


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
#' @keywords internal
#' @noRd
aux_var <- function(x){
  return(cpp_variance(x))
}
