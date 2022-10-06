# AUXILIARY FUNCTIONS FOR SIMPLEX CASE
# (01) dirichlet_loglkd_general   : logLKD for general   Dirichlet distribution
# (02) dirichlet_loglkd_symmetric : logLKD for symmetric Dirichlet distribution
# (03) dirichlet_mle_general      : MLE for general   Dirichlet with Newton's method.
# (04) dirichlet_mle_symmetric    : MLE for symmetric Dirichlet with Newton's method.
# (05) dirichlet_mle_symmauto     : MLE for symmetric Dirichlet with R's routine.




# (01) dirichlet_loglkd_general -------------------------------------------
#' @keywords internal
#' @noRd
dirichlet_loglkd_general <- function(X, alpha_vec){
  # get parameters
  N = base::nrow(X)
  K = base::ncol(X)
  vec_ci = base::colSums(base::log(X))
  
  # compute
  term1 = N*base::lgamma(base::sum(alpha_vec))
  term2 = -N*base::sum(base::lgamma(alpha_vec))
  term3 = base::sum((alpha_vec-1)*vec_ci)
  return(term1+term2+term3)
}

# (02) dirichlet_loglkd_symmetric -----------------------------------------
#' @keywords internal
#' @noRd
dirichlet_loglkd_symmetric <- function(X, alpha){
  # get parameters
  N = base::nrow(X)
  K = base::ncol(X)
  C = base::sum(base::log(X))
  
  # compute
  output = N*base::lgamma(alpha*K) - N*K*base::lgamma(alpha) + (alpha-1)*C
  return(output)
}

# (03) dirichlet_mle_general ----------------------------------------------
#' @keywords internal
#' @noRd
dirichlet_mle_general <- function(X){
  # get parameters
  N = base::nrow(X)
  K = base::ncol(X)
  
  # iteration criteria
  par_abstol  = sqrt(.Machine$double.eps)
  par_maxiter = 100
  
  # others
  mat_ones  = outer(rep(1,K), rep(1,K))
  log_pbar  = base::colMeans(base::log(X))
  alprob    = base::colMeans(X)
  tmp_p2    = mean(X[,1]^2)
  tmp_xsi   = (alprob[1]-tmp_p2)/(tmp_p2-(alprob[1]^2) )
  alpha_old = as.vector(tmp_xsi*alprob)
  
  # iterate
  for (it in 1:par_maxiter){
    # Newton : prep
    iter_g = N*base::digamma(base::sum(alpha_old)) - N*base::digamma(alpha_old) + N*log_pbar
    iter_z = N*base::trigamma(base::sum(alpha_old))
    iter_H = base::diag(-N*base::trigamma(alpha_old)) + iter_z
    
    # Newton : compute and adjust
    alpha_new = alpha_old - base::solve(iter_H, iter_g)
    alpha_new[alpha_new < .Machine$double.eps] <- 1e-12
    
    # update
    alpha_inc = max(abs(alpha_new - alpha_old))
    alpha_old = alpha_new
    
    # break
    if (alpha_inc < par_abstol){
      break
    }
  }
  
  # return
  return(alpha_old)
}


# (04) dirichlet_mle_symmetric --------------------------------------------
#' @keywords internal
#' @noRd
dirichlet_mle_symmetric <- function(X){
  # get parameters
  N = base::nrow(X)
  K = base::ncol(X)
  
  # iteration criteria
  par_abstol  = sqrt(.Machine$double.eps)
  par_maxiter = 100
  
  # other stuffs
  scalar_C = base::sum(base::log(X))
  
  # iteration
  alpha_old = 1.0
  for (it in 1:par_maxiter){
    # compute the numerator & denominator
    term_num = N*K*base::digamma(alpha_old*K) - N*K*base::digamma(alpha_old) + scalar_C
    term_den = N*(K^2)*base::trigamma(alpha_old*K) - N*K*base::trigamma(alpha_old)
    
    # update
    alpha_new = alpha_old - (term_num/term_den)
    alpha_inc = abs(alpha_old-alpha_new)
    alpha_old = alpha_new
    
    # break
    if (alpha_inc < par_abstol){
      break
    }
  }
  
  # return
  return(alpha_old)
}


# (05) dirichlet_mle_symmauto ---------------------------------------------
#' @keywords internal
#' @noRd
dirichlet_mle_symmauto <- function(X){
  # get parameters
  N = base::nrow(X)
  K = base::ncol(X)
  C = base::sum(base::log(X))
  
  # cost function
  fun_opt <- function(alpha){
    return(N*lgamma(alpha*K) - N*K*lgamma(alpha) + (alpha-1)*C)
  }
  
  # solve & return
  sol_auto = stats::optimize(fun_opt, lower=1e-10, upper=1e+4, maximum = TRUE)$maximum
  return(sol_auto)
}




# # compare : dirichlet_mle_general
# vec1 <- c(5, 7, 1, 3)
# vec2 <- rep(1,3)
# vec3 <- rep(5,3)
# 
# # varying
# x <- Compositional::rdiri( 100, vec1)
# Compositional::diri.est(x)$param
# sirt::dirichlet.mle(x)$alpha
# dirichlet_mle_general(x)


