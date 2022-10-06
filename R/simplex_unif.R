#' Probability Simplex : Tests of Uniformity
#' 
#' Given a data \eqn{X \in \mathbb{R}{n\times p}} such that its rows are 
#' vectors in a probability simplex, i.e., \eqn{x \in \Delta_{p-1}
#' =\lbrace z \in \mathbb{R}^p~|~z_j > 0, \sum_{i=1}^p z_i = 1 \rbrace,
#' } test whether the data is uniformly distributed.
#' 
#' @param X an \eqn{(n\times p)} data matrix where each row is an observation.
#' @param method (\emph{case-insensitive}) name of the method to be used, including\describe{
#' \item{LRT}{likelihood-ratio test with the Dirichlet distribution.}
#' \item{LRTsym}{likelihood-ratio test using the symmetric Dirichlet distribution (default).}
#' }
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
#' \donttest{
#' ## pseudo-uniform data generation
#' N = 100
#' P = 4
#' X = matrix(stats::rnorm(N*P), ncol=P)
#' for (n in 1:N){
#'   x = X[n,]
#'   x = abs(x/sqrt(sum(x^2)))
#'   X[n,] = x^2
#' }
#' 
#' ## run the tests
#' simplex.uniform(X, "LRT")
#' simplex.uniform(X, "lrtsym")
#' }
#' 
#' @concept simplex
#' @export
simplex.uniform <- function(X, method){
  # -------------------------------------------------------------------------
  # Input
  # data name
  DNAME = deparse(substitute(X)) # borrowed from HDtest
  
  # data
  X = check_simplex(X, "simplex.uniform")

  # choice of a method
  method.all = c("lrt","lrtsym")
  if (missing(method)){
    method.now = "lrtsym"
  } else {
    method.now = match.arg(tolower(method), method.all)
  }
  
  # -------------------------------------------------------------------------
  # Method branching
  output = switch(method.now,
                  "lrtsym" = simplex_uniform_lrtsym(X),
                  "lrt"    = simplex_uniform_lrt(X))
  output["data.name"] = DNAME
  return(output)
}



# auxiliary : simplex_uniform_*** methods names ---------------------------
#' @keywords internal
#' @noRd
simplex_uniform_lrtsym <- function(X){
  # parameters
  par_mle  = dirichlet_mle_symmauto(X)
  par_null = 1
  
  # log-lkd
  llkd_null = dirichlet_loglkd_symmetric(X, par_null)
  llkd_mle  = dirichlet_loglkd_symmetric(X, par_mle)
  
  # specifics
  thestat = -2*(llkd_null - llkd_mle)
  pvalue  = stats::pchisq(thestat, df=1, lower.tail = FALSE)
  
  # wrapping & return
  hname   = "Test Uniformity on Simplex : LRT-Symmetric Dirichlet."
  Ha      = "data is not uniformly distributed."
  
  names(thestat) = "lambda"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname)
  class(res) = "htest"
  return(res)
}

#' @keywords internal
#' @noRd
simplex_uniform_lrt <- function(X){
  # parameters
  K = base::ncol(X)
  par_mle  = dirichlet_mle_general(X)
  par_null = rep(1, K)
  
  # log-lkd
  llkd_null = dirichlet_loglkd_general(X, par_null)
  llkd_mle  = dirichlet_loglkd_general(X, par_mle)
  
  # specifics
  thestat = -2*(llkd_null - llkd_mle)
  pvalue  = stats::pchisq(thestat, df=K, lower.tail = FALSE)
  
  # wrapping & return
  hname   = "Test Uniformity on Simplex : LRT-Dirichlet."
  Ha      = "data is not uniformly distributed."
  
  names(thestat) = "lambda"
  res   = list(statistic=thestat, p.value=pvalue, alternative = Ha, method=hname)
  class(res) = "htest"
  return(res)  
}




# #simple test
# nit = 1000
# vec_sym <- rep(0, nit)
# vec_lrt <- rep(0, nit)
# for (i in 1:nit){
#   X <- Compositional::rdiri(100, c(1, 1, 1))
#   vec_sym[i] = simplex.uniform(X, method = "lrtsym")$p.value
#   vec_lrt[i] = simplex.uniform(X, method = "lrt")$p.value
#   print(paste0("* iteration ",i," complete."))
#   utils::flush.console()
# }
# par(mfrow=c(1,2))
# hist(vec_sym, main="sym")
# hist(vec_lrt, main="lrt")