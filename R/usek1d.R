#' Apply k-sample tests for two univariate samples
#' 
#' Any \eqn{k}-sample method implies that it can be used for 
#' a special case of \eqn{k=2}. \code{usek1d} lets any \eqn{k}-sample tests 
#' provided in this package be used with two univariate samples \eqn{x} and \eqn{y}.
#' 
#' @param x a length-\eqn{n} data vector.
#' @param y a length-\eqn{m} data vector.
#' @param test.name character string for the name of k-sample test to be used.
#' @param ... extra arguments passed onto the function \code{test.name}.
#' 
#' @return a (list) object of \code{S3} class \code{htest} containing: \describe{
#' \item{statistic}{a test statistic.}
#' \item{p.value}{\eqn{p}-value \eqn{P(H_0|H_1)} under current setting.}
#' \item{alternative}{alternative hypothesis.}
#' \item{method}{name of the test.}
#' \item{data.name}{name(s) of provided sample data.}
#' }
#' 
#' @examples
#' \donttest{
#' ### compare two-means via anova and t-test
#' ### since they coincide when k=2
#' x = rnorm(50)
#' y = rnorm(50)
#' 
#' ### run anova and t-test
#' test1 = usek1d(x, y, "meank.anova")
#' test2 = mean2.ttest(x,y)
#' 
#' ### print the results of comparison
#' cat(paste("\n* Comparison of ANOVA and t-test'\n\n",
#' sprintf("* p-value from ANOVA   : %.4f\n",test1$p.value),
#' sprintf("*              t-test  : %.4f\n",test2$p.value),sep=""))
#' }
#' 
#' @export
usek1d <- function(x, y, test.name, ...){
  ##############################################################
  # PREPROCESSING
  check_1d(x)
  check_1d(y)
  if (!is.character(test.name)){
    stop("* usek1d : provide your 'test.name' as a valid character string.")
  }
  mydata = list()
  mydata[[1]] = x
  mydata[[2]] = y
  
  ##############################################################
  # RUN
  funcget = get(test.name)
  tmpout = funcget(mydata,...)
  
  ##############################################################
  # RETURN
  DNAME = paste(deparse(substitute(x))," and ",deparse(substitute(y)),sep="") # borrowed from HDtest
  tmpout$data.name = DNAME
  return(tmpout)
}