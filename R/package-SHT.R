#' Statistical Hypothesis Testing Toolbox
#' 
#' Testing statistical hypotheses is one of the most significant inference problems in statistics. 
#' As time goes by, the community has witnessed surge of unorthodox settings such as high-dimensionality 
#' as well as atypical problems galore. We, though not complete nor ambitious, try to gather up 
#' some of frequently appearing tests as many as possible. Entire list of available tests 
#' can be seen at README file. 
#'
#' @docType package
#' @name SHT
#' @aliases SHT-package
#' @import Rdpack
#' @importFrom pracma pinv
#' @importFrom foreach "%dopar%" foreach registerDoSEQ
#' @importFrom parallel detectCores stopCluster makeCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom fastclime fastclime
#' @importFrom utils packageVersion
#' @importFrom stats pt sd cov pchisq pf median aov cov2cor pnorm
#' @importFrom Rcpp evalCpp
#' @useDynLib SHT
NULL

