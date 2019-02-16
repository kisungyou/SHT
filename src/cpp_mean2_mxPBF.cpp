#include "RcppArmadillo.h"
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double mean2_get_sigmasq(arma::vec x, double gamma){
  int n = x.n_elem;
  
  arma::vec vv(n, fill::ones);
  
  double term1 = arma::dot(x,x);
  double term2 = std::pow(arma::dot(x,vv), 2.0)/(arma::dot(vv,vv)*(1.0+gamma));
  
  double output = (term1-term2)/static_cast<double>(n);
  return(output);
}
  


// [[Rcpp::export]]
arma::vec cpp_mean2_mxPBF_single(arma::mat X, arma::mat Y, double a0, double b0, double gamma){
  // 1. get some parameters
  int n1 = X.n_rows; double nn1 = static_cast<double>(n1);
  int n2 = Y.n_rows; double nn2 = static_cast<double>(n2); 
  int p  = X.n_cols; 
  int n  = (n1+n2);  double nn  = static_cast<double>(n);
  
  // 2. pre-compute sigmasq
  arma::vec sigsqX(p, fill::zeros);
  arma::vec sigsqY(p, fill::zeros);
  arma::vec sigsqZ(p, fill::zeros);
  for (int i=0;i<p;i++){
    sigsqX(i) = mean2_get_sigmasq(X.col(i), gamma);
    sigsqY(i) = mean2_get_sigmasq(Y.col(i), gamma);
    sigsqZ(i) = mean2_get_sigmasq(arma::join_vert(X.col(i),Y.col(i)), gamma);
  }
  
  // 3. main computation
  arma::vec logBFvec(p,fill::zeros);
  double log_gammas = std::log((gamma/(1.0+gamma)));
  double term1, term2;
  for (int i=0;i<p;i++){
    term1 = 2.0*b0 + nn*sigsqZ(i);
    term2 = 2.0*b0 + nn1*sigsqX(i) + nn2*sigsqY(i);
    logBFvec(i) = 0.5*log_gammas + ((nn/2.0)+a0)*std::log((term1/term2));
  }
  return(logBFvec);
}

// [[Rcpp::export]]
arma::vec cpp_mean2_mxPBF_multiple(arma::mat X, arma::mat Y, double a0, double b0, double gamma, int nCores){
  // 1. get some parameters
  int n1 = X.n_rows; double nn1 = static_cast<double>(n1);
  int n2 = Y.n_rows; double nn2 = static_cast<double>(n2); 
  int p  = X.n_cols; 
  int n  = (n1+n2);  double nn  = static_cast<double>(n);
  
  // 2. pre-compute sigmasq
  arma::vec sigsqX(p, fill::zeros);
  arma::vec sigsqY(p, fill::zeros);
  arma::vec sigsqZ(p, fill::zeros);

  #pragma omp parallel for num_threads(nCores) shared(p,sigsqX,sigsqY,sigsqZ,X,Y,gamma)
  for (int i=0;i<p;i++){
    sigsqX(i) = mean2_get_sigmasq(X.col(i), gamma);
    sigsqY(i) = mean2_get_sigmasq(Y.col(i), gamma);
    sigsqZ(i) = mean2_get_sigmasq(arma::join_vert(X.col(i),Y.col(i)), gamma);
  }

  // 3. main computation
  arma::vec logBFvec(p,fill::zeros);
  double log_gammas = std::log((gamma/(1.0+gamma)));
  double term1, term2;
  for (int i=0;i<p;i++){
    term1 = 2.0*b0 + nn*sigsqZ(i);
    term2 = 2.0*b0 + nn1*sigsqX(i) + nn2*sigsqY(i);
    logBFvec(i) = 0.5*log_gammas + ((nn/2.0)+a0)*std::log((term1/term2));
  }
  return(logBFvec);
}