#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat adjust_clime(arma::mat omega){
  // parameters and output
  int p = omega.n_rows;
  arma::mat output(p,p,fill::zeros);
  
  // iterate
  double a = 0.0;
  double b = 0.0;
  
  double a2 = 0.0;
  double b2 = 0.0;
  
  double theval = 0.0;
  for (int i=0;i<(p-2);i++){
    for (int j=(i+1);j<(p-1);j++){
      a = omega(i,j); a2 = a*a;
      b = omega(j,i); b2 = b*b;
      
      if (a2<=b2){
        theval = a;
      } else {
        theval = b;
      }
      output(i,j) = a;
      output(j,i) = a;
    }
  }
  for (int i=0;i<p;i++){
    output(i,i) = omega(i,i);
  }
  return(output);
}