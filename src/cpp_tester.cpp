#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


//////////////////////////////////////// test 1. lgamma of arma : works fine
// [[Rcpp::export]]
arma::vec testcpp_lgamma(arma::vec x){
  
  int n = x.n_elem;
  arma::vec output(n, fill::zeros);
  for (int i=0;i<n;i++){
    output(i) = Rf_lgammafn(x(i));
  }
  return(output);
}
