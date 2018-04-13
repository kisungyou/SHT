// 1. cpp_variance     : compute univariate sample variance

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


/*
 * 1. cpp_variance : compute univariate sample variance
 */
// [[Rcpp::export]]
double cpp_variance(const arma::vec& x){
  int n = x.n_elem;
  double nn = static_cast<double>(n);
  double output = 0.0;
  double xbar = mean(x);
  double diff = 0.0;
  
  for (int i=0;i<n;i++){
    diff = x(i)-xbar;
    output += (diff*diff);
  }
  output /= (nn-1);
  return(output);
}
