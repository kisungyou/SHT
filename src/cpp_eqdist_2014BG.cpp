#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double cpp_eqdist_2014BG_statistic(arma::mat DX, arma::mat DY, arma::mat DXY){
  // 1. parameters
  int m = DXY.n_rows; double mm = static_cast<double>(m);
  int n = DXY.n_cols; double nn = static_cast<double>(n);
  
  // 2. compute sums
  double sumXX = 0.0;
  double sumYY = 0.0;
  double sumXY = 0.0;
  for (int i=0;i<(m-1);i++){
    for (int j=(i+1);j<m;j++){
      sumXX += DX(i,j);
    }
  }
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      sumXY += DXY(i,j);
    }
  }
  for (int i=0;i<(n-1);i++){
    for (int j=(i+1);j<n;j++){
      sumYY += DY(i,j);
    }
  }
  sumXX /= (mm*(mm-1.0))/2.0;
  sumYY /= (nn*(nn-1.0))/2.0;
  sumXY /= (mm*nn);
  
  // 3. compute statistic (just as the squared format)
  double output = (sumXX-sumXY)*(sumXX-sumXY) + (sumXY-sumYY)*(sumXY-sumYY);
  return(output);
}

// [[Rcpp::export]]
double cpp_eqdist_2014BG_computeS(arma::mat D){
  // 1. parameter
  int m = D.n_rows; double mm = static_cast<double>(m);
  
  // 2. two summations
  // 2-1. triplets
  double term1 = 0.0;
  for (int i=0;i<(m-2);i++){
    for (int j=(i+1);j<(m-1);j++){
      for (int k=(j+1);k<m;k++){
        term1 += D(i,j)*D(i,k);
      }
    }
  }
  double term1fin = term1/((mm*(mm-1.0)*(mm-2.0))/6.0);
  
  
  // term1 /= ((mm*(mm-1.0)*(mm-2.0))/6.0);
  // 2-2. pair with squared
  double term2 = 0.0;
  for (int i=0;i<(m-1);i++){
    for (int j=(i+1);j<m;j++){
      term2 += D(i,j);
    }
  }
  double term2fin = term2/((mm*(mm-1.0))/2.0);
  // term2 /= (mm*(mm-1.0))/2.0;
  
  // 3. return output
  // double output = term1-(term2*term2);
  double output = term1fin - (term2fin*term2fin);
  return(output);
}