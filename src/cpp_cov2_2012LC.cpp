#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// OPTION : unbiased=FALSE =====================================================
// [[Rcpp::export]]
double cov2_2012LC_A(arma::mat &X){
  // parameters
  int n = X.n_rows;
  double nh = static_cast<double>(n);
  double tmpval = 0.0;
  
  // first term.
  double sum1 = 0.0;
  double denom1 = nh*(nh-1.0);
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i!=j){
        tmpval = arma::dot(X.row(i), X.row(j));
        sum1 += tmpval*tmpval;
      }
    }
  }
  
  // second term
  double sum2 = 0.0;
  double denom2 = nh*(nh-1.0)*(nh-2.0)/2.0;
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      for (int k=0; k<n; k++){
        if ((i!=j)&&(j!=k)&&(i!=k)){
          sum2 += arma::dot(X.row(i), X.row(j))*arma::dot(X.row(j), X.row(k)); // (i,j,k)
        }
      }
    }
  }
  
  // third term
  double sum3 = 0.0;
  double denom3 = nh*(nh-1.0)*(nh-2.0)*(nh-3.0);
  for (int i=0;i<n;i++){
    for (int j=0;j<n;j++){
      for (int k=0;k<n;k++){
        for (int l=0;l<n;l++){
          if ((i!=j)&&(i!=k)&&(i!=l)&&(j!=k)&&(j!=l)&&(k!=l)){
            sum3 += arma::dot(X.row(i), X.row(j))*arma::dot(X.row(k), X.row(l)); // (i,j,k,l)
          }
        }
      }
    }
  }

  // return
  double output = (sum1/denom1) - (sum2/denom2) + (sum3/denom3);
  return(output);
}
// [[Rcpp::export]]
double cov2_2012LC_C(arma::mat &X, arma::mat &Y){
  // FOUR TERMS TO BE COMPUTED
  int n1 = X.n_rows;
  int n2 = Y.n_rows;

  // 1. (i,j)
  double tmpval = 0.0;
  double sum1 = 0.0;
  for (int i=0;i<n1;i++){
    for (int j=0;j<n2;j++){
      tmpval = arma::dot(X.row(i), Y.row(j));
      sum1 += tmpval*tmpval;
    }
  }
  double term1 = sum1/static_cast<double>(n1*n2);
  
  
  // 2. (i,k,j) to be subtracted
  double sum2 = 0.0;
  for (int i=0;i<n1;i++){
    for (int k=0;k<n1;k++){
      for (int j=0;j<n2;j++){
        if (i!=k){
          sum2 += arma::dot(X.row(i), Y.row(j))*arma::dot(Y.row(j), X.row(k));
        }
      }
    }
  }
  double term2 = sum2/static_cast<double>(n1*n2*(n1-1));

  // 3. (i,k,j) to be subtracted
  double sum3 = 0.0;
  for (int i=0;i<n1;i++){
    for (int k=0;k<n1;k++){
      for (int j=0;j<n2;j++){
        if (i!=k){
          sum3 += arma::dot(Y.row(i), X.row(j))*arma::dot(X.row(j), Y.row(k));
        }
      }
    }
  }
  double term3 = sum3/static_cast<double>(n1*n2*(n2-1));

  // 4. (i,k,j,l)
  double sum4 = 0.0;
  for (int i=0;i<n1;i++){
    for (int k=0;k<n1;k++){
      if (i!=k){
        for (int j=0;j<n2;j++){
          for (int l=0;l<n2;l++){
            if (j!=l){
              sum4 += arma::dot(X.row(i), Y.row(j))*arma::dot(X.row(k), Y.row(l));
            }
          }
        }
      }
    }
  }
  double term4 = sum4/static_cast<double>(n1*n2*(n1-1)*(n2-1));

  // return
  double output = term1-term2-term3+term4;
  return(output);
}



// OPTION : unbiased=TRUE ======================================================
// [[Rcpp::export]]
double cov2_2012LC_A_no_bias(arma::mat &X){
  int n = X.n_rows;
  double nh = static_cast<double>(n);
  
  // iterate over
  double sumh = 0.0;
  double tmpval = 0.0;
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      if (i!=j){
        tmpval = arma::dot(X.row(i), X.row(j));
        sumh  += tmpval*tmpval;
      }
    }
  }
  
  // return
  double output = sumh/static_cast<double>(nh*(nh-1.0));
  return(output);
}
// [[Rcpp::export]]
double cov2_2012LC_C_no_bias(arma::mat &X, arma::mat &Y){
  int n1 = X.n_rows;
  int n2 = Y.n_rows;
  
  // iterate
  double sumh = 0.0;
  double tmpval = 0.0;
  
  for (int i=0; i<n1; i++){
    for (int j=0; j<n2; j++){
      tmpval = arma::dot(X.row(i), Y.row(j));
      sumh  += tmpval*tmpval;
    }
  }
  
  // return
  double output = sumh/static_cast<double>(n1*n2);
  return(output);
}




