#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double cpp_cov2_2012LC_biased_computeA(arma::mat X){
  int n = X.n_rows;
  int p = X.n_cols;
  double nh = static_cast<double>(n);
  
  double term1 = 1.0/(nh*(nh-1));
  double term2 = 0;
  
  arma::colvec Xi(p, fill::zeros);
  arma::colvec Xj(p, fill::zeros);
  double tmpval = 0.0;
  for (int i=0;i<p;i++){
    Xi = X.col(i);
    for (int j=0;j<p;j++){
      Xj = X.col(j);
      if (i!=j){
        tmpval = arma::dot(Xi, Xj);
        term2 += (tmpval*tmpval);
      }
    }
  }
  
  double A = term1*term2;
  return(A);
}

// [[Rcpp::export]]
double cpp_cov2_2012LC_biased_computeC(arma::mat X, arma::mat Y){
  int n1 = X.n_rows;
  int n2 = Y.n_rows;
  double n1n2 = static_cast<double>(n1*n2);
  int p  = X.n_cols;
  
  arma::colvec Xi(p,fill::zeros);
  arma::colvec Yj(p,fill::zeros);
  double tmpval = 0.0;
  double innsum = 0.0;
  for (int i=0;i<p;i++){
    Xi = X.col(i);
    for (int j=0;j<p;j++){
      Yj = Y.col(j);
      tmpval = arma::dot(Xi, Yj);
      innsum += (tmpval*tmpval);
    }
  }
  
  double C = innsum/n1n2;
  return(C);
}


// [[Rcpp::export]]
double cpp_cov2_2012LC_computeA(arma::mat X){
  int n = X.n_rows;
  int p = X.n_cols;
  double nh = static_cast<double>(n);
  
  // common things
  double tmpval = 0.0;
  
  // first term.
  double term1 = 0.0;
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if (i!=j){
        tmpval = arma::dot(X.col(i), X.col(j));
        term1 += tmpval*tmpval;
      }
    }
  }
  term1 /= nh*(nh-1.0);
  
  // second term to be subtracted
  double term2 = 0.0;
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      for (int k=0;k<p;k++){
        if ((i!=j)&&(j!=k)&&(i!=k)){
          term2 += arma::dot(X.col(i), X.col(j))*arma::dot(X.col(j), X.col(k)); // (i,j,k)
        }
      }
    }
  }
  term2 /= nh*(nh-1.0)*(nh-2.0)/2;
  
  // third term to be added
  double term3 = 0.0;
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      for (int k=0;k<p;k++){
        for (int l=0;l<p;l++){
          if ((i!=j)&&(i!=k)&&(i!=l)&&(j!=k)&&(j!=l)&&(k!=l)){
            term3 += arma::dot(X.col(i), X.col(j))*arma::dot(X.col(k), X.col(l)); // (i,j,k,l)
          }
        }
      }
    }
  }
  term3 /= nh*(nh-1.0)*(nh-2.0)*(nh-3.0);
  
  // return output
  double A = term1-term2+term3;
  return(A);
}

// [[Rcpp::export]]
double cpp_cov2_2012LC_computeC(arma::mat X, arma::mat Y){
  // FOUR TERMS TO BE COMPUTED
  int n1 = X.n_rows;
  int n2 = Y.n_rows;
  int p  = X.n_cols;
  
  // 1. (i,j)
  double tmpval = 0.0;
  double term1 = 0.0;
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      tmpval = arma::dot(X.col(i), Y.col(j));
      term1 += tmpval*tmpval;
    }
  }
  term1 /= static_cast<double>(n1*n2);
  
  // 2. (i,k,j) to be subtracted
  double term2 = 0.0;
  for (int i=0;i<p;i++){
    for (int k=0;k<p;k++){
      for (int j=0;j<p;j++){
        if (i!=k){
          term2 += arma::dot(X.col(i), Y.col(j))*arma::dot(Y.col(j), X.col(k));
        }
      }
    }
  }
  term2 /= static_cast<double>(n1*n2*(n1-1));
  
  // 3. (i,k,j) to be subtracted
  double term3 = 0.0;
  for (int i=0;i<p;i++){
    for (int k=0;k<p;k++){
      for (int j=0;j<p;j++){
        if (i!=k){
          term3 += arma::dot(Y.col(i), X.col(j))*arma::dot(X.col(j), Y.col(k));
        }
      }
    }
  }
  term3 /= static_cast<double>(n1*n2*(n2-1));
  
  // 4. (i,k,j,l)
  double term4 = 0.0;
  for (int i=0;i<p;i++){
    for (int k=0;k<p;k++){
      if (i!=k){
        for (int j=0;j<p;j++){
          for (int l=0;l<p;l++){
            if (j!=l){
              term4 += arma::dot(X.col(i), Y.col(j))*arma::dot(X.col(k), Y.col(l));
            }
          }
        }
      }
    }
  }
  term4 /= static_cast<double>(n1*n2*(n1-1)*(n2-1));

  // return
  double output = term1-term2-term3+term4;
  return(output);
}
