#include "RcppArmadillo.h"
#include <omp.h>

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

double cov2_get_tau(arma::vec xi, arma::vec xj, double gamma){
  int nx = xi.n_elem;
  double nnx = static_cast<double>(nx);
  arma::mat Hj = (xj*xj.t())/arma::dot(xj,xj);
  
  return((arma::dot(xi,xi) - (arma::dot(Hj*xi,xi)/(1.0+gamma)))/nnx);
}

// [[Rcpp::export]]
arma::mat cpp_cov2_mxPBF(arma::mat X, arma::mat Y, arma::mat Z, double a0, double b0, double gamma, int nCores){
  // 1. get some parameters
  int n1 = X.n_rows; double nn1 = static_cast<double>(n1);
  int n2 = Y.n_rows; double nn2 = static_cast<double>(n2); 
  int p  = X.n_cols; 
  int n  = (n1+n2); double nn = static_cast<double>(n1+n2);
  
  // 2. prepare
  arma::vec Xi(n1,fill::zeros);
  arma::vec Yi(n2,fill::zeros);
  arma::vec Zi(n,fill::zeros);
  
  arma::vec Xj(n1,fill::zeros);
  arma::vec Yj(n2,fill::zeros);
  arma::vec Zj(n,fill::zeros);
  
  double term_const = 0.5*log(gamma/(1.0+gamma)) + Rf_lgammafn((nn1/2.0) + a0) + Rf_lgammafn((nn2/2.0) + a0) - Rf_lgammafn(nn/2.0 + a0) + a0*log(b0) - Rf_lgammafn(a0);
  
  // 3. iterate !
  arma::mat logBFmat(p,p,fill::zeros);

  #pragma omp parallel for num_threads(nCores) collapse(2) shared(X,Y,Z,gamma,a0,b0,nn,nn1,nn2,logBFmat)
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      Xi = X.col(i);
      Yi = Y.col(i);
      Zi = Z.col(i);
      Xj = X.col(j);
      Yj = Y.col(j);
      Zj = Z.col(j);
      
      if (i!=j){
        double term1 = (nn1/2.0 + a0)*log(b0 + (nn1/2.0)*cov2_get_tau(Xi,Xj,gamma));
        double term2 = (nn2/2.0 + a0)*log(b0 + (nn2/2.0)*cov2_get_tau(Yi,Yj,gamma));
        double term3 = (nn/2.0  + a0)*log(b0 + (nn/2.0)*cov2_get_tau(Zi,Zj,gamma));
        
        logBFmat(i,j) = term_const-(term1+term2)+term3;
      }
    }
  }
  return(logBFmat);
}
  
