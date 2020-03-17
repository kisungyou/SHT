#ifdef _OPENMP
  #include <omp.h>
#endif
#include "RcppArmadillo.h"
#include "cpp_extras.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

double cov2_get_tau(arma::vec xi, arma::vec xj, double gamma){
  int n = xi.n_elem; double nn = static_cast<double>(n);

  double dotij = arma::dot(xi,xj);
  double term1 = arma::dot(xi,xi);
  double term2 = (dotij*dotij)/((1.0+gamma)*arma::dot(xj,xj));

  double output = (term1-term2)/nn;
  return(output);
}


// double cov2_get_tau(arma::vec xi, arma::vec xj, double gamma){
//   int nx = xi.n_elem;
//   double nnx = static_cast<double>(nx);
//   arma::mat Hj = (xj*xj.t())/arma::dot(xj,xj);
//   
//   return((arma::dot(xi,xi) - (arma::dot(Hj*xi,xi)/(1.0+gamma)))/nnx);
// }

// [[Rcpp::export]]
arma::mat cpp_cov2_mxPBF_single(arma::mat X, arma::mat Y, double a0, double b0, double gamma){
  // 1. get some parameters
  int n1 = X.n_rows; double nn1 = static_cast<double>(n1);
  int n2 = Y.n_rows; double nn2 = static_cast<double>(n2); 
  int p  = X.n_cols; 
  int n  = (n1+n2);  double nn  = static_cast<double>(n);
  
  // 2. prepare
  arma::vec Xi(n1,fill::zeros);
  arma::vec Yi(n2,fill::zeros);
  arma::vec Zi(n,fill::zeros);
  
  arma::vec Xj(n1,fill::zeros);
  arma::vec Yj(n2,fill::zeros);
  arma::vec Zj(n,fill::zeros);
  
  
  double term_const = 0.5*std::log(static_cast<float>(gamma/(1.0+gamma))) + Rf_lgammafn((nn1/2.0) + a0) + Rf_lgammafn((nn2/2.0) + a0) - Rf_lgammafn(nn/2.0 + a0) + a0*std::log(static_cast<float>(b0)) - Rf_lgammafn(a0);

  // 3. iterate !
  arma::mat logBFmat(p,p,fill::zeros);
  
  for (int i=0;i<p;i++){
    Xi = X.col(i);
    Yi = Y.col(i);
    Zi = arma::join_vert(Xi,Yi);
    for (int j=0;j<p;j++){
      Xj = X.col(j);
      Yj = Y.col(j);
      Zj = arma::join_vert(Xj,Yj);
      
      if (i!=j){
        double term1 = (nn1/2.0 + a0)*mylog(b0 + (nn1/2.0)*cov2_get_tau(Xi,Xj,gamma));
        double term2 = (nn2/2.0 + a0)*mylog(b0 + (nn2/2.0)*cov2_get_tau(Yi,Yj,gamma));
        double term3 = (nn/2.0  + a0)*mylog(b0 + (nn/2.0)*cov2_get_tau(Zi,Zj,gamma));
        
        logBFmat(i,j) = term_const-(term1+term2)+term3;
      }
    }
  }
  return(logBFmat);
}


// [[Rcpp::export]]
arma::mat cpp_cov2_mxPBF_multiple(arma::mat X, arma::mat Y, double a0, double b0, double gamma, int nCores){
  // 1. get some parameters
  int n1 = X.n_rows; double nn1 = static_cast<double>(n1);
  int n2 = Y.n_rows; double nn2 = static_cast<double>(n2); 
  int p  = X.n_cols; 
  int n  = (n1+n2);  double nn  = static_cast<double>(n);
  
  // 2. prepare
  // arma::vec Xi(n1,fill::zeros);
  // arma::vec Yi(n2,fill::zeros);
  // arma::vec Zi(n,fill::zeros);
  // 
  // arma::vec Xj(n1,fill::zeros);
  // arma::vec Yj(n2,fill::zeros);
  // arma::vec Zj(n,fill::zeros);
  
  
  double term_const = 0.5*mylog(gamma/(1.0+gamma)) + Rf_lgammafn((nn1/2.0) + a0) + Rf_lgammafn((nn2/2.0) + a0) - Rf_lgammafn(nn/2.0 + a0) + a0*mylog(b0) - Rf_lgammafn(a0);
  
  // 3. iterate !
  // 3-1. openmp implementation for tau's
  arma::mat taumatX(p,p,fill::zeros);
  arma::mat taumatY(p,p,fill::zeros);
  arma::mat taumatZ(p,p,fill::zeros);  
  
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nCores) collapse(2) shared(X,Y,p,gamma,taumatX,taumatY,taumatZ)
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if (i!=j){
        arma::vec Xi = X.col(i);
        arma::vec Yi = Y.col(i);
        arma::vec Zi = arma::join_vert(Xi,Yi);
        
        arma::vec Xj = X.col(j);
        arma::vec Yj = Y.col(j);
        arma::vec Zj = arma::join_vert(Xj,Yj);
        
        taumatX(i,j) = cov2_get_tau(Xi,Xj,gamma);
        taumatY(i,j) = cov2_get_tau(Yi,Yj,gamma);
        taumatZ(i,j) = cov2_get_tau(Zi,Zj,gamma);  
      }
    }
  }
  #else
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if (i!=j){
        arma::vec Xi = X.col(i);
        arma::vec Yi = Y.col(i);
        arma::vec Zi = arma::join_vert(Xi,Yi);
        
        arma::vec Xj = X.col(j);
        arma::vec Yj = Y.col(j);
        arma::vec Zj = arma::join_vert(Xj,Yj);
        
        taumatX(i,j) = cov2_get_tau(Xi,Xj,gamma);
        taumatY(i,j) = cov2_get_tau(Yi,Yj,gamma);
        taumatZ(i,j) = cov2_get_tau(Zi,Zj,gamma);  
      }
    }
  }
  #endif
  // 3-2. gather up !
  double term1, term2, term3;
  arma::mat logBFmat(p,p,fill::zeros);
  for (int i=0;i<p;i++){
    for (int j=0;j<p;j++){
      if (i!=j){
        term1 = (nn1/2.0 + a0)*mylog(b0 + (nn1/2.0)*taumatX(i,j));
        term2 = (nn2/2.0 + a0)*mylog(b0 + (nn2/2.0)*taumatY(i,j));
        term3 = (nn/2.0  + a0)*mylog(b0 + (nn/2.0)*taumatZ(i,j));  
        
        logBFmat(i,j) = term_const-(term1+term2)+term3;
      }
    }
  }
  
  // 3-3. return 
  return(logBFmat);
}