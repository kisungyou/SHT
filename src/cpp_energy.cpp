#ifdef _OPENMP
  #include <omp.h>
#endif
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double energy_distance(arma::mat Xi, arma::mat Xj, double alpha, int nCores){
  // get parameters
  int ni = Xi.n_rows; double nni = static_cast<double>(ni);
  int nj = Xj.n_rows; double nnj = static_cast<double>(nj);

  // compute 1. Mij
  double Mij = 0.0;
  // 
  // #pragma omp parallel for num_threads(nCores) collapse(2) shared(Mij,Xi,Xj,ni,nj,alpha) reduction(+: Mij)
  // for (int i=0;i<ni;i++){
  //   for (int j=0;j<nj;j++){
  //     arma::rowvec xdiff = Xi.row(i) - Xj.row(j);
  //     Mij += std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
  //   }
  // }
  // 
  // Mij /= (nni*nnj);
  
  arma::mat matMij(ni,nj,fill::zeros);
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nCores) collapse(2) shared(matMij,Xi,Xj,ni,nj,alpha)
  for (int i=0;i<ni;i++){
    for (int j=0;j<nj;j++){
      arma::rowvec xdiff = Xi.row(i) - Xj.row(j);
      matMij(i,j) = std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
    }
  }
  #else
  for (int i=0;i<ni;i++){
    for (int j=0;j<nj;j++){
      arma::rowvec xdiff = Xi.row(i) - Xj.row(j);
      matMij(i,j) = std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
    }
  }
  #endif
  Mij = arma::accu(matMij)/(nni*nnj);
  
  // compute 2. Mii
  // double Mii = 0.0;
  // #pragma omp parallel for num_threads(nCores) collapse(2) shared(Mii,Xi,ni,alpha)
  // int i,j;
  // #pragma omp parallel for private(i) shared(Mii,Xi,ni,alpha,j) schedule(dynamic)
  // for (i=0;i<(ni-1);i++){
  //   #pragma omp parallel for private(j) shared(Mii,Xi,ni,alpha) schedule(dynamic)
  //   for (int j=0;j<ni;j++){
  //   
  //       arma::rowvec xdiff = Xi.row(i)-Xi.row(j);
  //       Mii += 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
  //   }
  // }
  // Mii /= (nni*nni);

  // arma::mat matMii(ni,ni,fill::zeros);
  // #pragma omp parallel for num_threads(nCores) collapse(2) shared(matMii,Xi,ni,alpha)
  // for (int i=0;i<ni;i++){
  //   for (int j=0;j<ni;j++){
  //     if (i<j){
  //       arma::rowvec xdiff = Xi.row(i)-Xi.row(j);
  //       matMii(i,j) = 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
  //     }
  //   }
  // }
  // Mii = arma::accu(matMii)/(nni*nni);
  
  double Mii = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nCores) collapse(2) shared(Xi,ni,alpha) reduction(+: Mii)
  for (int i=0;i<ni;i++){
    for (int j=0;j<ni;j++){
      if (i<j){
        arma::rowvec xdiff = Xi.row(i)-Xi.row(j);
        Mii += 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
      }
    }
  }
  #else
  for (int i=0;i<ni;i++){
    for (int j=0;j<ni;j++){
      if (i<j){
        arma::rowvec xdiff = Xi.row(i)-Xi.row(j);
        Mii += 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
      }
    }
  }
  #endif
  Mii /= (nni*nni);
  
  
  
  // compute 3. Mjj
  double Mjj = 0.0;
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nCores) collapse(2) shared(Xj,nj,alpha) reduction(+: Mjj)
  for (int i=0;i<nj;i++){
    for (int j=0;j<nj;j++){
      if (i<j){
        arma::rowvec xdiff = Xj.row(i)-Xj.row(j);
        Mjj += 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
      }
    }
  }
  #else
  for (int i=0;i<nj;i++){
    for (int j=0;j<nj;j++){
      if (i<j){
        arma::rowvec xdiff = Xj.row(i)-Xj.row(j);
        Mjj += 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
      }
    }
  }
  #endif
  Mjj /= (nnj*nnj);
  
  // arma::mat matMjj(nj,nj,fill::zeros);
  // #pragma omp parallel for num_threads(nCores) collapse(2) shared(matMjj,Xj,nj,alpha)
  // for (int i=0;i<nj;i++){
  //   for (int j=0;j<nj;j++){
  //     if (i<j){
  //       arma::rowvec xdiff = Xj.row(i)-Xj.row(j);
  //       matMjj(i,j) = 2.0*std::pow(arma::dot(xdiff,xdiff), alpha/2.0);
  //     }
  //   }
  // }
  // Mjj = arma::accu(matMjj)/(nnj*nnj);
  
  // finalize
  double output = ((nni*nnj)/(nni+nnj))*(2.0*Mij-Mii-Mjj);
  return(output);
}