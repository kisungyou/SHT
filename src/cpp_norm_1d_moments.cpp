/*
 * LIST OF CURRENTLY AVAILABLE C-FUNCTIONS
 * 1. 1987JB 
 * 2. 2008RJB
 * 3. 1996AJB
 */

#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// extra 1 : compute k-th sample moment
double norm_samplemoment(arma::vec x, int order){
  int n = x.n_elem;          // number of elements
  double xm = arma::mean(x); // mean value
  double dod = static_cast<double>(order);
  
  double output = 0.0;
  for (int i=0;i<n;i++){
    output += std::pow(x(i)-xm, dod);
  }
  output /= static_cast<double>(n);
  return(output); 
}
// extra 2 : sample skewness "sqrt(b1)"
double norm_skewness(arma::vec x){
  double term1 = norm_samplemoment(x, 3);
  double term2 = norm_samplemoment(x, 2);
  
  double output = term1/std::pow(term2, 1.5);
  return(output);
}
// extra 3 : sample kurtosis "b2"
double norm_kurtosis(arma::vec x){
  double term1 = norm_samplemoment(x, 4);
  double term2 = norm_samplemoment(x, 2);
  
  double output = term1/(term2*term2);
  return(output);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1. 1987JB
//    norm_1987JB_single
//    norm_1987JB_mcarlo
// [[Rcpp::export]]
double norm_1987JB_single(arma::vec x){
  int n = x.n_elem; 
  
  double sskew = norm_skewness(x);
  double skurt = norm_kurtosis(x);
  
  double JB = static_cast<double>(n)*((sskew*sskew/6.0) + ((skurt-3.0)*(skurt-3.0)/24));
  return(JB);
}
// [[Rcpp::export]]
Rcpp::List norm_1987JB_mcarlo(arma::vec x, int nreps){
  // 1. compute the statistic
  double statistic = norm_1987JB_single(x);
  
  // 2. monte carlo simulation
  int counter = 0;
  int n = x.n_elem;
  arma::vec mcsample(n, fill::zeros); 
  double tmpstat = 0.0;
  for (int i=0;i<nreps;i++){
    mcsample.randn();  // 2-1. fill with Gaussian random numbers
    tmpstat = norm_1987JB_single(mcsample);
    if (tmpstat > statistic){
      counter += 1;
    }
  }
  
  // 3. report statistic and number of acceptances.
  return Rcpp::List::create(Rcpp::Named("statistic")=statistic,
                            Rcpp::Named("counts")=counter);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2. 2008RJB
//    norm_2008RJB_single
//    norm_2008RJB_mcarlo
// [[Rcpp::export]]
double norm_2008RJB_single(arma::vec x, double C1, double C2){
  double mypi = 3.14159265358979323846264338327950288419716939937510;
  int n = x.n_elem; 
  double dn = static_cast<double>(n);
  
  double mu3 = norm_samplemoment(x, 3);
  double mu4 = norm_samplemoment(x, 4);
  
  double Jn = 0.0;
  double Mn = arma::median(x);
  double diff = 0.0;
  for (int i=0;i<n;i++){
    diff = x(i)-Mn;
    if (diff >= 0){
      Jn += diff;
    } else {
      Jn -= diff;
    }
  }
  Jn *= std::sqrt(mypi/2.0)/dn;
  
  double RJB = (dn/C1)*std::pow(mu3/std::pow(Jn, 3.0),2.0) + (dn/C2)*std::pow((mu4/std::pow(Jn,4.0))-3.0, 2.0);
  return(RJB);
}
// [[Rcpp::export]]
Rcpp::List norm_2008RJB_mcarlo(arma::vec x, int nreps, double C1, double C2){
  // 1. compute the statistic
  double statistic = norm_2008RJB_single(x, C1, C2);
  
  // 2. monte carlo simulation
  int counter = 0;
  int n = x.n_elem;
  arma::vec mcsample(n, fill::zeros); 
  double tmpstat = 0.0;
  for (int i=0;i<nreps;i++){
    mcsample.randn();  // 2-1. fill with Gaussian random numbers
    tmpstat = norm_2008RJB_single(mcsample, C1, C2);
    if (tmpstat > statistic){
      counter += 1;
    }
  }
  
  // 3. report statistic and number of acceptances.
  return Rcpp::List::create(Rcpp::Named("statistic")=statistic,
                            Rcpp::Named("counts")=counter);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3. 1996AJB
//    norm_1996AJB_single
//    norm_1996AJB_mcarlo
// [[Rcpp::export]]
double norm_1996AJB_single(arma::vec x){
  int n = x.n_elem;
  double nn = static_cast<double>(n);
  
  // two terms
  double term1 = std::pow(norm_skewness(x), 2.0)*(nn+1.0)*(nn+3.0)/(6.0*(nn-2.0));
  double term21 = std::pow(norm_kurtosis(x) - (3.0*(nn-1.0)/(nn+1.0)), 2.0);
  double term22 = (24.0*nn*(nn-2)*(nn-3))/((nn+1)*(nn+1)*(nn+3)*(nn+5));
  
  double AJB = term1 + (term21/term22);
  return(AJB);
}
// [[Rcpp::export]]
Rcpp::List norm_1996AJB_mcarlo(arma::vec x, int nreps){
  // 1. compute the statistic
  double statistic = norm_1996AJB_single(x);
  
  // 2. monte carlo simulation
  int counter = 0;
  int n = x.n_elem;
  arma::vec mcsample(n, fill::zeros); 
  double tmpstat = 0.0;
  for (int i=0;i<nreps;i++){
    mcsample.randn();  // 2-1. fill with Gaussian random numbers
    tmpstat = norm_1996AJB_single(mcsample);
    if (tmpstat > statistic){
      counter += 1;
    }
  }
  
  // 3. report statistic and number of acceptances.
  return Rcpp::List::create(Rcpp::Named("statistic")=statistic,
                            Rcpp::Named("counts")=counter);
}

