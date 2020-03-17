#include <RcppArmadillo.h>
#include "cpp_extras.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// 1. log for a numeric value
double mylog(double val){
  return(static_cast<double>(std::log(static_cast<float>(val))));
}

