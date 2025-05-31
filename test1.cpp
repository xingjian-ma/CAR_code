#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat test(){
  mat A(4, 5, fill::randu);
  
  mat B = resize(A, 7, 6);
  return B;
}