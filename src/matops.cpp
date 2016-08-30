#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat MATAVE2(arma::mat A, arma::mat B) {
  return (A + B)/2;
}

// [[Rcpp::export]]
arma::mat MATAVE3(arma::mat A, arma::mat B, arma::mat C) {
  return (A + B + C)/3;
}
