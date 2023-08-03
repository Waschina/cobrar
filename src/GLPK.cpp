#include <iostream>
#include <RcppArmadillo.h>
#include <glpk.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

/*
 * Get glpk version number
 */
// [[Rcpp::export]]
Rcpp::CharacterVector getGLPKVersion() {
  return glp_version();
}


