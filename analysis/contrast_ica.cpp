#include <RcppArmadillo.h>
#include <math.h>
#include <chol_gaussc.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double contrast_ica(const arma::mat & X, double sigma, double eta, double kappa) {

  int i, N, m;
  arma::mat Rkappa;
  arma::mat sizes;

  N = X.n_cols;
  m = X.n_rows;

  Rcpp::NumericMatrix x = conv_to<NumericMatrix>::from(X);

  for (i = 0; i < m; ++i){
    Rcpp::List res = chol_gaussc(x.row(i) / sigma, 1, N * eta)
  }


  // return Rcpp::List::create(
  //   Rcpp::Named("G") = eig_sym(X));

  return 2;
}

