#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List symm_eigen(const arma::mat & X) {
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, X, "standard");



  return Rcpp::List::create(
    Rcpp::Named("values") = S,
    Rcpp::Named("vectors") = V
  );
}

// [[Rcpp::export]]
Rcpp::List symm_eigen2(const arma::mat & X) {
  int k = X.n_rows - 1;
  sp_mat sX = conv_to<sp_mat>::from(X);
  arma::mat V;
  arma::vec S;
  arma::eigs_sym(S, V, sX, k);



  return Rcpp::List::create(
    Rcpp::Named("values") = S,
    Rcpp::Named("vectors") = V
  );
}
