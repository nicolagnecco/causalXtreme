#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
arma::mat center_cols(const arma::mat & G){
  int N;
  arma::mat G2;

  N = G.n_rows;
  G2 = G - arma::repmat(arma::mean(G), N, 1);

  return G2;
}

// [[Rcpp::export]]
arma::mat center_rows(const arma::mat & G){
  int N;
  arma::mat G1;
  arma::mat G2;

  G1 = G.t();

  N = G1.n_rows;
  G2 = G1 - arma::repmat(arma::mean(G1), N, 1);

  return G2.t();
}

arma::mat scale_cols(const arma::mat & G){
  int N;
  arma::mat G2;

  N = G.n_rows;
  G2 = G / arma::repmat(arma::stddev(G), N, 1);

  return G2;
}

arma::mat scale_rows(const arma::mat & G){
  int N;
  arma::mat G1;
  arma::mat G2;

  G1 = G.t();

  N = G1.n_rows;
  G2 = G1 / arma::repmat(arma::stddev(G1), N, 1);

  return G2.t();
}

double mentapprc(const arma::rowvec & x){
  // define variables
  int n;
  arma::rowvec x_;
  double xstd, k1, k2, gamma, gaussianEntropy, negentropy, entropy;

  // compute size of vector
  n = x.size();

  // standardize vector
  x_ = scale_rows(center_rows(x));

  // constants we need
  xstd = arma::stddev(x);
  k1 = 36 / (8 * sqrt(3) - 9);
  gamma = 0.37457;
  k2 = 79.047;
  gaussianEntropy = log(2 * M_PI) / 2 + 0.5;

  // this is negentropy
  negentropy = k2 * pow((arma::mean(arma::log(arma::cosh(x_))) - gamma), 2)
    + k1 * pow(arma::mean(x_ % arma::exp(- pow(x_, 2) / 2)), 2);

  // this is entropy
  entropy = gaussianEntropy - negentropy + log(xstd);

  // return result
  return entropy;
}

// [[Rcpp::export]]
double pwlingc(const arma::mat & X) {
  // define variables
  int p, n;
  arma::mat X_;
  arma::mat C;
  arma::rowvec res1;
  arma::rowvec res2;
  double LR;

  // compute size of matrix
  p = X.n_rows;
  n = X.n_cols;

  // standardize matrix
  X_ = scale_rows(center_rows(X));

  // compute covariance matrix
  C = cov(X_.t());


  // general entropy-based method, for variables of any distribution
  res1 = X_.row(0) - C(0, 1) * X_.row(1);
  res2 = X_.row(1) - C(1, 0) * X_.row(0);
  LR = mentapprc(X_.row(0)) - mentapprc(X_.row(1))
    - mentapprc(res1) + mentapprc(res2);

  // Rcpp::Rcout << res1 << "\n"
              // << res2 << "\n";

  // return result
  return pow(std::min(0.0, LR), 2);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
*/
