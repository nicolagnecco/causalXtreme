#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> cpp_stl_vec(int n) {
  std::vector<double> diagG(n);

  for (int i = 0; i < n; ++i){
    diagG[i] = i;
  }
  return diagG;
}

// [[Rcpp::export]]
Rcpp::NumericVector rcpp_stl_vec(int n) {
  Rcpp::NumericVector diagG(n);

  for (int i = 0; i < n; ++i){
    diagG[i] = i;
  }
  return diagG;
}


// [[Rcpp::export]]
NumericMatrix gibbs_cpp(int N, int thin) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < thin; j++) {
      x = rgamma(1, 3, 1 / (y * y + 4))[0];
      y = rnorm(1, 1 / (x + 1), 1 / sqrt(2 * (x + 1)))[0];
    }
    mat(i, 0) = x;
    mat(i, 1) = y;
  }

  return(mat);
}
