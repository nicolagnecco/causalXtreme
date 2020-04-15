#include <RcppArmadillo.h>
#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// [[Rcpp::export]]
Rcpp::List chol_gaussc(Rcpp::NumericVector x, double sigma, double tol) {
  // Variable definition
  int m, n, nmax, i, iter, j, jast;
  double residual, a, b, c, maxdiagG;

  // Variable instantiation
  // m = x.nrow(); /* dimension of input space might be greater than 1*/
  n = x.size(); /* number of samples */
  nmax = 20 * 3  / 2;

  // Dynamic vectors
  Rcpp::NumericVector diagG(n);
  // std::vector<double> diagG(n);
  Rcpp::NumericVector G(nmax * n);
  Rcpp::IntegerVector pp(n);


  // Computaions
  for (i = 0; i < n; ++i){
    diagG[i] = 1;
  }

  for (i = 0; i < n; ++i){
    pp[i] = i;
  }

  iter=0;
  residual=n;
  jast=0;


  while (residual > tol){

    if (iter==(nmax-1)){
      /* need to reallocate memory to G */
      nmax += nmax / 2;
      Rcpp::NumericVector Gbis(nmax * n);
      for (i = 0; i < iter * n; ++i){
        Gbis[i] = G[i];
      }
      G = Gbis;
    }


    /* switches already calculated elements of G and order in pp */
    if (jast!=iter){
      i = pp[jast];
      pp[jast] = pp[iter];
      pp[iter] = i;

      for (i = 0;i <= iter; i++){
        a = G[jast + n * i];
        G[jast + n * i] = G[iter + n * i];
        G[iter + n * i] = a;
      }
    }

    G[iter*(n+1)]= sqrt(diagG[jast]);
    a=-.5/sigma/sigma;

    for (i=iter+1; i<=n-1; i++){
      if (m<=1){
        b=(x[pp[iter]]-x[pp[i]])*(x[pp[iter]]-x[pp[i]]);
      }else{
        b=0.0;
        for (j=0;j<=m-1;j++){
          c=x[j+m*pp[iter]]-x[j+m*pp[i]];
          b+=c*c;
        }
      }
      G[i+n*iter]=exp(a*b);
    }

    if (iter>0)
      for (j=0; j<=iter-1; j++)
        for (i=iter+1; i<=n-1; i++) G[i+n*iter]-=G[i+n*j]*G[iter+n*j];

    for (i=iter+1; i<=n-1; i++)
    {
      G[i+n*iter]/=G[iter*(n+1)];
    }

    residual=0.0;
    jast=iter+1;
    maxdiagG=0;
    for (i=iter+1; i<=n-1; i++)
    {
      b=1.0;
      for (j=0;j<=iter;j++)
      {
        b-=G[i+j*n]*G[i+j*n];
      }
      diagG[i]=b;
      if (b>maxdiagG)
      {
        jast=i;
        maxdiagG=b;
      }
      residual+=b;
    }

    iter++;
  }

  // collect objects G and Pvec
  Rcpp::NumericMatrix g_out(n, iter, G.begin());
  Rcpp::NumericMatrix pvec_out(1, n, pp.begin());

  // return list
  return Rcpp::List::create(
    Rcpp::Named("G") = g_out,
    Rcpp::Named("Pvec") = pvec_out
  );
}

// [[Rcpp::export]]
Rcpp::List contrast_ica(const arma::mat & X, double sigma, double eta, double kappa) {

  int i, N, m;
  arma::mat Rkappa;
  arma::mat sizes;
  Rcpp::NumericVector x;
  Rcpp::List res;

  N = X.n_cols;
  m = X.n_rows;

  // convert from armadillo to rcpp
  NumericMatrix mymat = NumericMatrix(m, N, X.begin());

  for (i = 0; i < m; ++i){

    // subset matrix
    x = mymat(i, _);

    // pass it to chol_gaussc
    res = chol_gaussc(x / sigma, 1, N * eta);

  }

  return Rcpp::List::create(
    Rcpp::Named("m") = m,
    Rcpp::Named("N") =  N,
    Rcpp::Named("v") = x,
    Rcpp::Named("res") = res,
    Rcpp::Named("i") = i
  );

  // return 2;
}
