#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>

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
Rcpp::List chol_gaussc(const Rcpp::NumericVector & x, double sigma, double tol) {
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

// [[Rcpp::export]]
Rcpp::List symm_eigen(const arma::mat & X) {
  int N;
  arma::mat U, V_temp, V;
  arma::vec S_temp, S;
  arma::vec ids;

  N = X.n_rows;
  arma::svd(U, S_temp, V_temp, X, "standard");
  V = arma::fliplr(V_temp);
  S = arma::flipud(S_temp);

  return Rcpp::List::create(
    Rcpp::Named("D") = S,
    Rcpp::Named("A") = V
  );
}

// [[Rcpp::export]]
Rcpp::List which2(Rcpp::NumericVector x, double lb) {
  Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
  Rcpp::NumericVector x_temp = x[!Rcpp::is_na(x)];
  Rcpp::IntegerVector v_temp = v[!Rcpp::is_na(x)];

  return Rcpp::List::create(
    Rcpp::Named("x") = x_temp[x_temp > lb],
    Rcpp::Named("id") = v_temp[x_temp > lb]
  );
}

// [[Rcpp::export]]
Rcpp::List contrast_icac(const arma::mat & X, double sigma, double eta, double kappa) {

  int i, N, m;
  arma::mat Rkappa;
  arma::mat sizes;
  Rcpp::NumericVector x;
  Rcpp::List res;

  N = X.n_cols;
  m = X.n_rows;

  // convert from armadillo to rcpp
  Rcpp::NumericMatrix mymat = Rcpp::NumericMatrix(m, N, X.begin());
  Rcpp::NumericVector D_filtered;

  for (i = 0; i < m; ++i){

    // subset matrix
    x = mymat(i, Rcpp::_);

    // pass it to chol_gaussc
    res = chol_gaussc(x / sigma, 1, N * eta);

    // take out values
    arma::mat G_temp = Rcpp::as<arma::mat>(res["G"]);
    arma::vec Pvec_temp = Rcpp::as<arma::vec>(res["Pvec"]);
    arma::uvec Pvec = arma::sort_index(Pvec_temp);

    // center cols of G
    arma::mat G = center_cols(G_temp);

    // regularization
    Rcpp::List res_eigen = symm_eigen(G.t() * G);
    Rcpp::NumericVector D = Rcpp::as<Rcpp::NumericVector>(res_eigen["D"]);
    Rcpp::List res_which = which2(D, .1);

    D_filtered = res_which["x"];
    Rcpp::IntegerVector indexes = res_which["id"];
    // !!!

    Rcpp::Rcout << D_filtered << "\n";
    Rcpp::Rcout << indexes << "\n";

  }

  return Rcpp::List::create(
    Rcpp::Named("m") = m,
    Rcpp::Named("N") =  N,
    Rcpp::Named("D_filtered") = D_filtered,
    Rcpp::Named("res") = res,
    Rcpp::Named("i") = i
  );

  // return 2;
}
