// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
NumericMatrix ris_recursion(int r, int NT, arma::mat BCH, arma::mat rand_mat, arma::vec V) {

  mat ln_prob(NT, r, fill::zeros);
  for (int j = 0; j < r; ++j) {
    rowvec nu(NT, fill::zeros);
    for (int z = 0; z < NT; ++z) {
      int zz = NT - z - 1;
      double sumterm = 0;
      if (zz != NT - 1) {
        rowvec b = BCH.row(zz);
        sumterm = dot(nu, b);
      }
      double nu0 = (1/BCH.at(zz,zz))*(V.at(zz) - sumterm);
      ln_prob.at(zz,j) = R::pnorm5(nu0, 0, 1, 1, 1);
      nu.at(zz) = R::qnorm5(rand_mat.at(zz,j), 0, 1, 1, 0)*R::pnorm5(nu0, 0, 1, 1, 0);
    }
  }

  NumericMatrix out = wrap(ln_prob);
  return out;
}
