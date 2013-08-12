#include <Rcpp.h>
#include <RcppArmadillo.h>

// y: starting point
// Z: null-space matrix
// n: number of outputs
// skiplength: number to skip
// discard: burninlength
RccpExport SEXP hnr.loop(SEXP y_s, SEXP Z_s, SEXP n_s, SEXP skiplength_s, SEXP discard_s) {
BEGIN_RCCP
  Rcpp::NumericVector y(y_s);
  Rcpp::NumericMatrix Z(Z_s);
  // can't do optimized matrix mult with Rcpp objects, so do a copyless (fast)
  // assignment to armadillo mat object
  arma::mat yarm(y.begin(), y.size(), false);
  arma::mat Zarm(Z.begin(), Z.nrow(), Z.ncol(), false);
  int n = Rcpp::as<int>(n_s);
  int skiplength = Rcpp::as<int>(skiplength_s);
  int discard = Rcpp::as<int>(discard);

  Rcpp::RNGScope scope;

  while(tmin ==0 && tmax ==0) {
    // r is a random unit vector in with basis in Z
    r = rnorm(ncol(Z));
    r = r/sqrt(sum(r^2));
    // copyless assignment to armadillo
    arma::colvec r_arm(r.begin(), r.size(), false);
    // u is a unit vector in the appropriate k-plane pointing in a
    // random direction Z %*% r is the same as in mirror
    u = Zarm*r_arm;
    
    //unboundedness
    if(any(u==0)) {
      throw std::runtime_error("Problem is unbounded")
    }
    c = yarm/u;
    // determine intersections of x + t*u with walls
    // the limits on how far you can go backward and forward
    // i.e. the maximum and minimum ratio y_i/u_i for negative and positive u.
    tmin = max(-c[u>0]); tmax = min(-c[u<0]);

    // if stuck on boundary point
    if(tmin==0 && tmax ==0) {
      runs++;
      if(runs >= 1000) {
        throw std::runtime_error("hitandrun found can't find feasible direction, cannot generate points");
      }
    }
  }

  y <- y + (tmin + (tmax - tmin)*Rcpp::runif(1))*u;

  if(i % skiplength == 0) {
    X(,index) = y;
    index = index + 1;
  }
  // return matrix of samples in columns
  return Rcpp::wrap(X(,Range(discard, X.ncol()-1)));
END_RCPP
}