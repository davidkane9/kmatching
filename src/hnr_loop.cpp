//' Loop for hitandrun
//' 
//' Runs the hitandrun algorithm in a loop
//' @name hnr_loop
//' @param y_s starting point
//' @param Z_s null-space matrix
//' @param n_s number of outputs
//' @param skiplength_s number to skip
//' @param discard_s burninlength, will also be used to generate samples
//' @param achr_s "accelerated convergence hit and run"
//' 
//' @export
//#include "hnr_loop.h"
#include <RcppArmadillo.h>

using namespace Rcpp;
// [[Rcpp:interfaces(r,cpp)]]
// [[Rcpp::export]]
SEXP hnr_loop(SEXP y_s, SEXP Z_s, SEXP n_s, SEXP skiplength_s, SEXP discard_s, SEXP achr_s) {
BEGIN_RCPP
  arma::vec y = as<arma::vec>(y_s);
  arma::mat Z = as<arma::mat>(Z_s);
  int n = Rcpp::as<int>(n_s);
  int skiplength = Rcpp::as<int>(skiplength_s);
  int discard = Rcpp::as<int>(discard_s);
  bool achr = Rcpp::as<bool>(achr_s);
  Rcpp::NumericVector r(Z.n_cols);
  arma::vec u(Z.n_rows);
  arma::vec c(Z.n_rows);
  arma::mat X(Z.n_rows, n+discard);
  arma::vec spans(discard);
  
  
  double tmin; double tmax; int runs; int index = 0;
  double rdistance; bool accepted; double univar;
  int spanindex;
  
  for(int i = 0; i < n*skiplength + discard;i++) {
    
   tmin = 0; tmax = 0; runs = 0; accepted = false; spanindex =0;
  
   Rcpp::RNGScope scope;
  
   while((tmin ==0 && tmax ==0) || !accepted) {
     // r is a random unit vector in with basis in Z
     r = rnorm(Z.n_cols);
     r = r/sqrt(sum(r*r));
     // copyless assignment to armadillo
     arma::vec r_arm(r.begin(), r.size(), false);
     // u is a unit vector in the appropriate k-plane pointing in a
     // random direction Z %*% r is the same as in mirror
     u = Z*r_arm;
    
     //unboundedness
     if(sum(u==0) >=1) {
       throw std::runtime_error("Problem is unbounded");
     }
     c = y/u;
     // determine intersections of x + t*u with walls
     // the limits on how far you can go backward and forward
     // i.e. the maximum and minimum ratio y_i/u_i for negative and positive u.
    
     tmin = max(-c.elem(find(u>0))); 
    
     tmax = min(-c.elem(find(u<0)));
     
     spans[spanindex] = tmax - tmin;

     // if stuck on boundary point
     if(tmin==0 && tmax ==0) {
       runs++;
       if(runs >= 1000) {
         throw std::runtime_error("hitandrun found can't find feasible direction, cannot generate points");
       }
     }
     
     if(!achr) {
       accepted = true;
     } else {
       if(index < discard) {
         accepted = true;
       } else {
         // compare this point to the distribution of spans
         // if a random uniform variate is < quantile of span, accept
         univar = runif(1)[1];
         if(sum(spans[spanindex] >= spans)/discard > univar) {
           accepted = true;
         }
       }
     }
     spanindex = (spanindex + 1) % discard;
   }
   rdistance =  (double) runif(1)[0];
   y = y + (tmin + (tmax - tmin)*rdistance)*u;

   if(i % skiplength == 0 || index< discard) {
     for(unsigned int k=0; k < X.n_rows; k++) {
       X(k,index) = y[k];
     }
     index = index + 1;
   }
  }
   // return matrix of samples in columns
  arma::mat ret = X.submat( arma::span(0, X.n_rows - 1) , arma::span(discard, X.n_cols-1));
  return Rcpp::wrap(ret);
END_RCPP
}