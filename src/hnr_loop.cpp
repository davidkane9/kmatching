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
  // initialize Armadillo objects and other variables
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
  arma::vec center(y.begin(), y.size());

  
  
  double tmin; double tmax; int runs; int index = 0;
  double rdistance; int pickindex;
 
  
  for(int i = 0; i < n*skiplength + discard;i++) {
    
   tmin = 0; tmax = 0; runs = 0;
  
   Rcpp::RNGScope scope;
  
   while(tmin ==0 && tmax ==0) {
     
     if(achr && (index >= discard)) {
       // randomly select from the previous indices
       pickindex = (int) floor(runif(1)[0]*index);
       // u is unit vector from center to 'pickindex' point
       u = X.col(pickindex) - center;
       u = u/sqrt(dot(u,u));
       // note that this skips matrix mult, which dominated running time
     } else {
     
      // r is a random unit vector in with basis in Z
      r = rnorm(Z.n_cols);
      r = r/sqrt(sum(r*r));
      // copyless assignment to armadillo
      arma::vec r_arm(r.begin(), r.size(), false);
      // u is a direction vector in the appropriate k-plane pointing in a
      // random direction Z %*% r is the same as in mirror
      u = Z*r_arm;
     }
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

     // if stuck on boundary point
     if(tmin==0 && tmax ==0) {
       runs++;
       if(runs >= 1000) {
         throw std::runtime_error("hitandrun found can't find feasible direction, cannot generate points");
       }
     }

   }
   
   // distance to travel along segment
   rdistance =  (double) runif(1)[0];
   // new y is a random point on segment
   y = y + (tmin + (tmax - tmin)*rdistance)*u;

   if(i % skiplength == 0 || index< discard) {
     for(unsigned int k=0; k < X.n_rows; k++) {
       X(k,index) = y[k];
     }
     center = ((index+ 1) * center + y)/(index+2);
     index = index + 1;
   }
  }
   // return matrix of samples in columns
  arma::mat ret = X.submat( arma::span(0, X.n_rows - 1) , arma::span(discard, X.n_cols-1));
  return Rcpp::wrap(ret);
END_RCPP
}