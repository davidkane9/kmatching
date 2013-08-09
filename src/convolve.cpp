#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int fibonacci(const int x) {
  if(x<2)
    return x;
  else 
    return fibonacci(x-1) + fibonacci(x-2);
}

// [[Rcpp::export]]
double sumC(NumericVector x) {
    int n = x.size();
    double total = 0;
    for(int i = 0; i < n; ++i) {
      total += x[i];
    }
    return total;
}