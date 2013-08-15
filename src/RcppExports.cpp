// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// hnr_loop
SEXP hnr_loop(SEXP y_s, SEXP Z_s, SEXP n_s, SEXP skiplength_s, SEXP discard_s, SEXP achr_s);
RcppExport SEXP kmatching_hnr_loop(SEXP y_sSEXP, SEXP Z_sSEXP, SEXP n_sSEXP, SEXP skiplength_sSEXP, SEXP discard_sSEXP, SEXP achr_sSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        SEXP y_s = Rcpp::as<SEXP >(y_sSEXP);
        SEXP Z_s = Rcpp::as<SEXP >(Z_sSEXP);
        SEXP n_s = Rcpp::as<SEXP >(n_sSEXP);
        SEXP skiplength_s = Rcpp::as<SEXP >(skiplength_sSEXP);
        SEXP discard_s = Rcpp::as<SEXP >(discard_sSEXP);
        SEXP achr_s = Rcpp::as<SEXP >(achr_sSEXP);
        SEXP __result = hnr_loop(y_s, Z_s, n_s, skiplength_s, discard_s, achr_s);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
