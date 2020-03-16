// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// speciation_r_cpp
double speciation_r_cpp(const NumericVector& tm, const NumericMatrix& tree, const NumericVector& pars, double soc);
RcppExport SEXP _emphasis_speciation_r_cpp(SEXP tmSEXP, SEXP treeSEXP, SEXP parsSEXP, SEXP socSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type tm(tmSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< double >::type soc(socSEXP);
    rcpp_result_gen = Rcpp::wrap(speciation_r_cpp(tm, tree, pars, soc));
    return rcpp_result_gen;
END_RCPP
}
// sum_speciation_rate
double sum_speciation_rate(double cbt, const NumericMatrix& tree, const NumericVector& pars, double soc);
RcppExport SEXP _emphasis_sum_speciation_rate(SEXP cbtSEXP, SEXP treeSEXP, SEXP parsSEXP, SEXP socSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type cbt(cbtSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type tree(treeSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< double >::type soc(socSEXP);
    rcpp_result_gen = Rcpp::wrap(sum_speciation_rate(cbt, tree, pars, soc));
    return rcpp_result_gen;
END_RCPP
}
// augment_cpp
List augment_cpp(NumericVector brts_in, NumericVector pars, int soc);
RcppExport SEXP _emphasis_augment_cpp(SEXP brts_inSEXP, SEXP parsSEXP, SEXP socSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type brts_in(brts_inSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pars(parsSEXP);
    Rcpp::traits::input_parameter< int >::type soc(socSEXP);
    rcpp_result_gen = Rcpp::wrap(augment_cpp(brts_in, pars, soc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_emphasis_speciation_r_cpp", (DL_FUNC) &_emphasis_speciation_r_cpp, 4},
    {"_emphasis_sum_speciation_rate", (DL_FUNC) &_emphasis_sum_speciation_rate, 4},
    {"_emphasis_augment_cpp", (DL_FUNC) &_emphasis_augment_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_emphasis(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
