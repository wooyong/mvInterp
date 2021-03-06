// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// which_mins_cpp
NumericVector which_mins_cpp(NumericVector x, int n);
RcppExport SEXP _mvInterp_which_mins_cpp(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(which_mins_cpp(x, n));
    return rcpp_result_gen;
END_RCPP
}
// subsetEqual_cpp
NumericVector subsetEqual_cpp(NumericVector x, NumericMatrix xGrid, int dimension);
RcppExport SEXP _mvInterp_subsetEqual_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(subsetEqual_cpp(x, xGrid, dimension));
    return rcpp_result_gen;
END_RCPP
}
// subsetNearestPoints_cpp
NumericVector subsetNearestPoints_cpp(NumericVector x, NumericMatrix xGrid, int dimension);
RcppExport SEXP _mvInterp_subsetNearestPoints_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(subsetNearestPoints_cpp(x, xGrid, dimension));
    return rcpp_result_gen;
END_RCPP
}
// uniqueGridPoints_cpp
NumericVector uniqueGridPoints_cpp(NumericMatrix xGrid, int dimension);
RcppExport SEXP _mvInterp_uniqueGridPoints_cpp(SEXP xGridSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(uniqueGridPoints_cpp(xGrid, dimension));
    return rcpp_result_gen;
END_RCPP
}
// uvInterpConstant_cpp
double uvInterpConstant_cpp(double x, NumericVector xGrid, NumericVector yGrid, int lowerBound, int returnValue);
RcppExport SEXP _mvInterp_uvInterpConstant_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP, SEXP lowerBoundSEXP, SEXP returnValueSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yGrid(yGridSEXP);
    Rcpp::traits::input_parameter< int >::type lowerBound(lowerBoundSEXP);
    Rcpp::traits::input_parameter< int >::type returnValue(returnValueSEXP);
    rcpp_result_gen = Rcpp::wrap(uvInterpConstant_cpp(x, xGrid, yGrid, lowerBound, returnValue));
    return rcpp_result_gen;
END_RCPP
}
// mvInterpConstant_cpp
double mvInterpConstant_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension, int lowerBound);
RcppExport SEXP _mvInterp_mvInterpConstant_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP, SEXP dimensionSEXP, SEXP lowerBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yGrid(yGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    Rcpp::traits::input_parameter< int >::type lowerBound(lowerBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(mvInterpConstant_cpp(x, xGrid, yGrid, dimension, lowerBound));
    return rcpp_result_gen;
END_RCPP
}
// mvInterpConstant_fast_cpp
double mvInterpConstant_fast_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension, int lowerBound);
RcppExport SEXP _mvInterp_mvInterpConstant_fast_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP, SEXP dimensionSEXP, SEXP lowerBoundSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yGrid(yGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    Rcpp::traits::input_parameter< int >::type lowerBound(lowerBoundSEXP);
    rcpp_result_gen = Rcpp::wrap(mvInterpConstant_fast_cpp(x, xGrid, yGrid, dimension, lowerBound));
    return rcpp_result_gen;
END_RCPP
}
// uvInterpLinear_cpp
double uvInterpLinear_cpp(double x, NumericVector xGrid, NumericVector yGrid);
RcppExport SEXP _mvInterp_uvInterpLinear_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yGrid(yGridSEXP);
    rcpp_result_gen = Rcpp::wrap(uvInterpLinear_cpp(x, xGrid, yGrid));
    return rcpp_result_gen;
END_RCPP
}
// mvInterpLinear_cpp
double mvInterpLinear_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension);
RcppExport SEXP _mvInterp_mvInterpLinear_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yGrid(yGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(mvInterpLinear_cpp(x, xGrid, yGrid, dimension));
    return rcpp_result_gen;
END_RCPP
}
// mvInterpLinear_fast_cpp
double mvInterpLinear_fast_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension);
RcppExport SEXP _mvInterp_mvInterpLinear_fast_cpp(SEXP xSEXP, SEXP xGridSEXP, SEXP yGridSEXP, SEXP dimensionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xGrid(xGridSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type yGrid(yGridSEXP);
    Rcpp::traits::input_parameter< int >::type dimension(dimensionSEXP);
    rcpp_result_gen = Rcpp::wrap(mvInterpLinear_fast_cpp(x, xGrid, yGrid, dimension));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mvInterp_which_mins_cpp", (DL_FUNC) &_mvInterp_which_mins_cpp, 2},
    {"_mvInterp_subsetEqual_cpp", (DL_FUNC) &_mvInterp_subsetEqual_cpp, 3},
    {"_mvInterp_subsetNearestPoints_cpp", (DL_FUNC) &_mvInterp_subsetNearestPoints_cpp, 3},
    {"_mvInterp_uniqueGridPoints_cpp", (DL_FUNC) &_mvInterp_uniqueGridPoints_cpp, 2},
    {"_mvInterp_uvInterpConstant_cpp", (DL_FUNC) &_mvInterp_uvInterpConstant_cpp, 5},
    {"_mvInterp_mvInterpConstant_cpp", (DL_FUNC) &_mvInterp_mvInterpConstant_cpp, 5},
    {"_mvInterp_mvInterpConstant_fast_cpp", (DL_FUNC) &_mvInterp_mvInterpConstant_fast_cpp, 5},
    {"_mvInterp_uvInterpLinear_cpp", (DL_FUNC) &_mvInterp_uvInterpLinear_cpp, 3},
    {"_mvInterp_mvInterpLinear_cpp", (DL_FUNC) &_mvInterp_mvInterpLinear_cpp, 4},
    {"_mvInterp_mvInterpLinear_fast_cpp", (DL_FUNC) &_mvInterp_mvInterpLinear_fast_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvInterp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
