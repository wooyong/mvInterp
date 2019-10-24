#ifndef MVINTERP_H
#define MVINTERP_H

#include <Rcpp.h>
using namespace Rcpp;

NumericVector which_mins_cpp(NumericVector x, int n);

NumericVector subsetEqual_cpp(NumericVector x, NumericMatrix xGrid, int dimension);

NumericVector subsetNearestPoints_cpp(NumericVector x, NumericMatrix xGrid, int dimension);

NumericVector uniqueGridPoints_cpp(NumericMatrix xGrid, int dimension);

double uvInterpConstant_cpp(double x, NumericVector xGrid, NumericVector yGrid, int lowerBound, int returnValue);

double mvInterpConstant_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension, int lowerBound);

double mvInterpConstant_fast_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension, int lowerBound);

double uvInterpLinear_cpp(double x, NumericVector xGrid, NumericVector yGrid);

double mvInterpLinear_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension);

double mvInterpLinear_fast_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension);

#endif