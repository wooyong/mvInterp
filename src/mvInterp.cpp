#include "mvInterp.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector which_mins_cpp(NumericVector x, int n) {

    // prepare auxiliary variables
    NumericVector index(1), indices(n);
    NumericVector y = clone(x);

    // sort grid up to nth minimum
    std::nth_element(y.begin(), y.begin()+n-1, y.end());
    std::sort(y.begin(), y.begin()+n-1);

    // record indices
    for(int k=0; k<n; k++) {
        index = match(NumericVector::create(y(k)), x) - 1;
        indices(k) = index(0);
    }

    // return indices array
    return indices;
}

// [[Rcpp::export]]
NumericVector subsetEqual_cpp(NumericVector x, NumericMatrix xGrid, int dimension) {

    // prepare labels array
    NumericVector labels;

    // prepare auxiliary variables
    int count;
    int d;

    // add label if equal
    for(int k=0; k<xGrid.nrow(); k++) {
        count = 0;
        for(d=0; d<dimension; d++) {
            if(x(d) == xGrid(k,d)) {
                count++;
            }
        }
        if(count == d) {
            labels.push_back(k);
        }
    }

    // return result
    return labels;
}

// [[Rcpp::export]]
NumericVector subsetNearestPoints_cpp(NumericVector x, NumericMatrix xGrid, int dimension) {

    // read dimension
    int gridLength = xGrid.nrow();

    // prepare counter
    NumericVector count(gridLength);

    // prepare auxiliary variables
    NumericVector coordinate(gridLength);
    NumericVector coordinateSubset;
    NumericVector countTicker(gridLength);
    double lb, ub;

    // pick nearest points at each dimension
    for(int d=0; d<dimension; d++) {

        // read coordinate grid
        coordinate = unique(xGrid(_,d));

        // if out of minimum, pick up the first two smallest values
        if(x(d) < min(coordinate)) {
            NumericVector indices = which_mins_cpp(coordinate, 2);
            lb = coordinate(indices(0));
            ub = coordinate(indices(1));

        // if out of maximum, pick up the first two largest values
        } else if(max(coordinate) < x(d)) {
            NumericVector indices = which_mins_cpp((-1) * coordinate, 2); // negate vector to search for max using min
            ub = coordinate(indices(0));
            lb = coordinate(indices(1));

        // if in the grid, pick up the largest lower bound and the smallest upper bound
        } else {
            coordinateSubset = coordinate[coordinate <= x(d)];
            lb = max(coordinateSubset);
            coordinateSubset = coordinate[x(d) <= coordinate];
            ub = min(coordinateSubset);
        }

        // increase count if lb or ub
        countTicker = (lb <= xGrid(_,d)) * (xGrid(_,d) <= ub);
        count = count + countTicker;
    }

    // prepare labels array
    NumericVector labels;

    // add label if equal
    for(int k=0; k<xGrid.nrow(); k++) {
        if(count(k) == dimension) {
            labels.push_back(k);
        }
    }

    // return result
    return labels;
}

// [[Rcpp::export]]
NumericVector uniqueGridPoints_cpp(NumericMatrix xGrid, int dimension) {

    // read dimension
    int gridLength = xGrid.nrow();

    // initialize unique entries indicator
    NumericVector uniqueEntriesIndicator = rep(1.0, gridLength);
    NumericVector uniqueEntries;

    // define auxiliary variables
    int j, k, l, count;

    // pick up unique entries
    for(j=0; j<gridLength; j++) {
        // if indicator is on, add to uniqueEntries and turn off indicators of all duplicates
        if(uniqueEntriesIndicator(j) > 0) {
            // add to uniqueEntries
            uniqueEntries.push_back(j);
            // turn off indicators of all duplicates
            for(k=j+1; k<gridLength; k++) {
                count = 0;
                for(l=0; l<dimension; l++) {
                    if(xGrid(j,l) == xGrid(k,l)) {
                        count++;
                    }
                }
                if(count == dimension) {
                    uniqueEntriesIndicator(k) = 0;
                }
            }
        }
    }

    // return unique entry indices
    return uniqueEntries;
}

// [[Rcpp::export]]
double uvInterpConstant_cpp(double x, NumericVector xGrid, NumericVector yGrid, int lowerBound, int returnValue) {

    // initialize index of the nearest point
    int boundIndex;

    if(lowerBound == 1) {

        // if lowerBound = 1, pick the largest lower bound
        boundIndex = which_min(xGrid);
        for(int k=0; k<yGrid.length(); k++) {
            if(xGrid(boundIndex) < xGrid(k) && xGrid(k) <= x) {
                boundIndex = k;
            }
        }

    } else {

        // if lowerBound = 0, pick the smallest upper bound
        boundIndex = which_max(xGrid);
        for(int k=0; k<yGrid.length(); k++) {
            if(x <= xGrid(k) && xGrid(k) < xGrid(boundIndex)) {
                boundIndex = k;
            }
        }

    }

    // return value
    if(returnValue == 1) {
        return yGrid(boundIndex);
    } else {
        return boundIndex;
    }

}

// [[Rcpp::export]]
double mvInterpConstant_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension, int lowerBound) {

    // if dimension is 1, return nearest value
    if(dimension == 1) {

        return uvInterpConstant_cpp(x(0), xGrid(_,0), yGrid, lowerBound, 1);

    // if dimension is not 1, interpolate using only the first dimension and call the function recursively for the rest of the dimensions
    } else {

        // read length of grid
        int gridLength = yGrid.length();

        // create dimension index range excluding the first dimension
        Range dimensionRange(1, dimension-1);

        // prepare interpolated grid
        NumericVector yGridNew(gridLength);

        // prepare auxiliary variables
        NumericVector xCoordinate(dimension-1);
        NumericVector xGridFirstDimensions = xGrid(_,0);

        for(int k=0; k<gridLength; k++) {

            // read k-th coordinate
            xCoordinate = xGrid(k,_);

            // read grid points that have the same values in the dimensions other than the first dimension
            NumericVector equalLabels =  subsetEqual_cpp(xCoordinate[dimensionRange], xGrid(_,dimensionRange), dimension-1);

            // interpolate according to the first dimension using points that are equal for the rest of the dimension
            yGridNew(k) = uvInterpConstant_cpp(x(0), xGridFirstDimensions[equalLabels], yGrid[equalLabels], lowerBound, 1);

        }

        // recursively call the function with removing the first dimension
        return mvInterpConstant_cpp(x[dimensionRange], xGrid(_,dimensionRange), yGridNew, dimension-1, lowerBound);

    }
}

// [[Rcpp::export]]
double mvInterpConstant_fast_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension, int lowerBound) {

    // read labels
    NumericVector nearestPoints = subsetNearestPoints_cpp(x, xGrid, dimension);

    // create subset x and y grids
    int subsetLength = nearestPoints.length();
    NumericMatrix xGridSubset(subsetLength, dimension);
    NumericVector yGridSubset(subsetLength);

    // fill in the subset grids
    for(int k=0; k<subsetLength; k++) {
        xGridSubset(k,_) = xGrid(nearestPoints(k),_);
        yGridSubset(k) = yGrid(nearestPoints(k));
    }

    // call interpolation function
    return mvInterpConstant_cpp(x, xGridSubset, yGridSubset, dimension, lowerBound);
}

// [[Rcpp::export]]
double uvInterpLinear_cpp(double x, NumericVector xGrid, NumericVector yGrid) {

    // declare lb and ub indices
    int lbIndex, ubIndex;

    // if out of minimum, set lb and ub to be the first two smallest values
    if(x < min(xGrid)) {
        NumericVector indices = which_mins_cpp(xGrid, 2);
        lbIndex = indices(0);
        ubIndex = indices(1);

    // if out of maximum, set lb and ub to be the first two largest values
    } else if(x > max(xGrid)) {
        NumericVector indices = which_mins_cpp((-1) * xGrid, 2); // negate vector to search for max using min
        lbIndex = indices(0);
        ubIndex = indices(1);

    // if in the grid, set lb and ub to be nearest points
    } else {
        // compute index of nearest points
        lbIndex = uvInterpConstant_cpp(x, xGrid, yGrid, 1, 0);
        ubIndex = uvInterpConstant_cpp(x, xGrid, yGrid, 0, 0);
    }

    // compute linear interpolation. If out of bound, return boundary value
    if(lbIndex == ubIndex) {
        return yGrid(lbIndex);
    } else {
        double weight = (x - xGrid(lbIndex)) / (xGrid(ubIndex) - xGrid(lbIndex));
        return (1-weight) * yGrid(lbIndex) + weight * yGrid(ubIndex);
    }
}

// [[Rcpp::export]]
double mvInterpLinear_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension) {

    // if dimension = 1, return linear interpolation
    if(dimension == 1) {

        // retrieve unique entries and apply linear interpolation
        NumericVector xGridVector = xGrid(_,0);
        NumericVector uniqueEntries = uniqueGridPoints_cpp(xGrid, 1);
        return uvInterpLinear_cpp(x(0), xGridVector[uniqueEntries], yGrid[uniqueEntries]);

    // if dimension is not 1, interpolate using only the first dimension and call the function recursively for the rest of the dimensions
    } else {

        // read length of grid
        int gridLength = yGrid.length();

        // create dimension index range excluding the first dimension
        Range dimensionRange(1, dimension-1);

        // prepare interpolated grid
        NumericVector yGridNew(gridLength);

        // prepare auxiliary variables
        NumericVector xCoordinate(dimension-1);
        NumericVector xGridFirstDimensions = xGrid(_,0);

        for(int k=0; k<gridLength; k++) {

            // read k-th coordinate
            xCoordinate = xGrid(k,_);

            // read grid points that have the same values in the dimensions other than the first dimension
            NumericVector equalLabels =  subsetEqual_cpp(xCoordinate[dimensionRange], xGrid(_,dimensionRange), dimension-1);

            // interpolate according to the first dimension using points that are equal for the rest of the dimension
            yGridNew(k) = uvInterpLinear_cpp(x(0), xGridFirstDimensions[equalLabels], yGrid[equalLabels]);

        }

        // recursively call the function with removing the first dimension
        return mvInterpLinear_cpp(x[dimensionRange], xGrid(_,dimensionRange), yGridNew, dimension-1);

    }
}

// [[Rcpp::export]]
double mvInterpLinear_fast_cpp(NumericVector x, NumericMatrix xGrid, NumericVector yGrid, int dimension) {

    // read labels
    NumericVector nearestPoints = subsetNearestPoints_cpp(x, xGrid, dimension);

    // create subset x and y grids
    int subsetLength = nearestPoints.length();
    NumericMatrix xGridSubset(subsetLength, dimension);
    NumericVector yGridSubset(subsetLength);

    // fill in the subset grids
    for(int k=0; k<subsetLength; k++) {
        xGridSubset(k,_) = xGrid(nearestPoints(k),_);
        yGridSubset(k) = yGrid(nearestPoints(k));
    }

    // call interpolation function
    return mvInterpLinear_cpp(x, xGridSubset, yGridSubset, dimension);
}