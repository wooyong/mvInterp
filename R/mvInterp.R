
mvInterpLinear = function(x, xGrid, yGrid) {
  mvInterpLinear_fast_cpp(x, xGrid, yGrid, ncol(xGrid))
}

subsetNearestPoints = function(x, xGrid) {
  subsetNearestPoints_cpp(x, xGrid, ncol(xGrid))
}

uniqueGridPoints = function(xGrid) {
  # add one since uniqueGridPoints_cpp returns zero-base labels
  uniqueGridPoints_cpp(xGrid, ncol(xGrid)) + 1
}
