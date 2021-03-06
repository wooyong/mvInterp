# mvInterp - multivariate linear interpolation

### Introduction

This package performs linear interpolation of multivariate functions of arbitrary dimension by recursively applying bilinear interpolation.

### Installation

To install **mvInterp**, type

```r
install.packages("https://wooyong.github.io/data/packages/mvInterp_1.0.0.tar.gz", repos=NULL, type="source")
```

Alternatively, to install **mvInterp** directly from source on Github, type

```r
# if devtools package is not installed, install it
install.packages("devtools")

# install mvInterp from Github using devtools command
library(devtools)
install_github("wooyong/mvInterp")
```

To load package, simply type

```r
library(mvInterp)
```

### Usage

The following example interpolates a function `f(x1, x2, x3) = x1 - x2 + x3^2` of three variables, by recursively applying bilinear interpolation.

A rectangular grid is required for `mvInterp`. Linear extrapolation is performed outside of the grid.

```r
require(mvInterp)

# create a rectangular grid [1,2,3]^3.
xGrid = as.matrix(expand.grid(rep(list(c(1,2,3)),3)))

# evaluate function values over the grid
yGrid = xGrid[,1] - xGrid[,2] + xGrid[,3]^2

# exact evaluation
mvInterpLinear(c(  2, 2,   2), xGrid, yGrid)
[1] 4

# linear interpolation
mvInterpLinear(c(1.5, 2, 2.5), xGrid, yGrid)
[1] 6

# linear extrapolation
mvInterpLinear(c(  4, 2,   4), xGrid, yGrid)
[1] 16
```
