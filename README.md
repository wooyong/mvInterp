# mvInterp - multivariate linear interpolation

### Introduction

This package performs linear interpolation of multivariate functions of arbitrary dimension.

This is done by recursively applying bilinear interpolation.

Linear extrapolation is performed for evaluation points out of the grid.

### Installation

To install **mvInterp**, type

```r
install.packages("https://wooyong.github.io/data/packages/mvInterp_2019.10.23.tar.gz", repos=NULL, type="source")
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

A rectangular grid is required for interpolation.

A grid is represented by a matrix whose rows represent evaluation points.

The corresponding function values are stored in another vector.

The following example interpolates a function `f(x1, x2, x3) = x1 - x2 + x3^2`.

```r
require(mvInterp)

# create f(x1, x2, x3) over the grid [1,2,3]^3.
xGrid = as.matrix(expand.grid(rep(list(c(1,2,3)),3)))
yGrid = xGrid[,1] - xGrid[,2] + xGrid[,3]^2

# exact evaluation
mvInterpLinear(c(  2, 2,   2), xGrid, yGrid)

# interpolation
mvInterpLinear(c(1.5, 2, 2.5), xGrid, yGrid)

# extrapolation
mvInterpLinear(c(  4, 2,   4), xGrid, yGrid)
```
