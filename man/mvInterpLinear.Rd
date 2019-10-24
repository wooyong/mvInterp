\name{mvInterpLinear}
\alias{mvInterpLinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Linear Interpolation of Multivariate Functions
}
\description{
Linear interpolation of functions defined over arbitrary dimension.
}
\usage{
mvInterpLinear(x, xGrid, yGrid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of coordinates at which the interpolated value is computed.
}
  \item{xGrid}{
a matrix representing a rectangular grid over the domain of the function, each row representing an evaluation point.
}
  \item{yGrid}{
a vector of function values evaluated at each row of \code{xGrid}.
}
}
\details{
This function recursively applies bilinear interpolation to evaluate a function defined over arbitrary dimension. The grid must be a rectangular grid. If the evaluation point is out of grid, extrapolation is performed using boundary values.
}
\value{
the function value evaluated at \code{x} via linear interpolation.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Wooyong Lee
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
require(mvInterp)
xGrid = as.matrix(expand.grid(rep(list(c(1,2,3)),3)))
yGrid = xGrid[,1] - xGrid[,2] + xGrid[,3]^2
mvInterpLinear(c(  2, 2,   2), xGrid, yGrid) # exact evaluation
mvInterpLinear(c(1.5, 2, 2.5), xGrid, yGrid) # interpolation
mvInterpLinear(c(  4, 2,   4), xGrid, yGrid) # extrapolation
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
