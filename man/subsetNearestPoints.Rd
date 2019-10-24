\name{subsetNearestPoints}
\alias{subsetNearestPoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset Grid Points Near a Given Point
}
\description{
Subset grid points that are relevant for linear interpolation, namely grid points that are nearest to a given point.
}
\usage{
subsetNearestPoints(x, xGrid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a vector of coordinates.
}
  \item{xGrid}{
a matrix representing a rectangular grid over the domain of the function, each row representing an evaluation point.
}
}
\details{
This function picks up two nearest points (one if exactly on \code{xGrid}) for each coordinate of \code{x}, which are the points used for linear interpolation. If \code{x} is out of grid, two smallest or largest coordinates are selected for extrapolation.
}
\value{
a subset of \code{xGrid} containing nearest points to \code{x}.
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

# exact point in the grid is subsetted
subsetNearestPoints(c(  2, 2,   2), xGrid)

# two nearest points for first and third dimension are subsetted
subsetNearestPoints(c(1.5, 2, 2.5), xGrid)

# two largest points for first and third dimension are subsetted
subsetNearestPoints(c(  4, 2,   4), xGrid)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
