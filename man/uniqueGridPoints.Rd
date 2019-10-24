\name{uniqueGridPoints}
\alias{uniqueGridPoints}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Indices for Unique Grid Points
}
\description{
Create a vector of indices that are unique evaluation points in the grid, which can be used to remove duplicate points in the grid.
}
\usage{
uniqueGridPoints(xGrid)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xGrid}{
a matrix representing a grid over the domain of the function, each row representing an evaluation point.
}
}
\details{
This function creates a vector of indices that correspond to unique evaluation points in the grid. The grid need not be a rectangular grid and this function can be applied to any matrix to retrieve unique row indices.
}
\value{
a vector of indices that correspond to unique row entries of \code{xGrid}.
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
xGrid = rbind(xGrid, xGrid)           # make duplicates
xGrid = xGrid[sample(1:nrow(xGrid)),] # sort rows
uniqueGridPoints(xGrid)               # returns row indices of unique entries
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
