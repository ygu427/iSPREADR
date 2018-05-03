\name{artifacts}
\alias{artifacts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate some artifacts for a multi-dimensional image-like data}
\description{
  Generate multivariate Gaussian artifacts for a multi-dimensional
  image-like data (array object).  
}
\usage{
artifacts(grids, center=rep(0, length(grids)), std=0.3, radius=1.0, magnitude=1.0, xscale=diag(length(grids)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grids}{A list of monotonically increasing coordinate grids.  For
    example, \code{grids <- list(grid.x=seq(0, 100, 2),
    grid.y=exp(seq(1, 5, .2))).}}
  \item{center}{The center of the artifact.  Default value: the
    origin. Example: \code{center=c(3,5)} defines a center for 2d data.}
  \item{std}{The standard deviation of the (unscaled) Gaussian
    density is set to \code{std*radius}. Default value: 0.5.}
  \item{radius}{The radius of the (unscaled) artifact. Generated value
    outside of this radius is zero. Default to 10.}
  \item{maxval}{Generated Gaussian artifacts will be scaled so that the
    maxium value equals \code{maxval}. Default to 1.0}
  \item{minval}{Generated Gaussian artifacts will be set to zero if they
    are less than\code{minval}. Default to 0.2.}
  \item{xscale}{A scaling matrix for the shape of the artifact. Use this to
    generate elliptic artifacts. Default to the identity matrix.}
}
\details{
  The unscaled artifacts is a multivariate Gaussian density centered at
  \code{center} with standard deviation \code{std}.  Its support is a
  ball of radius \code{radius}.  A user specified \code{xscale} matrix
  can be used to stretch/rotate this ball.
}
\value{An array defined on \code{grids} which represents an artifact. It
  can be added to an existing array with the same grids.}
\author{Xing Qiu}
\examples{
## A3 is a 3x3 scaling matrix.
A3 <- matrix(c(1.5, 0.5, 0.0,
               0.5, 1.5, 0.0,
               0.0, 0.0, 1.0), nrow=3)
## From a DTI related simulation, 128x128x3 voxels, each 2x2x2mm large.
grid.x <- seq(1, 255, 2); grid.y <- seq(1, 255, 2); grid.z <- seq(1,5,2)
N.x <- length(grid.x); N.y <- length(grid.y); N.z <- length(grid.z)
grids <- list(grid.x, grid.y, grid.z)

##Sharp effects,  60 voxels, Max: 10. 
Z1 <- artifacts(grids, center=c(50,80,3), xscale=A3, radius=4, maxval=10.0, minval=3) 

## Diffused effects, 516 voxels. Max: 1.2.
Z2 <- artifacts(grids, center=c(50,80,3), maxval=1.2, radius=15)

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
