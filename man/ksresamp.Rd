\name{iSPREADR}
\alias{iSPREADR}
\docType{package}
\title{
  Spatial data analysis based on kernel smoothing and resampling
}
\description{
This package provides a multi-dimensional kernel smoother
and some tools for hypothesis testing based on resampling techniques.
Although this package is designed with diffusion tensor imaging data
as a "model data", it may be used to detect significantly different
regions for other types of multi-dimensional, large scale spatial
data.  Since this package is designed for analyzing large data (about
2M data points) by resampling method, computational efficiency has
very high priority.  The current implementation depends on the
implicit parallelization provided by foreach and doMC packages.
}
\details{
\tabular{ll}{
Package: \tab iSPREADR\cr
Type: \tab Package\cr
Version: \tab 0.7\cr
Date: \tab 2011-05-06\cr
License: \tab GPL (>=2.0)\cr
LazyLoad: \tab yes\cr
}
Please see the help files of individual functions.
}
\author{
  Xing Qiu
Maintainer: Xing Qiu <xqiu@bst.rochester.edu>
}
\references{
  Zhu, Hu, Qiu radiology paper, maybe another paper in statistical journal.
}
\keyword{package}
\examples{
# See examples of each function's help file
}
