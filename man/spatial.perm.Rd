\name{spatial.perm}
\alias{spatial.perm}
\alias{flip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random permutation of the scan-label of each voxel while keeping the spatial information intact.}
\description{
  \code{spatial.perm()} takes a list of scans (arrays), randomly permute
  the scan labels of each voxel while keeping the spatial information
  intact.  This function is used in \code{rep.test()}.  \code{flip()} is
  a shortcut for pre/post analysis.  \code{diff.permed <- flip(pre - post)} is
  equivalent to \code{XY <- spatial.perm(pre, post); diff.permed <-
    XY[[1]] - XY[[2]]}.
}
\usage{
spatial.perm(XYlist)

flip(diff.array)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{XYlist}{Lists of multivariate samples to be spatially permuted.
    As an example, it can be a list of FA maps sampled from both
    groups. The elements in this list must be either R vectors or
    arrays. }
  \item{diff.array}{An array/vector of pre/post differences.}
}
\details{
  Spatial permutations.
}
\value{A list of N-distances using different kernels (\code{Ndist()}) or
    an K by N dimensional matrix of N-distances (\code{Ndist.perm()}),
    where each row represents one permutation and each column a kernel.}
\author{Xing Qiu}
\seealso{\code{\link{sample}}, \code{\link{combn}}}
\examples{
## Generate some FA maps
N.x <- 64; N.y <- 64; N.z <- 3
for (i in 1:8) assign(paste("FA",i,sep=""), array(rnorm(N.x*N.y*N.z), c(N.x,N.y,N.z)))
Xlist <- list(FA1, FA2, FA3, FA4); Ylist <- list(FA5, FA6, FA7, FA8)

diff.orig <- FA2 -FA1 
diff.permed <- flip(diff.orig)

FAs.permed <- spatial.perm(c(Xlist, Ylist))
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
