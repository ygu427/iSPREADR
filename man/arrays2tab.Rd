\name{arrays2tab}
\alias{arrays2tab}
\alias{tab2arrays}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Table format to/from array format converters}

\description{ A list of 2D or 3D images \code{X1}, \code{X2}, ... can be
  stored in two popular formats.
  \itemize{
    \item{Array}{A list of multi-dimensional arrays such that
      \code{X1[i,j,k]} is the value of the corresponding pixel/voxel. The
      advantage of this format is that it is intuitive, enforces a
      "rectangular" shape of the figure, which is a prerequisite of the
      spatial kernel smoother used in this package.  The actual spatial
      location (coordinates) must be supplied separately, such as the
      \code{grids} list used by several other functions in this package.}
    \item{Table}{A 2D table with each row represents a pixel/voxel.  The
      first 2 or 3 columns store the spatial coordinates and each other
      column an image.  The main advantage of this format is that it more
      flexible.  The images don't have to be rectangular.  You can omit the
      background of brain images altogether, resulting in more compact
      object.  It is also easier to be exported to spreadsheet-like files.
      The \pkg{fields} package uses this format.}
  }
  These two functions provides convertion between these two formats.
}
\usage{
arrays2tab(grids, XYlist, zero.rm=FALSE, EPSILON=0.1^6)
tab2arrays(grids, tab, default.val=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grids}{A list of monotonically increasing coordinate grids.  For
    example, \code{grids <- list(grid.x=seq(0, 100, 2),
    grid.y=exp(seq(1, 5, .2))).}}

  \item{XYlist}{Lists of scans in array format.}
  \item{tab}{Scans in table format.}
  \item{zero.rm}{Should we remove voxels with constant zero
    values across all arrays?  Default to FALSE.}
  \item{EPSILON}{If the absolute value of an entry is less than EPSILON
    it will be reset to zero.  You can also use this to do simple
    background thresholding.}
  \item{default.val}{The value to be filled at blank voxels. Default to 0.}
}

\details{\code{arrays2tab()} and \code{tab2arrays()} are roughly the
  inverse of each other.  \code{arrays2tab()} is relative easy to
  implement and fast to run. but \code{tab2arrays()} is neither easy nor
  fast since it is hard to make an arbitrary shaped image without any
  particular ordering into a rectangular one.  \code{tab2arrays()}
  assumes the first few columns of the table are coordinates.
  \code{grids} still need to be specified manually because 1. otherwise
  there might be places that should be in the grid but not detectable
  automatically; 2. it provides a safe guard against some common
  mistakes.  If the first couple of columns have entirely different
  values (perhaps by a mistake), you may end up generating a 2M by 2M by
  2M array ---- it will eat all your memory and bring down your box in
  just seconds. }

\value{Either a table containing both spatial coordinates and scans or a
list of scans in array format.}

\author{Xing Qiu}
\seealso{\code{\link{rep.test}}}
\examples{
## Generate some FA maps in array format
N.x <- 8; N.y <- 8; N.z <- 3
for (i in 1:8) {
  FA.i <- array(rnorm(N.x*N.y*N.z), c(N.x,N.y,N.z))
  ## make some zeros and very small numbers
  FA.i[1:3, 1, ] <- 0; FA.i[1:3, 2, ] <- 0.00001
  assign(paste("FA",i,sep=""), FA.i)
}
grids <- list(2*(1:N.x)-1, 2*(1:N.y)-1, 2*(1:N.z)-1)
Xlist <- list(FA1, FA2, FA3, FA4); Ylist <- list(FA5, FA6, FA7, FA8)

tab1 <- arrays2tab(grids, Xlist) # dim(tab1) = c(N.x*N.y*N.z, 3+length(Xlist))
tab2 <- arrays2tab(grids, Xlist, zero.rm=TRUE) # 9 zeros are removed
# Now another 9 voxels with very small intensities are removed.
tab3 <- arrays2tab(grids, Xlist, zero.rm=TRUE, EPSILON=0.0001)

## Xlist2 is almost identical to Xlist except that those voxels with
## very small values are replaced by zeros
Xlist2 <- tab2arrays(grids, tab3) 
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}
