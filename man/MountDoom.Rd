\name{MountDoom}
\alias{MountDoom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A crude algorithm to remove boundary voxels and reduce the
  ringing effects. }
\description{
  It takes an image arrary; removes boundary voxels within certain
  numbers of voxels.
}
\usage{
MountDoom(img, minval=0, thickness=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{An image presented as a 3D array.}
  \item{thickness}{The "thickness" of boundary, measured in number of
    voxels, to be removed.}
  \item{minval}{Any voxels with value less or equal to minval will be
    considered as the initial background.}
}
\details{
  Some text to be written.
}
\value{A mask file (0=background; 1=foreground) that can be used to remove boundary points. }
\author{Xing Qiu}
\examples{
mask1 <- MountDoom(img1)

## we can use multiplication to get the trimmed image
img1b <- mask1*img1

## A more aggressive boundary removal
mask1b <- MountDoom(img1, thickness=3)

}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
