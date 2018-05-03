\name{Smooth}
\alias{Smooth}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Multi-dimensional smoothing function.}
\description{
  Multi-dimensional smoothing function.
}
\usage{Smooth(grids, img, sm.method="ksmooth", ...)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grids}{A list of monotonically increasing coordinate grids.  For
    example,
    \code{grids <- list(grid.x=seq(0, 100, 2),
    grid.y=exp(seq(1, 5, .2))).}}
  \item{img}{An m-dim array defined on \code{grids} to be smoothed.}
  \item{sm.method}{Smoothing method. Valid choices: 1. \code{ksmooth} (kernel
    smoothing with a gauss kernel, the default choice), 2. \code{tps}
    (thin plate spline), 3. \code{anisodiff} (anisotropic smoothing
    based on Perona-Malik equation)}
  \item{...}{Other options passed to different smoothers such as
    \itemize{
      \item{\code{kernel}:}{ A choice of kernel (for kernel smoother;
        see help of function \code{ksmooth}) for more details.}
      \item{bandwidth:}{ For options ksmooth and Tps. The smoothing
        bandwidth. The kernels are scaled so that their 
        quartiles (viewed as probability densities) are at +/-
        '0.25*bandwidth'.}
      \item{niter:}{ For option \code{anisodiff} only. Number of iterations. }
      \item{kappa.param:}{ For option \code{anisodiff} only. A trade-off
        parameter that controls the estimation of Kappa. (the diffusion
        contrast parameter). It is expressed as 'lambda' in Equation (8)
      of the iSPREAD paper (Liu et al Phys. Med. Biol. 2016).}
      \item{lambda:}{ For option \code{anisodiff} only. This is 'lamdba'
        in Equation (7), Perona and Malik IEEE/TPAMI 1990.  This
        parameter  controls the speed of numerical scheme and must be in
        [0, 1/4] to ensure the stability of the algorithm.}
      \item{option:}{ For option \code{anisodiff} only. 1. Use Gaussian
        p.d.f. as the diffusivity function. 2. Use Charbonier
        conductivity function as the diffusivity function. The second
        one is the default option because it gives more stable results.}
      There are other options for function \code{Tps()} implemented in
      \pkg{fields}.}
  }
}
\details{
  This function is a wrapper of different smoothing methods.  The
  default method, "ksmooth", is a multi-dimensional generalization of
  \code{ksmooth()}.  The second method, "anisodiff", implements an
  anisotropic smoother based on Perona-Malik equation. Care must be
  taken on the boundary points because   currently there is no mechanism
  to take care of them automatically.  I may add some code to deal
  with the boundary points more gracefully in a future release.  Another
  method, "tps", is the thin plate spline smoother, depends on package
  \pkg{fields}.  This function is too slow for very large scale data, it
  is included primarily as an example. 
}
\value{An array of the same dimension as \code{img}.}
\references{[kernel smoothing, multi-dim]}
\author{Xing Qiu}
\seealso{\code{\link{ksmooth}}, \code{\link{Tps}}}
\examples{
grid.x <- seq(1, 255, 2); grid.y <- seq(1, 255, 2); grid.z <- seq(1,5,2)
grids <- list(grid.x, grid.y, grid.z)
Y <- array(rnorm(128*128*3), c(128,128,3))
myfit1 <- Smooth(grids, Y, bandwidth=5.0, kernel="normal")
myfit2 <- Smooth(grids, Y, bandwidth=5.0, kernel="box")
myfit3 <- Smooth(grids, Y, kappa.param=0.2, sm.method="anisodiff")
## The following command depends on package fields.
\dontrun{myfit3 <- Smooth(grids, Y, bandwidth=5.0, sm.method="tps")}
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
