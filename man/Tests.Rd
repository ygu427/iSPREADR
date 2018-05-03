\name{Tests}
\alias{pre.post.test}
\alias{rep.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Detecting abnormal region from either one pre and one post images
  or several repeated pre and post images.}
\description{
  These two functions provides user level interface to conduct
  hypothesis testing based on spatial smoothing.  Briefly, this
  procedure can be described as:
  1. kernel smoothing; 2. compute the mean diff image; 3. compute the
  original N-stat; 4. use permutation to generate null distributions;
  5. global hypothesis testing based on N-stats; 6. local hypothesis
  testing based on mean diff image, both FWER and FDR adjusted
  p-values. Please see the pre-print for more details.
}
\usage{
pre.post.test(grids, pre, post, perms=1000, balanced=TRUE,
norm=c("L1","L2","Linf"), method=c("null", "normal"), MTP="BH", ...)

rep.test(grids, pre.list, post.list, perms=50,
rand.comb=0, balanced=TRUE, norm=c("L1","L2","Linf"),
method=c("null", "normal"), MTP="BH", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grids}{A list of monotonically increasing coordinate grids.  For
    example, \code{grids <- list(grid.x=seq(0, 100, 2),
    grid.y=exp(seq(1, 5, .2))).}}
  \item{pre}{An m-dim array defined on \code{grids}, the "pre" scan.}
  \item{pre.list}{A list of m-dim array defined on \code{grids}, the
    "pre" scans.}
  \item{post}{An m-dim array defined on \code{grids}, the "post" scan.}
  \item{post.list}{A list of m-dim array defined on \code{grids}, the
    "post" scans. pre/post lists may have different number of scans.}
  \item{perms}{Number of spatial permutations to run.  See 'Details' for
  more explanations.}
  \item{rand.comb}{The number of random scan permutations. See 'Details'
  for an explanation.  The default value is 0, which means
  \code{rep.test()} will exhaust all combinations.}
  \item{balanced}{A logical variable used by the built-in multiple
    testing procedure which controls FWER by simultaneous confidence
    band.  Balanced: I assume the difference is symmetrically
    distributed around zero, so \code{abs(diff)} is used in the
    hypothesis testing.  Unbalanced: conduct two separate one-sided test
    for \code{diff} to obtain \code{p.up}, \code{p.down}; the two sided
    p-value is then defined as \code{min(1.0, 2*min(p.up, p.down))}.
    This p-value is conservative and asymptotically exact.  Default:
    balanced=TRUE.}
  \item{norm}{A list of functional norms to be used as summary
    statistics for pre/post test and as kernels in N-distance for
    repeated scans test.  Implemented norms: L1, L2, and Linf.  Since
    the computational cost is very small, by default these functions
    compute all three of them.}
  \item{method}{Smoothing method used by the generalized bootstrap
    p-value estimation.  Valid methods: gld (genral lambda distribution,
    depends on package \code{gld} and is quite slow), sn (skew-normal,
    depends on package \code{sn}), normal, and null (the usual
    nonparametric estimate). Default value: \code{c("null", "normal")}.}
  \item{MTP}{Specify the multiple testing procedure used to compute
    \code{fdr.adj.p}. Default to \code{BH}.  See help of
    \code{p.adjust()} for more options.}
  \item{...}{Other options passed to different smoothers such as a choice
    of kernel (\code{kernel=box}) or various options of function
    \code{Tps()} implemented in \pkg{fields}.}
}

\details{For repeated scan test (\code{rep.test()}), the total number of
  permutations is the number of spatial permutations (controlled by
  \code{perms}) times the number of pre/post group combinations of scan
  permutations due to the symmetry of the underlying hypothesis testing
  procedure.  Using this symmetry property can greatly reduce the number
  of kernel smoothings, which is the most expensive computation.  

  Let there be \eqn{n_1}{n1} pre scans and \eqn{n_2}{n2} post scans,
  there are \eqn{n_1 + n_2}{n1+n2} choose \eqn{n1}{n1} of them,
  which can be large but not astronomically large since nobody gets
  scanned for more than 10 times in a study in practice and
  \code{rep.test()} can happily compute the case of 8 versus 8 scans
  (set \code{perm=1} of course) within minutes.  Unless \code{rand.comb}
  is set manually, \code{rep.test()} exhausts all the combinations
  because I employ a computation trick to do this very efficiently (see
  our paper about this).  Consequently you don't need to specify a very
  large number of spatial permutations.  \cr

  However, for the sake of writing a robust function, and just in case
  someone will test this program in a completely unintended environment,
  parameter \code{rand.comb} can be specified so this function will do
  \code{rand.comb} many random scan permutations for each spatial
  permutation (controlled by \code{perms}).

  For usual nonparametric p-value estimation I recommend using 1000
  total permutations, so you need to set \code{perms=1000} for
  \code{pre.post.test()}, or \code{perms=50} for \code{rep.test()} with
  3 pre scans and 3 post scans (because 6 choose 3 is 20, 20 times 50
  makes 1000).  Because of this symmetry property, processing many
  repeated scans by\code{rep.test()} is usually much faster than doing a
  pre/post test! \cr

  To further reduce computation cost, you may want to consider using
  smoothed (generalized bootstrap) p-value estimation, which reduces the
  variance (so that you can reduce permutations) at the cost of a small
  biasedness (so either the type I error rate may go up or the detection
  power may go down).  For these tests I recommend 200 total
  permutations.
}

\value{A list with components
  \item{p.global}{A list of global p-values global p-values computed by
    resampling the specified functional norms of the difference map
    (\code{pre.post.test()}) or the N-distances with specified norms
    (kernels).}
  \item{rawp}{A list of arrays of unadjusted per-voxel p-values computed
  from permutations.}
  \item{scb.adj.p}{A list of arrays of per-voxel p-values adjusted for
    familywise error rate by a simultaneous confidence band with different
    generalized bootstrap method.}
  \item{wy.adj.p}{An array of per-voxel p-values adjusted for
    familywise error rate by the Westfall-Young stepdown procedure.}
  \item{fdr.adj.p}{An array of per-voxel p-values adjusted for
    false discovery rate by the Benjamini-Hochberg procedure.}
  \item{fit.diff}{An array of fitted (smoothed) difference map.}
  \item{band.up}{A vector of upper bounds of fitted difference map
    for each permutation.}
  \item{band.down}{A vector of lower bounds of fittted difference map
    for each permutation.}
}
\references{[my own paper], [westfall young], [papers about general bootstrap]}
\author{Xing Qiu}
\seealso{\code{\link{genboot.test}}, \code{\link{Smooth}},
    \code{\link{Ndist}}, \code{spatial.perm}}
\examples{
## Generate a small data
Nx <- 32; Ny <- 32; Nz <- 3
for (i in 1:5) assign(paste("FA",i,sep=""), array(rnorm(Nx*Ny*Nz), c(Nx,Ny,Nz)))
grid.x <- (1:Nx)*2-1; grid.y <- (1:Ny)*2-1; grid.z <- (1:Nz)*2-1; 
grids <- list(grid.x, grid.y, grid.z)
Xlist <- list(FA1, FA2, FA3); Ylist <- list(FA4, FA5)

## takes about a minute on my PC (4x 2.33GHz), YMMV.
\dontrun{r1 <- pre.post.test(grids, FA1, FA2)}

## takes about 20 seconds on my PC (4x 2.33GHz), YMMV.
\dontrun{r2 <- rep.test(grids, Xlist, Ylist, method="normal")}
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{htest}% __ONLY ONE__ keyword per line
