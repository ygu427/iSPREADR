\name{genboot.test}
\alias{genboot.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized bootstrap (permutation) test.}
\description{
  Generalized bootstrap (permutation) test.
}
\usage{
genboot.test(nullvec, x, alternative="greater", method="normal", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nullvec}{A vector of statistics sampled under the null
    hypothesis by either bootstrap or permutation.}
  \item{x}{The original statistic(s) to be compared with the resampling
    generated null distribution.}
  \item{alternative}{Valid alternatives: greater, less, two.sided.
    Default to "greater", i.e., x is "significantly greater than the
    null distribution." I decide to use "greater" as the default option
    because in my study, most summary statistics are functional norms
    which warrant one sided test.}
  \item{method}{Smoothing method used by the generalized bootstrap
    p-value estimation.  Valid methods: gld (genral lambda
    distribution), sn (skew-normal), normal, and null (the usual
    nonparametric estimate).}
  \item{...}{other optional parameters passed to \code{starship()} or
    \code{sn.em()}, which are the parameter estimating function for
    general lambda distribution and skew-normal distribution, resp.}
}
\details{
  This function takes one (or a vector of) original statistic \code{x} and a
  vector of the same statistic calculated under the null hypothesis by
  resampling method (\code{nullvec}) to produce a p-value.  When
  \code{method="null"}, it uses the common nonparametric estimate, i.e.,
  just compute the proportion of members of \code{nullvec} greater or
  equal to \code{x}.  When method is "gld", "sn", or "normal", the null
  vector is used to fit a general lambda distribution (with four
  parameters), a skew-normal distribution (three parameters), or a
  normal distribution (two parameters); p-value is then computed from
  this estimated distribution.
}
\value{An estimated p-value or a vector of such p-values.}
\references{[Chernick 1992 book, 6.2.2.], [Dudewicz, 1992], and a
  reference to the sn package.}
\author{Xing Qiu}
\seealso{\code{\link{starship}}, \code{\link{sn.em}}}
\examples{
nv <- sort(round(rnorm(200),2))
## These two are very fast
genboot.test(nv, 1:3, method="null")
genboot.test(nv, 1:3, method="normal")

## The following commands require loading the sn and gld packages.
library(sn)
library(gld)
genboot.test(nv, 1:3, method="sn")

## This is quite slow. It take about 10 seconds to run on my PC (4x 2.33GHz), YMMV.
\dontrun{genboot.test(nv, 1:3, method="gld")}
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{smooth}
\keyword{htest}% __ONLY ONE__ keyword per line
