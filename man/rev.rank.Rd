\name{rev.rank}
\alias{rev.rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the reverse rank of one vector in another vector}

\description{ Computes the reverse rank of \code{x} in \code{x.perm}.
  This function can be used to compute permutation p-values
  quickly. Suppose \code{x} is a statistics of interest and
  \code{x.perm} is a vector of this statistic computed by a resampling
  method which simulates the null hypothesis, \code{rx <- rev.rank(x,
  x.perm)} computes number of \code{x.perm} which is greater or equal
  to \code{x}, so \code{rx/length(x.perm)} is the usual (nonparametric)
  permutation p-value for testing a one sided alternative hypothesis.
}

\usage{
rev.rank(x, x.perm)
}
\arguments{
  \item{x}{A number or a vector of numbers. It usually is a statistics
    of interest in a permutation test.}
  \item{x.perm}{A vector of numbers with which \code{x} is
    compared. Usually this is the same statistic computed by a resampling
    method which simulates the null hypothesis.}
  }

\details{
  This function is designed to compute permutation p-values quickly.
  \code{x} can be either a number or a vector, it is the second case for
  which this function (algorithm) is written: it is designed to take
  advantage of what has been ordered and ranked for \code{xo[1:k]} when
  it computes the rank of \code{xo[k+1]}, where \code{xo} is the ordered
  \code{o}. As a result it is much faster (more than 100 folds) than a
  naive serial algorithm such as \code{sapply(x, function(x.i) sum(x.i
  <= x.perm))}, and significantly (about 30\%) faster than the following
  pure R implementation: \code{length(x) - rank(c(x, perm.x),
  ties.method="min")[1:length(x)] + rank(xvec, ties.method="min")}
}
\value{An integer or a vector of the reverse rank of \code{x} in \code{x.perm}.}

\author{Xing Qiu}
\seealso{\code{\link{rank}}}
\examples{
## The two approaches are equivalent
xvec <- c(0.5, 2.5); x.perm <- 0:2
sapply(xvec, function(x) sum(x <=x.perm))
rev.rank(xvec, x.perm)

## Now let us generate two large vectors and compare the performance of
## different approaches

xvec1 <- rnorm(2*10^4); x.perm1 <- rnorm(3*10^4)
## On my computer, the following naive approach took about 10sec
\dontrun{system.time(p.i <- sapply(xvec1, function(x) sum(x <=x.perm1)))}
## it took only 0.017sec
system.time(p.i2 <- rev.rank(xvec1, x.perm1))

## Yet another approach
func2 <- function(xvec, x.perm){
  nx <- length(xvec)
  r1 <- rank(xvec, ties.method="min")
  r2 <- rank(c(xvec, x.perm), ties.method="min")[1:nx]
  return(nx - r2 + r1)
}
## about 0.022 sec
system.time(p.i3 <- func2(xvec1, x.perm1))
}                               % end examples.

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{utilities}% __ONLY ONE__ keyword per line
