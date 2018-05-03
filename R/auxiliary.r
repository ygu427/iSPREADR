#### Collections of useful auxiliary functions

## My crude attempt to remove ringing effects
MountDoom <- function(img, thickness=2, minval=0){
  dims <- dim(img); nx <- dims[1]; ny <- dims[2]; nz <- dims[3]
  ## Identify a crude background (==0) mask
  mask0 <- ifelse(img<=minval, 0, 1)
  ## enlarge the background mask by 26 neighbors. First, let us
  ## generate all 26 directions of increments (without lengths
  ## yet). Note that I need to remove c(0, 0, 0) from this list
  DD <- as.matrix(expand.grid(c(0, 1, -1), c(0, 1, -1), c(0, 1, -1)))[-1,]
  ## with length constraints
  DD2 <- data.frame()
  for (i in 1:nrow(DD)){
    Li <- sqrt(sum(DD[i,]^2))           # length
    if ( Li <= thickness){
      steps <- seq(1, floor(thickness / Li))
      DD2 <- rbind(DD2, steps %x% t(DD[i,]))
    }
  }
  ## Now use these directions to enlarge the background
  mask <- mask0
  if (nrow(DD2)>0){
    for (j in 1:nrow(DD2)){
      Dj <- as.numeric(DD2[j,])
      ## maskj is 2*thickness units larger than mask; so that it can always
      ## "contain" the shifted image.
      maskj <- array(0, dims+2*thickness)
      maskj[(seq(1,nx)+thickness+Dj[1]), (seq(1,ny)+thickness+Dj[2]),
      (seq(1,nz)+thickness+Dj[3])] <- mask0
      ## Now remove the extra space and take the difference; then
      ## divided by the length of differential to obtain Delta.j
      maskj2 <- maskj[seq(1,nx)+thickness,seq(1,ny)+thickness,seq(1,nz)+thickness]
      ## enlarge the original mask
      mask <- mask * maskj2
    }
  }
  return(mask)
}


## grids: a LIST of monotonically increasing coordinate grids.  For
## example, grids <- list(grid.x=seq(0, 100, 2), grid.y=exp(seq(1, 5,
## .2))).

## This function uses a pre-scribed threshold to detect
## foreground/background from a list of images (XYlist). It then tries
## to find the smallest rectangular prism to contain the foreground
foreground <- function(grids, XYlist, threshold){
  ## fg.min/max defines the cuboid that contains all foregrounds
  fg.min <- sapply(grids, max)          #starts with the worst case
  fg.max <- sapply(grids, min)          #starts with the worst case
  for (img in XYlist){
    fg.min.i <- apply(which(img >= threshold, arr.ind=TRUE), 2, min)
    fg.max.i <- apply(which(img >= threshold, arr.ind=TRUE), 2, max)
    fg.min <- pmin(fg.min, fg.min.i)
    fg.max <- pmax(fg.max, fg.max.i)
  }
  ## produce a mask first
  Ns <- sapply(grids, length); mask <- array(0, Ns)
  ## use a trick to write a general algorithm
  cmd <- paste0("mask[", paste(paste(fg.min, fg.max, sep=":"), collapse=","), "] <- 1")
  eval(parse(text = cmd))
  ## Now use fg information to construct smaller images and coordinate
  ## grids.
  grids.new <- lapply(1:length(grids), function(i) grids[[i]][fg.min[i]:fg.max[i]])
  XYlist.new <- lapply(1:length(XYlist), function(i) extract(XYlist[[i]], indices=lapply(1:length(fg.min), function(i) seq(fg.min[i], fg.max[i]))))
  return(list(grids=grids.new, XYlist=XYlist.new, mask=mask, fg.min=fg.min, fg.max=fg.max))
}

## A useful function that simulates some "true signal" in simulation
## analyses
artifacts <- function(grids, center = rep(0, length(grids)), std = .5, radius = 10, maxval = 1, minval=.2, xscale = diag(length(grids)))
{
  Ns <- sapply(grids, length)
  lambda1 <- max(abs(eigen(xscale, only.values = TRUE)$values))
  L1 <- lambda1 * radius
  grids2 <- NULL
  inds <- NULL
  for (n in 1:length(Ns)) {
    grids2[[n]] <- grids[[n]] - center[n]
    inds[[n]] <- seq(Ns[n])[grids2[[n]] >= -L1 & grids2[[n]] <=
                              L1]
    grids2[[n]] <- grids2[[n]][inds[[n]]]
  }
  coords <- solve(xscale) %*% t(expand.grid(grids2))
  func1 <- function(xvec) {
    if (sum(xvec^2) > radius^2) {
      s <- 0
    } else {
      s <- prod(dnorm(xvec, sd = std*radius))
    }
    return(s)
  }
  art1 <- array(apply(coords, 2, func1), sapply(inds, length))
  W <- array(0, Ns)
  cmd <- paste("W", "[", paste(inds, collapse = ","), "]",
               "<-", "art1", sep = "")
  eval(parse(text = cmd))
  ## normalize W so that the maximum == maxval
  W <- maxval*W/max(W)
  ## get rid of very small values
  W <- ifelse(W<minval, 0, W)
  return(W)
}

## compute the L1, L2, L-inf norm (over number of observations) of a
## fitted object.  Assuming EQUAL grid.
norms <- function(Y, norm=c("L1", "L2", "Linf")){
  rr <- NULL
  for (nn in norm){
    rr[[nn]] <- switch(nn, "L1"=mean(abs(Y)),
                       "L2"=sqrt(mean(Y^2)),
                       "Linf"=max(abs(Y)),
                       stop("Norms implemented: L1, L2, and Linf."))
  }
  return(rr)
}

## A C implementation of multivariate N-stat. This procedure is 30
## times faster than ksmooth, so a) I don't have to re-write it in C;
## b) I can afford to compute all 3 kernels (the default option of
## function norms()).  The C implementation
Ndist <- function(Xlist, Ylist, kernels=c("L1", "L2", "Linf")){
  Nx <- length(Xlist); Ny <- length(Ylist); N <- Nx + Ny
  XYlist <- c(Xlist,Ylist)
  dist.pairs <- array(0, c(length(kernels),N,N), dimnames=list(kernels,1:N,1:N))
  ##lower triangle only.
  for (i in 2:N){
    for (j in 1:(i-1)){
      d.ij <- norms(XYlist[[i]] - XYlist[[j]])
      dist.pairs[,j,i] <- d.ij;
    }}

  ndists <- 1:length(kernels); names(ndists) <- kernels
  for (kn in 1:length(kernels)){
    ndist.k <- 0.0
    ndists[kn] <- .C("ndist_from_distmat", as.double(dist.pairs[kn,,]),
                     as.integer(N), as.integer(Nx),
                     ndist=as.double(ndist.k), PACKAGE = "iSPREADR")$ndist
  }
  return(ndists)
}

## This is the permutation version of Ndist.  It employes a trick to
## GREATLY reduce the computing time: it first computes all pairwise
## distances, then use permutation of these NUMBERS to compute the
## permuted N-stats.  For flexibility combs is an external list of
## permutations/combinations which can be generated by either
## combn(N,Nx,function(x) c(x, setdiff(1:N, x)), simplify=FALSE) or
## foreach(icount(rand.comb)) %do% sample(N).  See rep.test()
## for more details.

## The C implementation can be 40~50 times faster than the R
## implementation for large number of combinations.
Ndist.perm <- function(Xlist, Ylist, combs, kernels=c("L1", "L2", "Linf")){
  Nx <- length(Xlist); Ny <- length(Ylist); N <- Nx + Ny;
  K <- length(combs); combsvec <- do.call("c", combs)-1
  XYlist <- c(Xlist,Ylist)
  dist.pairs <- array(0, c(length(kernels),N,N), dimnames=list(kernels,1:N,1:N))
  for (i in 2:N){
    for (j in 1:(i-1)){
      d.ij <- norms(XYlist[[i]] - XYlist[[j]])
      dist.pairs[,j,i] <- d.ij
    }}
  my.Ns <- matrix(0, nrow=length(combs), ncol=length(kernels))
  colnames(my.Ns) <- kernels
  for (kn in 1:length(kernels)){
    ndistvec.k <- rep(0,K)
    my.Ns[,kn] <- .C("ndist_from_distmat_perm", as.double(dist.pairs[kn,,]),
                     as.integer(combsvec), as.integer(K), as.integer(N),
                     as.integer(Nx), ndistvec=as.double(ndistvec.k),
                     PACKAGE = "iSPREADR")$ndistvec
  }
  return(my.Ns)
}

## flip: sort of wild bootstrap for an array.  Group action: G \cong 2^N.
flip <- function(diff.array){
  Y <- as.array(diff.array)
  g <- 2*array(rbinom(length(Y), 1, .5), dim(Y))-1
  return(g*Y)
}

## XYlist: a list of n arrays which are i.i.d. random vectors of length
## N under H0.  Permutation group action: G \cong n^N.
spatial.perm <- function(XYlist){
  n <- length(XYlist); Ns <- dim(XYlist[[1]])
  Amat <- matrix(unlist(XYlist), nrow=prod(Ns))
  Bmat <- t(apply(Amat, 1, sample))
  foreach(col=iter(Bmat)) %do% array(col, Ns)
}

## genboot.test().  Can be used in p-value evaluation when the number
## of permutations is small (<500). The generalized bootstrap is a
## smoothed bootstrap method (Chernick 1992 book, 6.2.2.).

## The default smoothing method is based on the generalized (Tukey)
## lambda distribution (Dudewicz, 1992).  Parameter estimation and
## probability computing is Implemented by R package "gld".

## Alternative methods include Gaussian (which is called "parametric
## bootstrap" by Efron and skew-normal distribution (3 parameter
## generalization of normal), implemented by R package "sn"; and may
## in the future include others.

## Efforts are made so to let this function deal gracefully with
## identical null vector case.  nullvec: typically a vector of summary
## statistics (like norms) generated by resampling under H0; x:
## typically the summary statistic computed from the original (could
## be H1) data.  Statistical validity: this summary statistic follows
## an asymptotic normal distribution under H0.

## The default alternative is "greater" because usually norm>0 and we
## want to show that x > nullvec.
.genboot.test.identical <- function(nullvec, x, alternative="greater"){
  ERR <- 10^-9
  x0a <- min(nullvec)-ERR; x0b <- max(nullvec)+ERR
  if (x > x0b){
    intv <- "high"
  } else if (x < x0a){
    intv <- "low"
  } else {
    intv <- "med"
  }
  pp <- switch(alternative,
               "greater"=switch(intv, high=0.0, med=0.5, low=1.0),
               "less"=switch(intv, high=1.0, med=0.5, low=0.0),
               "two.sided"=switch(intv, high=0.0, med=0.5, low=0.0),
               stop("Valid alternatives: greater, less, two.sided."))
  return(pp)
}

.genboot.test.gld <- function(nullvec, x, alternative="greater", ...){
  ## I choose this initgrid because it is way faster than the default.
  params <- starship(nullvec, initgrid=list(lcvect=(-2:2)/10, ldvect=(-2:2)/10), ...)$lambda
  prob <- do.call("pgl", c(list(x), params))
  pp <- switch(alternative,
               "greater"= 1-prob,
               "less"=prob,
               "two.sided"=2*(1-prob),
               stop("Valid alternatives: greater, less, two.sided."))
  return(pp)
}

.genboot.test.normal <- function(nullvec, x, alternative="greater"){
  params <- list("mean"=mean(nullvec), "sd"=sd(nullvec))
  prob <- do.call("pnorm", c(list(x), params))
  pp <- switch(alternative,
               "greater"= 1-prob,
               "less"=prob,
               "two.sided"=2*(1-prob),
               stop("Valid alternatives: greater, less, two.sided."))
  return(pp)
}

.genboot.test.sn <- function(nullvec, x, alternative="greater", ...){
  ## depends on package "sn"
  params <- sn.em(y=nullvec, ...)$dp
  prob <- do.call("psn", c(list(x), params))
  pp <- switch(alternative,
               "greater"= 1-prob,
               "less"=prob,
               "two.sided"=2*(1-prob),
               stop("Valid alternatives: greater, less, two.sided."))
  return(pp)
}

.genboot.test.null <- function(nullvec, x, alternative="greater"){
  ## A convenient function wrapper for nonparametric bootstrap test.
  prob <- sapply(x, function(xi) sum(xi>nullvec))/length(nullvec)
  pp <- switch(alternative,
               "greater"= 1-prob,
               "less"=prob,
               "two.sided"=2*(1-prob),
               stop("Valid alternatives: greater, less, two.sided."))
  return(pp)
}

genboot.test <- function(nullvec, x, alternative="greater", method="normal", ...){
  ERR <- 10^-9
  if (max(nullvec) - min(nullvec) < ERR) {
    pp <- .genboot.test.identical(nullvec, x, alternative)
  } else {                              #None identical vector case
    pp <- switch(method,
                 "null"=.genboot.test.null(nullvec, x, alternative),
                 "normal"=.genboot.test.normal(nullvec, x, alternative),
                 "sn"=.genboot.test.sn(nullvec, x, alternative, ...),
                 "gld"=.genboot.test.gld(nullvec, x, alternative, ...),
                 stop("Valid methods: gld (genral lambda distribution), sn (skew-normal), normal, and null (the usual nonparametric estimate)."))
  }
  return(pp)
}

## Computes the reverse rank of x in x.perm.  This function can be
## used to compute permutation p-values quickly.  x can be a vector.

rev.rank <- function(x, x.perm) {
  y <- sort(x.perm, decreasing=TRUE)
  o <- order(x, decreasing=TRUE); ro <- order(o)
  rv <- rep(-1,length(x))
  z <- .C("vec_rev_rank", as.double(x[o]), as.double(y),
          as.integer(length(x)), as.integer(length(y)),
          rv=as.integer(rv), PACKAGE = "iSPREADR")
  return(z$rv[ro])
}

#### converts table format <==> array format

## XYlist: a list of (at least 2d) arrays. zero.rm=TRUE: remove points
## with constant zero values across all arrays. EPSILON: if the
## absolute value of an entry is less than EPSILON it will be reset to
## zero.  You can use this to do simple thresholding.
arrays2tab <- function(grids, XYlist, zero.rm=FALSE, EPSILON=0.1^6){
  Ns <- dim(XYlist[[1]]); NN <- prod(Ns); Nd <- length(Ns)
  Nlist <- length(XYlist)
  coords <- expand.grid(grids)
  colnames(coords) <- paste("X", 1:Nd, sep="")
  vals <- matrix(0, NN, Nlist)
  colnames(vals) <- paste("Scan", 1:Nlist, sep="")
  for (j in 1:Nlist) vals[,j] <- Xlist[[j]]
  tab <- cbind(coords, vals)
  if (zero.rm){
    zero.ind <- apply(abs(vals)<EPSILON, 1, all)
    return (tab[!zero.ind, ])
  } else {
    return (tab)
  }
}

## This is roughly the inverse of arrays2tab. It assumes the first Nd
## columns are the coordinates.  grids still need to be specified
## because a) otherwise there might be places that should be in the
## grid but not detectable automatically; b) it provides a safe guard
## against generating a 2M by 2M by 2M array ---- it will eat all your
## memory and bring down your box in just seconds.  default.val: value
## to be filled at blank voxels.
tab2arrays <- function(grids, tab, default.val=0){
  Ns <- sapply(grids, length); NN <- prod(Ns); Nd <- length(Ns)
  Nlist <- dim(tab)[2] - Nd
  ## sanity check.  Ensure that the unique values of the coordinates
  ## supplied by tab is a subset of those supplied by grids.
  for (k in 1:Nd){
    coords.tab <- unique(tab[,k])
    sanity <- coords.tab %in% grids[[k]]
    if (!all(sanity)) stop(paste("Coordinates mismatch. The", k, "column of table has the following element(s) which are not in the grids:", paste(coords.tab[!sanity], collapse=", ")))
  }
  ## compute the correspondences just once.  I need to compute the
  ## indices in vector (1d array) form because there is no easy and
  ## fast way of slicing an array by a given index vector.
  W <- c(1,cumprod(Ns))                 #weights for each array
  inds <- 1 + foreach (k=1:Nd, .combine="+") %dopar% {
    (Ns[k] - rev.rank(tab[,k], grids[[k]]))*W[k]
  }
  ## Now assign values from tab to arrays
  Alist <- foreach (j=1:Nlist) %dopar% {
    X.k <- array(default.val, Ns)
    for (i in 1:dim(tab)[1]) X.k[inds[i]] <- tab[i,j+Nd]
    X.k
  }
  return(Alist)
}
