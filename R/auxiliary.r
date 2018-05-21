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


## flip: sort of wild bootstrap for an array.  Group action: G \cong 2^N.
flip <- function(diff.array){
  Y <- as.array(diff.array)
  g <- 2*array(rbinom(length(Y), 1, .5), dim(Y))-1
  return(g*Y)
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
