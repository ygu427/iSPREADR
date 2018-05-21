#### Collections of all smoothing related functions.

## This is a wrapper around the ksmooth() function.
.ksmooth.md <- function(grids, img, bandwidth=5.0, kernel="normal", ...){
  ## multi-dim ksmooth wrapper.  grids: a list of (marginal)
  ## coordinates. img: an m-dim array of values to be smoothed.
  ## dim(img) must be sapply(grids, length) otherwise this function
  ## won't work.
  ## Still needs some work for one dim array.
  Y <- as.array(img)
  if (any(dim(Y) != sapply(grids, length))){
    stop("Dimensions of coordinates (grids) and values (img) do not match.")
  } else K <- length(dim(Y))
  y.k <- Y
  for (k in 1:K){
    y.k <- aperm(apply(y.k, seq(K)[-k], function(y) ksmooth(grids[[k]], y, n.points=length(grids[[k]]), bandwidth=bandwidth, kernel=kernel, ...)$y), append(2:K, 1, k-1))
  }
  return(y.k)
}

## An uaxiliary function that estimates the gradient vector from a
## 27-point stencil by WLS. Used in AnisoDiff3.  I need to design a
## faster algorithm of this function (possibly combines with Div())
## and implement it in C.
Grad <- function(Img, voxel.spacing){
  dims <- dim(Img); nx <- dims[1]; ny <- dims[2]; nz <- dims[3]
  ## V0 is the set of 26 directions in voxel (not in mm)
  V0 <- as.matrix(expand.grid(c(0, 1, -1), c(0, 1, -1), c(0, 1, -1)))[-1,]
  ## V is measured in mm
  V <- V0 %*% diag(voxel.spacing)
  Delta <- array(0, c(dims, 26))
  for (k in 1:nrow(V0)) {
    ## Dk is one of the 26 directions (no length); dk is the length
    ## of the differential unit
    Dk <- V0[k,]; dk <- sqrt(sum((Dk*voxel.spacing)^2))
    ## Imgk is 2 units larger than Img; so that it can always
    ## "contain" the shifted image.
    Imgk <- array(0, dims+2)
    Imgk[(seq(1,nx)+1+Dk[1]), (seq(1,ny)+1+Dk[2]), (seq(1,nz)+1+Dk[3])] <- Img
    ## Now remove the extra space and take the difference; then
    ## divided by the length of differential to obtain Delta.k
    Delta[,,,k] <- (Imgk[seq(1,nx)+1,seq(1,ny)+1,seq(1,nz)+1] - Img)/dk
  }
  ## Now the WLS estimator
  W <- diag(26) - matrix(1, 26, 26)/27
  hatmat <- solve(t(V) %*% W %*% V) %*% t(V) %*% W
  grad.est <- apply(Delta, c(1,2,3), function(x) drop(hatmat %*% x))
  return(grad.est)
}

## An auxiliary function that computes numerical divergence from a
## 7-point stencil. Used in AnisoDiff3.  Here H is a (3, nx, ny, nz)
## array of weighted gradients.
Div <- function(H, voxel.spacing){
  dims <- dim(H)[-1]; nx <- dims[1]; ny <- dims[2]; nz <- dims[3]
  ## calculation is based on a 7-point stencil (3x2=6 directions)
  emat <- diag(3); divH <- array(0, c(nx,ny,nz))
  for (k in 1:nrow(emat)){
    ek <- emat[k,]; ek.neg <- -emat[k,]; dk <- voxel.spacing[k]
    ## Hk is 2 units larger than H, so that it always "contain" H
    Hk <- array(0, dims+2); Hk.neg <- array(0, dims+2)
    Hk[(seq(1,nx)+1+ek[1]), (seq(1,ny)+1+ek[2]), (seq(1,nz)+1+ek[3])] <- H[k,,,]
    Hk.neg[(seq(1,nx)+1-ek[1]), (seq(1,ny)+1-ek[2]), (seq(1,nz)+1-ek[3])] <- H[k,,,]
    ## remove the extra space and take the diff, then divided by voxel
    divH <- divH + (Hk[seq(1,nx)+1,seq(1,ny)+1,seq(1,nz)+1] - H[k,,,])/dk
  }
  return(divH)
}


## This is the anisotropic, nonlinear smoother based on Perona-Malik
## equation, modified by Xing. The crucial difference is that we
## estimate the gradient from 26 directions first; then calculate the
## diffusion coefficient g(|| Delta I ||), which is a more accurate
## reflection of the original PDE.


## wrapper function for different spatial smoothers
Smooth <- function(grids, img, sm.method="ksmooth", ...) {
  sm.method <- match.arg(sm.method)
  ys <- switch(sm.method,
               "ksmooth"=.ksmooth.md(grids, img, ...),
               stop("Valid smoothing methods ksmooth (kernel smoothing, the default choice), more methods are coming!))
  return(ys)
}
