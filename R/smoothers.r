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

## This is the thin plate spline smoother, depends on package
## "fields".  This function is too slow for large scale data, it is
## included as an example.
.tps.md <- function(grids, img, bandwidth=5.0, ...){
  x <- expand.grid(grids); Y <- as.vector(img)
  mod1 <- fastTps(x, Y, theta=bandwidth, ...)
  return(array(predict(mod1), dim(img)))
}

## This is the R implementation of the *original* algorithm based on
## the PM paper.  It uses a 7-point stencil to implement the 3D
## anisodiff.  Kappa is also the original kappa. See the Python
## reference implementation for details
AnisoDiff3.PM <- function(img, niter=3, kappa=.25, lambda=.25, voxel.spacing=c(1,1,1), option=2) {
  dims <- dim(img)
  if (length(dims) != 3) {
    stop("AnisoDiff only works for 3D grey-scale maps (such as FA and MD maps calculated from DTI analysis.")
  } else {
    nx <- dims[1]; ny <- dims[2]; nz <- dims[3]
  }
  ## initialize the diffused map by the original image
  imgOut <- img
  ## initialize some internal variables. Note that we only need one
  ## array for each direction (namely, one array for NS, EW, and UD,
  ## resp.
  NS <- EW <- UD <- deltaS <- deltaE <- deltaD <- array(0, dims)
  gS <- gE <- gD <- array(1, dims)
  ## the iterations
  for (k in 1:niter){
    ## calculate array differences
    deltaD[-nx, , ] <- imgOut[-1, , ] - imgOut[-nx, , ]
    deltaS[, -ny, ] <- imgOut[, -1, ] - imgOut[, -ny, ]
    deltaE[, , -nz] <- imgOut[, , -1] - imgOut[, , -nz]
    ## calculate g 
    if (option == 1) {
      gD <- exp(-(deltaD/kappa)^2)/voxel.spacing[1]
      gS <- exp(-(deltaS/kappa)^2)/voxel.spacing[2]
      gE <- exp(-(deltaE/kappa)^2)/voxel.spacing[3]
    } else if (option == 2){
      gD <- 1/(1 + (deltaD/kappa)^2)/voxel.spacing[1]
      gS <- 1/(1 + (deltaS/kappa)^2)/voxel.spacing[2]
      gE <- 1/(1 + (deltaE/kappa)^2)/voxel.spacing[3]
    } else {
      stop("Option can only take two values: 1 (the exponential P-M equation) and 2 (the reciprocal version).")
    }
    ## update matrices
    D <- gD*deltaD
    E <- gE*deltaE
    S <- gS*deltaS
    ## calculate the differences and shift Up/North/West by one voxel.
    UD <- D; UD[-1, , ] <- UD[-1, , ] - D[-nx, , ]
    NS <- S; NS[, -1, ] <- NS[, -1, ] - S[, -ny, ]
    EW <- E; EW[, , -1] <- EW[, , -1] - E[, , -nz]
    ## update the image
    imgOut <- imgOut + lambda*(UD+NS+EW)
  }
  return(imgOut)
}

## This is the anisotropic nonlinear smoother based on Perona-Malik
## equation, implemented by Bilan in MatLab and adapted by Xing to
## R. Note that I have changed kappa to the original scale for easy
## comparison. Bilan's original matlab code uses an estimation in
## Vocci et al. that doesn't work very well.
AnisoDiff3.BL <- function(img, niter=3, kappa=.25, lambda=.25, voxel.spacing=c(1,1,1), option=2) {
  dims <- dim(img)
  if (length(dims) != 3) {
    stop("AnisoDiff only works for 3D grey-scale maps (such as FA and MD maps calculated from DTI analysis.")
  } else {
    nx <- dims[1]; ny <- dims[2]; nz <- dims[3]
  }
  ## turn grids into voxel.spacing
  if (any(sapply(grids, function(xx) length(unique(diff(xx)))) != c(1, 1, 1))) {
    stop("Currently, AnisoDiff smoothing method only works for 3D equal spaced images.")
  }
  voxel.spacing <- sapply(grids, function(xx) unique(diff(xx)))
  ## initialize the diffused map by the original image
  imgD.next <- img
  ## generate all 26 directions of differentials (without lengths
  ## yet). Note that I need to remove c(0, 0, 0) from this list
  DD <- as.matrix(expand.grid(c(0, 1, -1), c(0, 1, -1), c(0, 1, -1)))[-1,]
  ## main iteration
  for (k in 1:niter) {
    ## imgD is the current diffused image; imgD.next is the updated imgD
    imgD <- imgD.next
    for (j in 1:nrow(DD)) {
      ## Dj is one of the 26 directions (no length); dj is the length
      ## of the differential unit
      Dj <- DD[j,]; dj <- sqrt(sum((Dj*voxel.spacing)^2))
      ## imgDj is 2 units larger than imgD; so that it can always
      ## "contain" the shifted image.
      imgDj <- array(0, dims+2)
      imgDj[(seq(1,nx)+1+Dj[1]), (seq(1,ny)+1+Dj[2]), (seq(1,nz)+1+Dj[3])] <- imgD
      ## Now remove the extra space and take the difference; then
      ## divided by the length of differential to obtain Delta.j
      Delta.j <- (imgDj[seq(1,nx)+1,seq(1,ny)+1,seq(1,nz)+1] - imgD)/dj
      if (option == 1) {
        cj <- exp(-(Delta.j/kappa)^2)
      } else if (option == 2){
        cj = 1/(1 + (Delta.j/kappa)^2)
      } else {
        stop("Option can only take two values: 1 (the exponential P-M equation) and 2 (the reciprocal version).")
      }
      ## We only implemented the adiabatic boundary condition here.
      cj[imgD == 0] <- 0
      ## update imgD.next
      imgD.next <- imgD.next + lambda*cj*Delta.j/(dj^2)
    }
  }
  return(imgD.next)
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
Smooth <- function(grids, img, sm.method=c("ksmooth", "anisodiff", "tps"), ...) {
  sm.method <- match.arg(sm.method)
  ys <- switch(sm.method,
               "ksmooth"=.ksmooth.md(grids, img, ...),
               "anisodiff"=.AnisoDiff3(grids, img, ...),
               "tps"=.tps.md(grids, img, ...),
               stop("Valid smoothing methods: 1. ksmooth (kernel smoothing, the default choice), 2. anisodiff (anisotropic smoother based on Perona-Malik equation). 2. tps (thin plate spline)."))
  return(ys)
}
