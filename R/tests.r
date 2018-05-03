## my own implementation of FWER controlling MTP based on SCB (equal
## variance version).  Future work: est. STD for each voxel and use it
## as a weight.  It has to be done in the pre.post.test() function
## because with weight, the rule of finding max/min will be different.
## [update 05/03/11] No, don't implement STD at this level at all.  Do
## variable transformations based on directional variance first; then
## do pre.post.test as usual and be done with the individual adjusted
## p.values. Optional: transform the bands back to the original scale
## if we have to report conf. band.
.p.scb <- function(y, band.up, band.down, balanced=TRUE, method=c("null", "normal"), ...){
  ## (FWER) adjusted p-value based on the simultaneous confidence
  ## band. Balanced: use one sided version for abs(diff); unbalanced:
  ## times 2 is conservative/asymptotically exact.
  Ns <- dim(y)
  if (balanced){
    bands <- pmax(abs(band.up), abs(band.down))
    ## bands <- abs(c(band.up, band.down))  this doesn't work ... needs more study!
    pscb <- foreach (m=method) %dopar% array(genboot.test(bands, abs(y), method=m, ...), Ns)
  } else {
    pscb <- foreach (m=method) %dopar% {
      p.up <- genboot.test(band.up, y, method=m, ...)
      p.down <- genboot.test(band.down, y, alternative="less", method=m, ...)
      array(pmin(1.0, 2*pmin(p.up, p.down)), Ns)
    }
  }
  names(pscb) <- method
  return(pscb)
}

## A customized combine function
.combfun <- function(x1, x2) {
  list("norms"=rbind(x1$norms, x2$norms),
       "band.up"=c(x1$band.up, x2$band.up),
       "band.down"=c(x1$band.down, x2$band.down),
       "wycounts"=x1$wycounts + x2$wycounts, #needs to be divided later
       "pcounts"=x1$pcounts + x2$pcounts) #needs to be divided later
}

## Wrapper for doing the global and pointwise inference for 1 pre and
## 1 post images.

## The stopping rule is currently not implemented due to a) a
## deficiency of foreach loop; b) I need to compute the distribution
## of p-value under H0 anyway.  If implemented, stop.p will be a
## (conservative) stopping rule based on cumulative nonparametric
## global p-values.  1.0 means no stop at all, which is necessary for
## est. p-value /distribution/ under H0.  Otherwise, use 0.05 can save
## computing time dramatically.
pre.post.test <- function(grids, pre, post, perms=1000, balanced=TRUE, norm=c("L1","L2","Linf"), method=c("null", "normal"), ...){
  DIFF <- post - pre                       #the difference map
  Ns <- dim(DIFF); NN <- prod(Ns)
  fit.orig <- Smooth(grids, DIFF, ...)
  norms.orig <- norms(fit.orig)
  ## for WY proc 
  S <- abs(fit.orig)                  #summary stats for pre/post test
  ord <- order(S); ro <- order(ord); Sord <- S[ord]
  ## The main permutation loop
  rr <- foreach (j=1:perms, .combine=.combfun) %dopar% {
    DIFF.j <- flip(DIFF)
    fit.diff <- Smooth(grids, DIFF.j, ...)
    Sk <- abs(fit.diff)
    Uk <- cummax(Sk[ord])
    ## return these values to the master node
    list("norms"=norms(fit.diff, norm=norm),
         "band.up"=max(fit.diff),
         "band.down"=min(fit.diff),
         "pcounts"=rev.rank(S, Sk),
         "wycounts"=as.integer(Uk >= Sord)
         )
  }
  rawp <- rr$pcounts/(perms*NN)
  wy.adj.p.ord <- array(cummax(rr$wycounts[NN:1])[NN:1],Ns)/perms
  wy.adj.p <- array(wy.adj.p.ord[ro], Ns)
  band.up <- rr$band.up; band.down <- rr$band.down
  ## Inference
  p.global <- foreach (n=norm, .combine="cbind") %:% foreach (m=method, .combine="c") %dopar% {
    genboot.test(rr$norms[,n], norms.orig[n], method=m)
  }; dimnames(p.global) <- list(method, norm)
  scb.adj.p <- .p.scb(fit.orig, band.up, band.down, balanced=balanced, method=method)
  return(list("p.global"=p.global,
              "norms.orig"= norms.orig,
              "norms.perm"= rr$norms,
              "rawp"=array(rawp,Ns),
              "scb.adj.p"=scb.adj.p,
              "wy.adj.p"=wy.adj.p,
              "fdr.adj.p"=array(p.adjust(rawp, "BH"), Ns),
              "fit.orig"=fit.orig,
              "band.up"=band.up,
              "band.down"=band.down))
}


.list.mean <- function(Alist){
  ## when I have time, I need to rewrite the below function in C
  ## because this is of order N choose Nx, and it will become the
  ## bottleneck when N choose Nx is > 175*Nx, or about 7 scans in
  ## each group.
  ss <- foreach(aa=Alist, .combine="+") %do% aa
  return(ss/length(Alist))
}

## rep.test.  pre.list, post.list are repeated pre/post scans from the
## same subject.  This is to ensure the exchangeability of spatial
## points.  The region is detected by mean diff, which is the best I
## can think of at this moment because N-dist itself is a *summary*
## statistic, it is not designed to make per-voxel inference.

## 05/23/2011.  Added the symmetry trick. See pre-print for more
## details.  Option rand.comp is the number of random combinations to
## be sampled (second layer resampling).  It defaults to not using
## random combinations because it is simply crazy to assume somebody
## have more than 20 pre/post scans!  However, just for the sake of
## writing a robust function, and just in case some morons will test
## this program in a completely unintended environment, I decide to
## add this layer of flexibility anyway.
rep.test <- function(grids, pre.list, post.list, perms=50, rand.comb=0, balanced=TRUE, norm=c("L1","L2","Linf"), method=c("null", "normal"), MTP="BH", ...){
  Nx <- length(pre.list); Ny <- length(post.list); N <- Nx+Ny
  Ns <- dim(pre.list[[1]]); NN <- prod(Ns)
  if (rand.comb == 0){ #enumerate ALL combinations. Works for N <= 20.
    combs <- combn(N, Nx, function(x) c(x,setdiff(1:N,x)), simplify=FALSE)
  } else { #random permutations are the SAME as random combinations, see pre-print for more details.
    combs <- foreach(icount(rand.comb)) %do% sample(N)
  }
  K <- length(combs)
  X.ks <- foreach (x=pre.list) %dopar% {
    Smooth(grids, x, ...)
  }
  Y.ks <- foreach (y=post.list) %dopar% {
    Smooth(grids, y, ...)
  }
  fit.orig <- .list.mean(Y.ks) - .list.mean(X.ks)
  norms.orig <- Ndist(X.ks, Y.ks, kernel=norm)
  ## WY proc
  S <- abs(fit.orig)
  ord <- order(S); ro <- order(ord); Sord <- S[ord]

  rr <- foreach (j=1:perms, .combine=.combfun) %dopar% {
    XY <- spatial.perm(c(pre.list, post.list))
    XY.ks <- foreach (x=XY) %do% {
      Smooth(grids, x, ...)
    }
    ## Now the loop for all combinations
    band.up <- rep(0,K); band.down <- rep(0,K); pcounts <- rep(0,NN)
    wycounts <- rep(0,NN)
    for (k in 1:K){
      X.ind <- combs[[k]][1:Nx]; Y.ind <- combs[[k]][(Nx+1):N]
      fit.diff <- .list.mean(XY.ks[X.ind]) - .list.mean(XY.ks[Y.ind])
      Sk <- abs(fit.diff)
      Uk <- cummax(Sk[ord])
      band.up[k] <- max(fit.diff); band.down[k] <- min(fit.diff)
      pcounts <- pcounts + rev.rank(S, Sk)
      wycounts <- wycounts + as.integer(Uk >= Sord)
    }
    ## return these values to the master node
    list("norms"=Ndist.perm(XY.ks[1:Nx], XY.ks[(Nx+1):N], combs, kernel=norm),
         "band.up"=band.up,
         "band.down"=band.down,
         "pcounts"=pcounts,
         "wycounts"=wycounts
         )
  }
  rawp <- rr$pcounts/(perms*K*NN)
  wy.adj.p.ord <- array(cummax(rr$wycounts[NN:1])[NN:1],Ns)/(perms*K)
  wy.adj.p <- wy.adj.p.ord[ro]
  band.up <- rr$band.up; band.down <- rr$band.down
  ## Inference
  p.global <- foreach (n=norm, .combine="cbind") %:% foreach (m=method, .combine="c") %dopar% {
    genboot.test(rr$norms[,n], norms.orig[n], method=m)
  }; dimnames(p.global) <- list(method, norm)
  scb.adj.p <- .p.scb(fit.orig, band.up, band.down, balanced=balanced, method=method)
  return(list("p.global"=p.global,
              "scb.adj.p"=scb.adj.p,
              "wy.adj.p"=wy.adj.p,
              "fdr.adj.p"=array(p.adjust(rawp, MTP), Ns),
              "fit.orig"=fit.orig,
              "band.up"=band.up,
              "band.down"=band.down))
}

## To be implemented: paired test and group test.

## This version doesn't detect a common diff. region.  Maybe in the
## future I can think of some reasonable way of detecting common
## diff. region.

