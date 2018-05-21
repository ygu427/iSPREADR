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


## To be implemented: paired test and group test.

## This version doesn't detect a common diff. region.  Maybe in the
## future I can think of some reasonable way of detecting common
## diff. region.

