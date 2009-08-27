##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## $Id: flipangle.R 112 2009-08-12 13:34:53Z bjw34032 $
##

#############################################################################
## dam() = double-angle method
#############################################################################

dam <- function(low, high, low.deg) {
  alpha <- acos(abs(high /(2*low)))
  (180/pi * alpha) / low.deg # radians to degrees
}

#############################################################################
## R10.lm() = estimate R1 using Levenburg-Marquardt
#############################################################################

R10.lm <- function(signal, alpha, TR, guess, nprint=0) {
  func <- function(x, y) {
    R1 <- x[1]
    m0 <- x[2]
    signal <- y[[1]]
    theta <- pi/180 * y[[2]]    # degrees to radians
    TR <- y[[3]]
    signal -
      m0 * sin(theta) * (1 - exp(-TR*R1)) / (1 - cos(theta) * exp(-TR*R1))
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=list(nprint=nprint, maxiter=150),
               y=list(signal, alpha, TR))
  list(R1=out$par[1], S0=out$par[2], info=out$info, message=out$message)
}

#############################################################################
## E10.lm() = estimate exp(-TR*R1) using Levenburg-Marquardt
#############################################################################

E10.lm <- function(signal, alpha, guess, nprint=0) {
  func <- function(x, signal, alpha) {
    E1 <- x[1]
    m0 <- x[2]
    theta <- pi/180 * alpha    # degrees to radians
    signal -
      m0 * sin(theta) * (1 - E1) / (1 - cos(theta) * E1)
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=list(nprint=nprint, maxiter=150),
               signal=signal, alpha=alpha)
  list(E10=out$par[1], m0=out$par[2], info=out$info, message=out$message)
}

#############################################################################
## setGeneric("R1.fast")
#############################################################################

setGeneric("R1.fast", function(flip, ...) standardGeneric("R1.fast"))
setMethod("R1.fast", signature(flip="array"),
          function(flip, flip.mask, fangles, TR, verbose=FALSE) 
	    dcemriWrapper("R1.fast", flip, flip.mask, fangles, TR, verbose))

#############################################################################
## R1.fast()
#############################################################################

.R1.fast <- function(flip, flip.mask, fangles, TR, verbose=FALSE) {

  if (length(dim(flip)) != 4)  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  if (!is.logical(flip.mask))  # Check flip.mask is logical
    stop("Mask must be logical.")
    
  X <- nrow(flip) ; Y <- ncol(flip) ; Z <- nsli(flip)
  nvoxels <- sum(flip.mask)
  
  if (verbose)
    cat("  Deconstructing data...", fill=TRUE)
  flip.mat <- matrix(flip[flip.mask], nrow=nvoxels)
  R10 <- M0 <- numeric(nvoxels)
  if (is.array(fangles))
    fangles.mat <- matrix(fangles[flip.mask], nrow=nvoxels)
  else
    fangles.mat <- matrix(fangles, nrow=nvoxels, ncol=length(fangles),
                          byrow=TRUE)

  if (verbose)
    cat("  Calculating R10 and M0...", fill=TRUE)
  for (k in 1:nvoxels) {
    fit <- E10.lm(flip.mat[k,], fangles.mat[k,],
                  guess=c(1, mean(flip.mat[k,])))
    if (fit$info == 1 || fit$info == 2 || fit$info == 3) {
      R10[k] <- log(fit$E10) / -TR
      M0[k] <- fit$m0
    } else {
      R10[k] <- M0[k] <- NA
    }
  }

  if (verbose)
    cat("  Reconstructing results...", fill=TRUE)
  R10.array <- M0.array <- array(NA, c(X,Y,Z,1))
  R10.array[flip.mask] <- R10
  M0.array[flip.mask] <- M0

  list(M0 = M0.array, R10 = R10.array)
}

#############################################################################
## setGeneric("CA.fast")
#############################################################################

setGeneric("CA.fast", function(dynamic, ...) standardGeneric("CA.fast"))
setMethod("CA.fast", signature(dynamic="array"),
	  function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
	      verbose=FALSE) 
	    dcemriWrapper("CA.fast", dynamic, dyn.mask, dangle, flip, fangles,
		TR, r1, verbose))

#############################################################################
## CA.fast() = estimate contrast-agent concentration and other stuff
#############################################################################

.CA.fast <- function(dynamic, dyn.mask, dangle, flip, fangles, TR,
                     r1=4, verbose=FALSE) {

  if (length(dim(flip)) != 4)  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  if (!is.logical(dyn.mask))  # Check dyn.mask is logical
    stop("Mask must be logical.")
  
  R1est <- R1.fast(flip, dyn.mask, fangles, TR, verbose)
  
  if (verbose)
    cat("  Calculating concentration...", fill=TRUE)
  theta <- dangle * pi/180
  B <- (1 - exp(-TR * R1est$R10)) / (1 - cos(theta) * exp(-TR * R1est$R10))
  A <- sweep(sweep(dynamic, 1:3, dynamic[,,,1], "-"),
             1:3, R1est$M0, "/") / sin(theta)
  R1t <- -(1/TR) * log((1 - sweep(A, 1:3, B, "+")) /
                       (1 - cos(theta) * sweep(A, 1:3, B, "+")))
  conc <- sweep(R1t, 1:3, R1est$R10, "-") / r1

  list(M0 = R1est$M0, R10 = R1est$R10, R1t = R1t, conc = conc)
}

#############################################################################
## setGeneric("CA.fast2")
#############################################################################

setGeneric("CA.fast2", function(dynamic, ...) standardGeneric("CA.fast2"))
setMethod("CA.fast2", signature(dynamic="array"),
	  function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
	      verbose=FALSE) 
	    dcemriWrapper("CA.fast2", dynamic, dyn.mask, dangle, flip, fangles,
		TR, r1, verbose))

#############################################################################
## CA.fast2()
#############################################################################

.CA.fast2 <- function(dynamic, dyn.mask, dangle, flip, fangles, TR, r1=4,
                     verbose=FALSE) {
  
  if (length(dim(flip)) != 4)  # Check flip is a 4D array
    stop("Flip-angle data must be a 4D array.")
  if (!is.logical(dyn.mask))  # Check dyn.mask is logical
    stop("Mask must be logical.")
  
  nangles <- length(fangles)
  nvoxels <- sum(dyn.mask)
  M <- nrow(flip)
  N <- ncol(flip)
  if(ntim(flip) != 2 || length(fangles) != 2)
    stop("Only two flip angles are allowed.")
  Z <- nsli(dynamic)
  W <- ntim(dynamic)
  
  if (verbose)
    cat("  Deconstructing data...", fill=TRUE)
  dyn.mat <- matrix(dynamic[dyn.mask], nvoxels)
  flip.mat <- matrix(flip[dyn.mask], nvoxels)
  R10 <- M0 <- numeric(nvoxels)

  if (verbose)
    cat("  Calculating R10 and M0...", fill=TRUE)
  for (k in 1:nvoxels) {
    x <- c(flip.mat[k,1] / tan(pi*fangles[1]/180),
           flip.mat[k,2] / tan(pi*fangles[2]/180))
    y <- c(flip.mat[k,1] / sin(pi*fangles[1]/180),
           flip.mat[k,2] / sin(pi*fangles[2]/180))
    fit <- lsfit(x, y)$coefficients
    R10[k] <- log(fit[2]) / -TR
    M0[k] <- fit[1] / (1 - fit[2])
  }

  if (verbose)
    cat("  Calculating concentration...", fill=TRUE)
  theta <- dangle * pi/180
  CD <- conc <- matrix(NA, nvoxels, W)
  B <- (1 - exp(-TR * R10)) / (1 - cos(theta) * exp(-TR * R10))
  A <- (dyn.mat - dyn.mat[,1]) / M0 / sin(theta)
  R1t <- -(1/TR) * log((1 - (A+B)) / (1 - cos(theta) * (A+B)))
  conc <- (R1t - R10) / r1

  if (verbose)
    cat("  Reconstructing results...", fill=TRUE)
  R10.array <- M0.array <- array(NA, c(M,N,Z,1))
  R10.array[dyn.mask] <- R10
  M0.array[dyn.mask] <- M0
  conc.array <- R1t.array <- array(NA, c(M,N,Z,W))
  ## conc.array <- array(NA, c(M,N,Z,W))
  mask4D <- array(dyn.mask, c(M,N,Z,W))
  conc.array[mask4D] <- unlist(conc)
  R1t.array[mask4D] <- unlist(R1t)
  ## list(M0 = M0.array, R10 = R10.array, conc = conc.array)
  list(M0 = M0.array, R10 = R10.array, R1t = R1t.array, conc = conc.array)
}
