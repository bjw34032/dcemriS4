##
##
## Copyright (c) 2009,2010 Brandon Whitcher and Volker Schmid
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
## $Id:$
##

#############################################################################
## T2.lm() = estimate exp(-TE/T2) using Levenburg-Marquardt
#############################################################################

T2.lm <- function(signal, TE, guess, nprint=0) {
  func <- function(x, y) {
    rho <- x[1]
    T2 <- x[2]
    signal <- y[[1]]
    TE <- y[[2]]
    signal - rho * exp(-TE/T2)
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=list(nprint=nprint, maxiter=150),
               y=list(signal, TE))
  list(rho=out$par[1], T2=out$par[2], info=out$info, message=out$message)
}

#############################################################################
## setGeneric("T2.fast")
#############################################################################

setGeneric("T2.fast", function(cpmg, ...) standardGeneric("T2.fast"))
setMethod("T2.fast", signature(cpmg="array"),
          function(cpmg, cpmg.mask, TE, verbose=FALSE) 
	    .dcemriWrapper("T2.fast", cpmg, cpmg.mask, TE, verbose))

#############################################################################
## T2.fast()
#############################################################################

.T2.fast <- function(cpmg, cpmg.mask, TE, verbose=FALSE) {

  if (length(dim(cpmg)) != 4) { # Check cpmg is a 4D array
    stop("CPMG data must be a 4D array.")
  }
  if (!is.logical(cpmg.mask)) { # Check cpmg.mask is logical
    stop("Mask must be logical.")
  }
    
  X <- nrow(cpmg)
  Y <- ncol(cpmg)
  Z <- nsli(cpmg)
  nvoxels <- sum(cpmg.mask)
  
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  cpmg.mat <- matrix(cpmg[cpmg.mask], nrow=nvoxels)
  T2 <- rho <- numeric(nvoxels)
  if (verbose) {
    cat("  Calculating T2 and rho...", fill=TRUE)
  }
  for (k in 1:nvoxels) {
    fit <- T2.lm(cpmg.mat[k,], TE, guess=c(0.75*cpmg.mat[k,1], 0.05))
    if (fit$info < 4) {
      T2[k] <- fit$T2
      rho[k] <- fit$rho
    } else {
      T2[k] <- rho[k] <- NA
    }
  }

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  T2.array <- rho.array <- array(NA, dim(cpmg)[1:3])
  T2.array[cpmg.mask] <- T2
  rho.array[cpmg.mask] <- rho

  list(rho = rho.array, T2 = T2.array)
}

