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

T2.lm <- function(signal, TE, guess, control=minpack.lm::nls.lm.control()) {
  func <- function(x, y) {
    rho <- x[1]
    T2 <- x[2]
    signal <- y[[1]]
    TE <- y[[2]]
    signal - rho * exp(-TE/T2)
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=control, 
                            y=list(signal, TE))
  list(rho=out$par[1], T2=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}

#############################################################################
## setGeneric("T2.fast")
#############################################################################

setGeneric("T2.fast", function(cpmg, ...) standardGeneric("T2.fast"))
setMethod("T2.fast", signature(cpmg="array"),
          function(cpmg, cpmg.mask, TE, 
                   control=minpack.lm::nls.lm.control(maxiter=150),
                   multicore=FALSE, verbose=FALSE)
          .dcemriWrapper("T2.fast", cpmg, cpmg.mask, TE, control, multicore,
                         verbose))

#############################################################################
## T2.fast()
#############################################################################

.T2.fast <- function(cpmg, cpmg.mask, TE, 
                     control=minpack.lm::nls.lm.control(maxiter=150),
                     multicore=FALSE, verbose=FALSE) {

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
  cpmg.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    cpmg.list[[k]] <- cpmg.mat[k,]
  }
  if (verbose) {
    cat("  Calculating T2 and rho...", fill=TRUE)
  }
  if (multicore && require("parallel")) {
    T2.list <- mclapply(cpmg.list, function(x) {
      T2.lm(x, TE, guess=c(0.75 * x[1], 0.05), control)
    })
  } else {
    T2.list <- lapply(cpmg.list, function(x) {
      T2.lm(x, TE, guess=c(0.75 * x[1], 0.05), control)
    })
  }
  rm(cpmg.list) ; gc()
  T2 <- rho <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  pb <- txtProgressBar()
  for (k in 1:nvoxels) {
    if (T2.list[[k]]$info > 0 && T2.list[[k]]$info < 5) {
      T2$par[k] <- T2.list[[k]]$T2
      rho$par[k] <- T2.list[[k]]$rho
      T2$error[k] <- sqrt(T2.list[[k]]$hessian[1,1])
      rho$error[k] <- sqrt(T2.list[[k]]$hessian[2,2])
    } else {
      T2$par[k] <- rho$par[k] <- T2$error[k] <- rho$error[k] <- NA
    }
    setTxtProgressBar(pb, k)
  }
  close(pb)
  rm(T2.list) ; gc()
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  T2.array <- rho.array <- T2.error <- rho.error <- array(NA, dim(cpmg)[1:3])
  T2.array[cpmg.mask] <- T2$par
  rho.array[cpmg.mask] <- rho$par
  T2.error[cpmg.mask] <- T2$error
  rho.error[cpmg.mask] <- rho$error

  list(rho = rho.array, T2 = T2.array, rho.error = rho.error,
       T2.error = T2.error)
}

