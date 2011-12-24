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
## $Id: dwi.R 332 2010-01-29 16:54:07Z bjw34032 $
## 

#############################################################################
## adc.lm() = estimate ADC using Levenburg-Marquardt
#############################################################################

adc.lm <- function(signal, b, guess, control=nls.lm.control()) {
  func <- function(x, y) {
    S0 <- x[1]
    D <- x[2]
    signal <- y[[1]]
    b <- y[[2]]
    signal - S0 * exp(-b*D)
  }
  require("minpack.lm") # Levenberg-Marquart fitting
  out <- nls.lm(par=guess, fn=func, control=control, y=list(signal, b))
  list(S0=out$par[1], D=out$par[2], hessian=out$hessian, info=out$info,
       message=out$message)
}

#############################################################################
## setGeneric("ADC.fast")
#############################################################################

setGeneric("ADC.fast", function(dwi, ...) standardGeneric("ADC.fast"))
setMethod("ADC.fast", signature(dwi="array"),
          function(dwi, bvalues, dwi.mask,
                   control=nls.lm.control(maxiter=150),
                   multicore=FALSE, verbose=FALSE)
          .dcemriWrapper("ADC.fast", dwi, bvalues, dwi.mask, control,
                         multicore, verbose))

.ADC.fast <- function(dwi, bvalues, dwi.mask,
                      control=nls.lm.control(maxiter=150),
                      multicore=FALSE, verbose=FALSE) {
  if (length(dim(dwi)) != 4) { # Check dwi is a 4D array
    stop("Diffusion-weighted data must be a 4D array.")
  }
  if (!is.logical(dwi.mask)) { # Check dyn.mask is logical
    stop("Mask must be logical.")
  }
  nvalues <- length(bvalues)
  nvoxels <- sum(dwi.mask)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  dwi.mat <- matrix(dwi[dwi.mask], nvoxels)
  dwi.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    dwi.list[[k]] <- dwi.mat[k,]
  }
  if (verbose) {
    cat("  Calculating S0 and D...", fill=TRUE)
  }
  if (multicore && require("parallel")) {
    fit.list <- mclapply(dwi.list, function(x) {
      adc.lm(x, bvalues, guess=c(0.75*x[1], 0.001), control)
    })
  } else {
    fit.list <- lapply(dwi.list, function(x) {
      adc.lm(x, bvalues, guess=c(0.75*x[1], 0.001), control)
    })
  }
  rm(dwi.list) ; gc()
  S0 <- D <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  for (k in 1:nvoxels) {
    if (fit.list[[k]]$info > 0 && fit.list[[k]]$info < 5) {
      S0$par[k] <- fit.list[[k]]$S0
      D$par[k] <- fit.list[[k]]$D
      S0$error[k] <- sqrt(fit.list[[k]]$hessian[1,1])
      D$error[k] <- sqrt(fit.list[[k]]$hessian[2,2])
    } else {
      S0$par[k] <- D$par[k] <- S0$error[k] <- D$error[k] <- NA
    }
  }
  rm(fit.list) ; gc()
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  S0.array <- D.array <- S0error <- Derror <- array(NA, dim(dwi)[1:3])
  S0.array[dwi.mask] <- S0$par
  D.array[dwi.mask] <- D$par
  S0error[dwi.mask] <- S0$error
  Derror[dwi.mask] <- D$error
  
  list(S0 = S0.array, D = D.array, S0.error = S0error, D.error = Derror)
}
