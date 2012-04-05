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
##
## $Id: dcemri_lm.R 329 2010-01-07 16:33:56Z bjw34032 $
##

#############################################################################
## setGeneric("dcemri.lm")
#############################################################################

setGeneric("dcemri.lm",
           function(conc,  ...) standardGeneric("dcemri.lm"))
setMethod("dcemri.lm", signature(conc="array"), 
	  function(conc,time,img.mask, model="extended", aif=NULL,
                   control=nls.lm.control(), user=NULL, guess=NULL,
                   multicore=FALSE, verbose=FALSE, ...)
          .dcemriWrapper("dcemri.lm", conc, time, img.mask, model, aif,
                         control, user, guess, multicore, verbose, ...))

.dcemri.lm <- function(conc, time, img.mask, model="extended", aif=NULL,
                       control=nls.lm.control(), user=NULL, guess=NULL,
                       multicore=FALSE, verbose=FALSE, ...) {
  require("minpack.lm")
  switch(model,
         weinmann = ,
         extended = {
           aif <- ifelse(is.null(aif), "tofts.kermode", aif)
           aif.names <- c("tofts.kermode","fritz.hansen","empirical")
           if (! aif %in% aif.names) {
             stop(sprintf("Only aif=\"%s\" or aif=\"%s\" or aif=\"%s\" are acceptable AIFs for model=\"weinmann\" or model=\"extended\"", aif.names[1], aif.names[2], aif.names[3]), call.=FALSE)
           }
         },
         kety.orton.exp = ,
         orton.exp = {
           aif <- ifelse(is.null(aif), "orton.exp", aif)
           if (! aif %in% c("orton.exp","user")) {
             stop("Only aif=\"orton.exp\" or aif=\"user\" are acceptable AIFs for model=\"orton.exp\" or model=\"kety.orton.exp\"", call.=FALSE)
           }
         },
         kety.orton.cos= ,
         orton.cos = {
           aif <- ifelse(is.null(aif), "orton.cos", aif)
           if (! aif %in% c("orton.cos","user")) {
             stop("Only aif=\"orton.cos\" or aif=\"user\" are acceptable AIFs for model=\"orton.cos\" or model=\"kety.orton.cos\"", call.=FALSE)
           }
         },
         stop(paste("Unknown model:", model), call.=FALSE))
  p <- aifParameters(aif, user)
  if (aif == "empirical") {
    ## paste model and "empirical" together
    model <- paste(model, aif, sep=".")
  }
  func.model <- compartmentalModel(model)
  func <- function(theta, signal, time, ...) {
    out <- signal - func.model(time, theta, p)
    out[!is.na(out)]
  }
  nvoxels <- sum(img.mask)
  switch(model,
         weinmann = ,
         weinmann.empirical = ,
         kety.orton.exp = ,
         kety.orton.cos = {
           if (is.null(guess)) {
             guess <- c("th1"=-1, "th3"=-1)
           } else {
             if (length(guess) != 2 || !all(names(guess) %in% c("th1","th3"))) {
               stop("Names of starting parameters must be \"th1\" and \"th3\"")
             }
           }
         },
         extended = ,
         extended.empirical = ,
         orton.exp = ,
         orton.cos = {
           if (is.null(guess)) {
             guess <- c("th0"=-3, "th1"=-1, "th3"=-1)
             vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
           } else {
             if (length(guess) != 3 ||
                 all(names(guess) %in% c("th0","th1","th3"))) {
               stop("Names of starting parameters must be \"th0\", \"th1\" and \"th3\"")
             }
           }
         },
         stop("Model/AIF combination is not supported."))
       
  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)
  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2)
      J <- K <- 1
  }
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels)
  for (k in 1:nvoxels) {
    conc.list[[k]] <- conc.mat[k,]
  }
  rm(conc.mat) ; gc()
  if (verbose) {
    cat("  Estimating the kinetic parameters...", fill=TRUE)
  }
  if (multicore && require("parallel")) {
    lm.list <- mclapply(conc.list, function(x) {
      nls.lm(par=guess, fn=func, control=control, signal=x, time=time, p=p)
    })
  } else {
    lm.list <- lapply(conc.list, function(x) {
      nls.lm(par=guess, fn=func, control=control, signal=x, time=time, p=p)
    })
  }
  rm(conc.list) ; gc()
  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sse <- rep(NA, nvoxels)
  for (k in 1:nvoxels) {
    if (lm.list[[k]]$info > 0 && lm.list[[k]]$info < 5) {
      ktrans$par[k] <- exp(lm.list[[k]]$par["th1"])
      kep$par[k] <- exp(lm.list[[k]]$par["th3"])
      ktrans$error[k] <- sqrt(lm.list[[k]]$hessian["th1","th1"])
      kep$error[k] <- sqrt(lm.list[[k]]$hessian["th3","th3"])
      sse[k] <- lm.list[[k]]$deviance
      if (model %in% c("extended", "orton.exp", "orton.cos",
                       "extended.empirical")) {
        vp$par[k] <- exp(lm.list[[k]]$par["th0"])
        vp$error[k] <- sqrt(lm.list[[k]]$hessian["th0","th0"])
      }
    }
  }
  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  R <- list(ktrans=A, ktranserror=B, time=time)
  rm(A,B)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  R$kep <- A
  R$keperror <- B
  R$ve <- R$ktrans / R$kep
  rm(A,B)
  if (model %in% c("extended", "orton.exp", "orton.cos",
                   "extended.empirical")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- vp$par
    B[img.mask] <- vp$error
    R$vp <- A
    R$vperror <- B
    rm(A,B)
  }
  A <- array(NA, c(I,J,K))
  A[img.mask] <- sse
  R$sse <- A
  rm(A)
  return(R)
}


