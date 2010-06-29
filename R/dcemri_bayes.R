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
## $Id: dcemri_bayes.R 332 2010-01-29 16:54:07Z bjw34032 $
##
#############################################################################
## setGeneric("dcemri.bayes")
#############################################################################

setGeneric("dcemri.bayes",
           function(conc, ...) standardGeneric("dcemri.bayes"))
setMethod("dcemri.bayes", signature(conc="array"), 
          function(conc, time, img.mask, model="extended",
                         aif=NULL, user=NULL, nriters=3500, thin=3,
                         burnin=500, tune=267, tau.ktrans=1,
                         tau.kep=tau.ktrans, ab.vp=c(1,19),
                         ab.tauepsilon=c(1,1/1000), samples=FALSE,
                         multicore=FALSE, verbose=FALSE, ...) .dcemriWrapper("dcemri.bayes", conc, time, img.mask, model,
                         aif, user, nriters, thin,
                         burnin, tune, tau.ktrans,
                         tau.kep, ab.vp,
                         ab.tauepsilon, samples,
                         multicore, verbose, ...))


.dcemri.bayes.single <- function(conc, time, nriters=3500, thin=3,
                                burnin=500, tune=267, tau.gamma=1,
                                tau.theta=1, ab.vp=c(1,19),
                                ab.tauepsilon=c(1,1/1000), aif.model=0,
                                aif.parameter=c(2.4,0.62,3,0.016), vp=1) {

  if (sum(is.na(conc)) > 0)
    return(NA)
  else {
    n <- floor((nriters - burnin) / thin)
    if (tune > (0.5*nriters))
      tune <- floor(nriters/2)
    singlerun <- .C("dce_bayes_run_single",
                    as.integer(c(nriters, thin, burnin, tune)),
                    as.double(conc),
                    as.double(tau.gamma),
                    as.double(tau.theta),
                    as.double(ab.vp),
                    as.double(ab.tauepsilon),
                    as.double(c(aif.model, aif.parameter)),
                    as.integer(vp),
                    as.double(time),
                    as.integer(length(time)),
                    as.double(rep(0,n)),
                    as.double(rep(0,n)),
                    as.double(rep(0,n)),
                    as.double(rep(0,n)), 
                    PACKAGE="dcemriS4")    
    return(list("ktrans"= singlerun[[11]], "kep"= singlerun[[12]],
                "vp"= singlerun[[13]], "sigma2"= 1/singlerun[[14]]))
  }
}

.dcemri.bayes <- function(conc, time, img.mask, model="extended",
                         aif=NULL, user=NULL, nriters=3500, thin=3,
                         burnin=500, tune=267, tau.ktrans=1,
                         tau.kep=tau.ktrans, ab.vp=c(1,19),
                         ab.tauepsilon=c(1,1/1000), samples=FALSE,
                         multicore=FALSE, verbose=FALSE, ...) {

  ## dcemri.bayes - a function for fitting 1-compartment PK models to
  ## DCE-MRI images using Bayes inference
  ##
  ## authors: Volker Schmid, Brandon Whitcher
  ##
  ## input:
  ##        conc: array of Gd concentration,
  ##        time: timepoints of aquisition,
  ##        img.mask: array of voxels to fit,
  ##        D(=0.1): Gd dose in mmol/kg,
  ##        model: AIF... "weinman" or "parker",
  ##
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror), samples if samples=TRUE
  ##

  extract.samples <- function(sample, I, J, K, NRI) {
    A <- array(NA, c(I,J,K,NRI))
    count <- -1
    for (k in 1:K) {
      for (j in 1:J) {
        for (i in 1:I) {
          if (img.mask[i,j,k]) {
            count <- count + 1
            A[i,j,k,] <- sample[(1:NRI) + count*NRI]
          }
        }
      }
    }
    return(A)
  }

  if (nriters < thin)
    stop("Please check settings for nriters")
  if (burnin < 0)
    stop("Please check settings for burnin")
  if (thin < 1)
    stop("Please check settings for thin")
  if (tune < 50)
    stop("Please check settings for tune")

  if (burnin < tune) {
    burnin <- tune
    nriters <- nriters + tune
  } else {
    nriters <- nriters + burnin
  }
  
  aif <- switch(model,
                weinmann = ,
                extended = {
                  if (is.null(aif)) {
                    "tofts.kermode"
                  } else {
                    switch(aif,
                           tofts.kermode = "tofts.kermode",
                           fritz.hansen = "fritz.hansen",
                           stop("Only aif=\"tofts.kermode\" or aif=\"fritz.hansen\" acceptable AIFs for model=\"weinmann\" or model=\"extended\"", call.=FALSE))
                  }
                },
                orton.exp = {
                  if (is.null(aif)) {
                    "orton.exp"
                  } else {
                    switch(aif,
                           orton.exp = "orton.exp",
                           user = "user",
                           stop("Only aif=\"orton.exp\" or aif=\"user\" acceptable aifs for model=\"orton.exp\""), call.=FALSE)
                  }
                },
                stop("Unknown model: " + model, call.=FALSE))

  mod <- model
  nvoxels <- sum(img.mask)
  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)
  
  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2) {
      J <- K <- 1
    }
  }

  img.mask <- array(img.mask,c(I,J,K))
 
  if (verbose) cat("  Deconstructing data...", fill=TRUE)
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0

  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
           aif.parameter <- c(D*a1, m1, D*a2, m2)
         },
         fritz.hansen = {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
           aif.parameter <- c(D*a1, m1, D*a2, m2)
         },
         orton.exp = {
           D <- 1; a1 <- 323; m1 <- 20.2; a2 <- 1.07; m2 <- 0.172
           aif.parameter <- c(D*a1, m1, D*a2, m2)
         },
         ## FIXME orton.cos does not seem to be implemented
         ##orton.cos = {
         ##  D <- 1; a1 <- 2.84; m1 <- 22.8; a2 <- 1.36; m2 <- 0.171
         ##  aif.parameter=c(D*a1, m1, D*a2, m2)
         ##},
         user = {
           if (verbose) cat("  User-specified AIF parameters...", fill=TRUE);
           D <- try(user$D); AB <- try(user$AB) 
           muB <- try(user$muB); AG <- try(user$AG)
           muG <- try(user$muG)
           aif.parameter <- c(D * AB, muB, D * AG, muG)
           ## aG, aB are probably related to orton.cos which isn't implemented
           ## aG <- try(user$aG); aB <- try(user$aB);
         },
         print("WARNING: AIF parameters must be specified!"))

  ## translate "model" to "aif.model" and "vp.do"
  switch(model,
         weinmann = { aif.model <- 0 ; vp.do <- 0 },
         extended = { aif.model <- 0 ; vp.do <- 1 },
         orton.exp = { aif.model <- 1 ; vp.do <- 1 },
         stop("Model is not supported."))

  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  if (samples) {
    sigma2.samples <- ktrans.samples <- kep.samples <- c()
    if (mod %in% c("extended", "orton.exp", "orton.cos"))
      Vp.samples <- c()
  }

  if (verbose) cat("  Estimating the kinetic parameters...", fill=TRUE)

  conc.list <- list()
  for (i in 1:nvoxels)
    conc.list[[i]] <- conc.mat[i,]

  if (!multicore) {
    fit <- lapply(conc.list, FUN=.dcemri.bayes.single,
                  time=time, nriters=nriters, thin=thin, burnin=burnin,
                  tune=tune, ab.vp=ab.vp, ab.tauepsilon=ab.tauepsilon,
                  aif.model=aif.model, aif.parameter=aif.parameter,
		  vp=vp.do)
  } else {
    require("multicore")
    fit <- mclapply(conc.list, FUN=.dcemri.bayes.single,
                  time=time, nriters=nriters, thin=thin, burnin=burnin,
                  tune=tune, ab.vp=ab.vp, ab.tauepsilon=ab.tauepsilon,
                  aif.model=aif.model, aif.parameter=aif.parameter, 
                  vp=vp.do)
  }

  if (verbose) cat("  Reconstructing results...", fill=TRUE)

  for (k in 1:nvoxels) {
    ktrans$par[k] <- median(fit[[k]]$ktrans)
    kep$par[k] <- median(fit[[k]]$kep)
    ktrans$error[k] <- sqrt(var(fit[[k]]$ktrans))
    kep$error[k] <- sqrt(var(fit[[k]]$kep))
    if (mod %in% c("extended", "orton.exp", "orton.cos")) {
      Vp$par[k] <- median(fit[[k]]$vp)
      Vp$error[k] <- sqrt(var(fit[[k]]$vp))
      if (samples)
	Vp.samples <- c(Vp.samples,fit[[k]]$v)
    }
    sigma2[k] <- median(fit[[k]]$sigma2)
    if (samples) {
      ktrans.samples <- c(ktrans.samples,fit[[k]]$ktrans)
      kep.samples <- c(kep.samples,fit[[k]]$kep)
      sigma2.samples <- c(sigma2.samples,fit[[k]]$sigma2)
    }
  }

  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  ktrans.out <- list(par=A, error=B)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  kep.out <- list(par=A, error=B)

  if (mod %in% c("extended", "orton.exp", "orton.cos")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par 
    B[img.mask] <- Vp$error
    Vp.out <- list(par=A, error=B)
  }

  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2
  sigma2.out <- A

  if (samples) {
    NRI <- length(ktrans.samples) / length(ktrans$par)
    ktrans.out <- list(par=ktrans.out$par,
                       error=ktrans.out$error,
                       samples=extract.samples(ktrans.samples,I,J,K,NRI))
    kep.out <- list(par=kep.out$par,
                    error=kep.out$error,
                    samples=extract.samples(kep.samples,I,J,K,NRI))
    if (mod %in% c("extended", "orton.exp", "orton.cos"))
      Vp.out <- list(par=Vp.out$par,
                     error=Vp.out$error,
                     samples=extract.samples(Vp.samples, I, J, K, NRI))
    sigma2.samples <- extract.samples(sigma2.samples, I, J, K, NRI)
  }

  returnable <- list(ktrans=ktrans.out$par, kep=kep.out$par,
                     ktranserror=ktrans.out$error, keperror=kep.out$error, 
                     ve=ktrans.out$par/kep.out$par, sigma2=sigma2.out,
                     time=time)

  if (mod %in% c("extended", "orton.exp", "orton.cos")) {
    returnable[["vp"]] <- Vp.out$par
    returnable[["vperror"]] <- Vp.out$error
    if (samples) {
      returnable[["vp.samples"]] <- Vp.out$samples
    } 
  } 
  if (samples) {
    returnable[["ktrans.samples"]] <- ktrans.out$samples
    returnable[["kep.samples"]] <- kep.out$samples
    returnable[["sigma2.samples"]] <- sigma2.samples
  }
  return(returnable)
}

