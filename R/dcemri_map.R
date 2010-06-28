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
## $Id: dcemri_map.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setGeneric("dcemri.map")
#############################################################################

setGeneric("dcemri.map", function(conc, ...) standardGeneric("dcemri.map"))
setMethod("dcemri.map", signature(conc="array"),
          function(conc, time, img.mask, model="extended",
                         aif=NULL, user=NULL, tau.ktrans=1,
                         tau.kep=tau.ktrans, ab.vp=c(1,19),
                         ab.tauepsilon=c(1,1/1000), samples=FALSE,
                         multicore=FALSE, verbose=FALSE, ...) .dcemriWrapper("dcemri.map", conc, time, img.mask, model,
                         aif, user, tau.ktrans,
                         tau.kep, ab.vp,
                         ab.tauepsilon, samples,
                         multicore, verbose, ...))


.dcemri.map.single <- function(conc, time, posterior, parameter,
                              transform, start, hyper, aif) {
  if (any(is.na(conc)))
    return(NA)
  else {    
    map <- optim(par=start, fn=posterior, conc=conc, time=time,
                 hyper=hyper, aif=aif, control=list("maxit"=25000))
    return.list <- list()
    for (i in 1:length(parameter))
      return.list[[parameter[i]]] <- transform[[i]](map$par[i])
    return(return.list)
  }
}

.dcemri.map <- function(conc, time, img.mask, model="extended",
                         aif=NULL, user=NULL, tau.ktrans=1,
                         tau.kep=tau.ktrans, ab.vp=c(1,19),
                         ab.tauepsilon=c(1,1/1000), samples=FALSE,
                         multicore=FALSE, verbose=FALSE, ...) {
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


	convterm <- function(kep, t, aif) {
	  aif[1] * (exp(-aif[2]*t) - exp(-kep*t)) / (kep-aif[2]) + aif[3] * (exp(-aif[4]*t) - exp(-kep*t)) / (kep-aif[4])
	}
	extraterm <- function(t, aif) {
	  aif[1] * (exp(-aif[2]*t)) + aif[3] * (exp(-aif[4]*t))
	}

  ## translate "model" to "aif.model" and "vp.do"
  switch(model,
         weinmann = {
           inverse <- function(x) 1/x
           parameter <- c("ktrans", "kep", "sigma2")
           transform <- c(exp, exp, inverse)
           start <- c(1, 1, 100)
           hyper <- c(tau.ktrans, tau.kep, ab.tauepsilon)
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             tauepsilon <- par[3]
             T <- length(time)
             p <- log(dnorm(gamma, 0, hyper[1]))
             p <- p + log(dnorm(theta, 0, hyper[2]))
             p <- p + log(dgamma(tauepsilon, hyper[3], rate=hyper[4]))
             conc.hat <- exp(gamma) * convterm(exp(theta), time, aif)
             p <- p + sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon))))
             if (is.na(p))
               p <- -1e-6
             return(-p)
           }
	 },		
         extended = {
           inverse <- function(x){return(1/x)}
           ident <- function(x){return(x)}
           parameter <- c("ktrans","kep","vp","sigma2")
           transform <- c(exp, exp, ident, inverse)
           start <- c(1, 1, 0.1, 100)
           hyper <- c(tau.ktrans, tau.kep, ab.vp, ab.tauepsilon)
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             vp <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             p <- log(dnorm(gamma, 0, hyper[1]))
             p <- p + log(dnorm(theta, 0, hyper[2]))
             p <- p + log(dgamma(tauepsilon, hyper[5], rate=hyper[6]))
             p <- p + log(dbeta(vp,hyper[3], hyper[4]))
             conc.hat <- (vp * extraterm(time, aif) +
                          exp(gamma) * convterm(exp(theta), time, aif))
             p <- p + sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon))))
             if (is.na(p))
               p <- -1e-6
             return(-p)
           }
	 },		
         stop("Model is not supported."))

  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  if (samples) {
    sigma2.samples <- ktrans.samples <- kep.samples <- c()
    if (mod %in% c("extended", "orton.exp", "orton.cos")) {
      Vp.samples <- c()
    }
  }

  if (verbose) cat("  Estimating the kinetic parameters...", fill=TRUE)

  conc.list <- list()
  for (i in 1:nvoxels)
    conc.list[[i]] <- conc.mat[i,]

  if (!multicore) {
    fit <- lapply(conc.list, FUN=.dcemri.map.single, time=time,
                  posterior=posterior, parameter=parameter,
		  transform=transform, start=start, hyper=hyper,
                  aif=aif.parameter)
  } else {
    require("multicore")
    fit <- mclapply(conc.list, FUN=.dcemri.map.single, time=time,
                    posterior=posterior, parameter=parameter,
                    transform=transform, start=start, hyper=hyper,
                    aif=aif.parameter)
  }

  if (verbose) cat("  Reconstructing results...", fill=TRUE)

  for (k in 1:nvoxels) {
    ktrans$par[k] <- fit[[k]]$ktrans
    kep$par[k] <- fit[[k]]$kep
    if (mod %in% c("extended", "orton.exp", "orton.cos")) {
      Vp$par[k] <- fit[[k]]$vp
    }
    sigma2[k] <- fit[[k]]$sigma2
  }

  A <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  ktrans.out <- list(par=A)
  A <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  kep.out <- list(par=A)

  if (mod %in% c("extended", "orton.exp", "orton.cos")) {
    A <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par 
    Vp.out <- list(par=A)
  }

  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2
  sigma2.out <- A

  returnable <- list(ktrans=ktrans.out$par, kep=kep.out$par,
                     ve=ktrans.out$par/kep.out$par, sigma2=sigma2.out,
                     time=time)

  if (mod %in% c("extended", "orton.exp", "orton.cos")) {
    returnable[["vp"]] <- Vp.out$par
  } 

  return(returnable)
}

