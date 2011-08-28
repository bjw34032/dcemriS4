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
## $Id: dcemri_map.R 332 2010-01-29 16:54:07Z bjw34032 $
##

#############################################################################
## setGeneric("dcemri.map")
#############################################################################

setGeneric("dcemri.map", function(conc, ...) standardGeneric("dcemri.map"))
setMethod("dcemri.map", signature(conc="array"),
          function(conc, time, img.mask, model="extended", aif=NULL,
                   user=NULL, ab.ktrans=c(0,1), ab.kep=ab.ktrans,
                   ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000), maxit=5000,
                   samples=FALSE, multicore=FALSE, verbose=FALSE, ...)
          .dcemriWrapper("dcemri.map", conc, time, img.mask, model,
                         aif, user, ab.ktrans, ab.kep, ab.vp,
                         ab.tauepsilon, maxit, samples, multicore, verbose,
                         ...))

.dcemri.map.single <- function(conc, time, posterior, parameter,
                               transform, start, hyper, aif, maxit,
                               verbose=FALSE) {
  if (any(is.na(conc))) {
    return(NA)
  } else {    
    map <- optim(par=start, fn=posterior, conc=conc, time=time,
                 hyper=hyper, aif=aif, control=list("maxit"=maxit))
    p <- length(parameter)
    return.list <- list()
    for (i in 1:p) {
      return.list[[parameter[i]]] <- transform[[i]](map$par[i])
    }
    return(return.list)
  }
}

.dcemri.map <- function(conc, time, img.mask, model="extended", aif=NULL,
                        user=NULL, ab.ktrans=c(0,1), ab.kep=ab.ktrans,
                        ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000),
                        maxit=5000, samples=FALSE, multicore=FALSE, 
                        verbose=FALSE, ...) {
  switch(model,
         weinmann = ,
         extended = {
           aif <- ifelse(is.null(aif), "tofts.kermode", aif)
           if (! aif %in% c("tofts.kermode","fritz.hansen"))
             stop("Only aif=\"tofts.kermode\" or aif=\"fritz.hansen\" are acceptable AIFs for model=\"weinmann\" or model=\"extended\"", call.=FALSE)
         },
         orton.exp = ,
         kety.orton.exp = {
           aif <- ifelse(is.null(aif), "orton.exp", aif)
           if (! aif %in% c("orton.exp","user"))
             stop("Only aif=\"orton.exp\" or aif=\"user\" are acceptable AIFs for model=\"orton.exp\" or model=\"kety.orton.exp\"", call.=FALSE)
         },
         orton.cos = ,
         kety.orton.cos = {
           aif <- ifelse(is.null(aif), "orton.cos", aif)
           if (! aif %in% c("orton.cos","user"))
             stop("Only aif=\"orton.cos\" or aif=\"user\" are acceptable AIFs for  model=\"orton.cos\" or model=\"kety.orton.cos\"", call.=FALSE)
         },
         stop(paste("Unknown model:",model), call.=FALSE))

  I <- nrow(conc)
  J <- ncol(conc)
  K <- nsli(conc)
  
  if (!is.numeric(dim(conc))) {
    I <- J <- K <- 1
  } else {
    if (length(dim(conc)) == 2) {
      J <- K <- 1
    }
    if (length(dim(conc)) == 3) {
      K <- 1
    }
  }

  aif.parameter <- aifParameters(aif, user)
  func.model <- compartmentalModel(model)
  inverse <- function(x) {
    return(1/x)
  }
  ident <- function(x) {
    return(x)
  }

  
  switch(model,
         weinmann =,
         kety.orton.exp = ,
         kety.orton.cos = {
           parameter <- c("ktrans", "kep", "sigma2")
           transform <- c(exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]), hyper[5] * hyper[6])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             tauepsilon <- par[3]
             conc.hat <- func.model(time, c(gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) +
                   log(dgamma(tauepsilon, hyper[5], rate=hyper[6])) +
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)

           }
	 },		
         extended = {
           inverse <- function(x) {
             1/x
           }
           ident <- function(x) {
             x
           }
           parameter <- c("ktrans", "kep", "vp", "sigma2")
           transform <- c(exp, exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.vp, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]),
                      log(hyper[5]/(hyper[5]+hyper[6])), hyper[7]*hyper[8])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             vp <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             conc.hat <- func.model(time, c(vp, gamma, theta), aif)
                         #ifelse(time > 0,
                         #       (exp(vp) * extraterm(time, aif) + exp(gamma) *
                         #        convterm(exp(theta), time, aif)),
                         #       0)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) +
                   log(dnorm(theta, hyper[3], hyper[4])) + 
                   log(dgamma(tauepsilon, hyper[7], rate=hyper[8])) + 
                   log(dbeta(exp(vp),hyper[5], hyper[6])) + 
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), -1e-6, p)
             return(-p)
           }
	 },		
         orton.exp = {
           inverse <- function(x) {
             1/x
           }
           ident <- function(x) {
             x
           }
           parameter <- c("ktrans", "kep", "vp", "sigma2")
           transform <- c(exp, exp, ident, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.vp, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]),
                      hyper[5]/(hyper[5]+hyper[6]), hyper[7]*hyper[8])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             vp <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             conc.hat <- func.model(time, c(vp, gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) + 
                   log(dnorm(theta, hyper[3], hyper[4])) + 
                   log(dgamma(tauepsilon, hyper[7], rate=hyper[8])) + 
                   log(dbeta(vp, hyper[5], hyper[6])) + 
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)
           }
	 },		
         kety.orton.exp = {
           inverse <- function(x) {
             1/x
           }
           parameter <- c("ktrans", "kep", "sigma2")
           transform <- c(exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.tauepsilon)
           start <- c(exp(hyper[1]), exp(hyper[3]), hyper[5]*hyper[6])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             tauepsilon <- par[3]
             T <- length(time)
             conc.hat <- func.model(time, c(gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) + 
                   log(dnorm(theta, hyper[3], hyper[4])) + 
                   log(dgamma(tauepsilon, hyper[5], rate=hyper[6])) + 
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)
           }
	 },		
        orton.cos = {
           inverse <- function(x) {
             1/x
           }
           ident <- function(x) {
             x
           }
           parameter <- c("ktrans", "kep", "vp", "sigma2")
           transform <- c(exp, exp, exp, inverse)
           hyper <- c(ab.ktrans, ab.kep, ab.vp, ab.tauepsilon)
           start <- c(hyper[1], hyper[3],
                      log(hyper[5] / (hyper[5] + hyper[6])),
                      hyper[7] * hyper[8])
           posterior <- function(par, conc, time, hyper, aif) {
             gamma <- par[1]
             theta <- par[2]
             theta0 <- par[3]
             tauepsilon <- par[4]
             T <- length(time)
             conc.hat <- func.model(time, c(theta0, gamma, theta), aif)
             p <- (log(dnorm(gamma, hyper[1], hyper[2])) + 
                   log(dnorm(theta, hyper[3], hyper[4])) + 
                   log(dgamma(tauepsilon, hyper[7], rate=hyper[8])) + 
                   log(dbeta(exp(theta0), hyper[5], hyper[6])) + 
                   sum(log(dnorm(conc, conc.hat, sqrt(1/tauepsilon)))))
             p <- ifelse(is.na(p), 1e-6, p)
             return(-p)
           }
	 },		
         stop("Model is not supported."))

  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  nvoxels <- sum(img.mask)
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels) # list()
  for (i in 1:nvoxels) {
    conc.list[[i]] <- conc.mat[i,]
  }

  if (verbose) {
    cat("  Estimating the kinetic parameters...", fill=TRUE)
  }
  if (multicore && require("multicore")) {
    fit <- mclapply(conc.list, FUN=.dcemri.map.single, time=time,
                    posterior=posterior, parameter=parameter,
                    transform=transform, start=start, hyper=hyper,
                    aif=aif.parameter, maxit=maxit, verbose=verbose)
  } else {
    fit <- lapply(conc.list, FUN=.dcemri.map.single, time=time,
                  posterior=posterior, parameter=parameter,
		  transform=transform, start=start, hyper=hyper,
                  aif=aif.parameter, maxit=maxit, verbose=verbose)
  }

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }
  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  for (k in 1:nvoxels) {
    try(ktrans$par[k] <- fit[[k]]$ktrans, silent=TRUE)
    try(kep$par[k] <- fit[[k]]$kep, silent=TRUE)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      try(vp$par[k] <- fit[[k]]$vp, silent=TRUE)
    }
    try(sigma2[k] <- fit[[k]]$sigma2, silent=TRUE)
  }
  A <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par

  returnable <- list(ktrans=A)
  A <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par

  returnable[["kep"]] <- A
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    A <- array(NA, c(I,J,K))
    A[img.mask] <- vp$par
    returnable[["vp"]] <- A
  }
  A <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2
  
  returnable[["sigma2"]] <- A
  returnable[["time"]] <- time
  returnable[["ve"]] <- returnable$ktrans / returnable$kep

  return(returnable)
}
