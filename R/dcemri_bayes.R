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
                         burnin=500, tune=267, ab.ktrans=c(0,1),
                         ab.kep=ab.ktrans, ab.vp=c(1,19),
                         ab.tauepsilon=c(1,1/1000), samples=FALSE,
                         multicore=FALSE, verbose=FALSE, dic=FALSE, ...)
          .dcemriWrapper("dcemri.bayes", conc, time, img.mask, model,
                        aif, user, nriters, thin, burnin, tune, ab.ktrans,
                        ab.kep, ab.vp, ab.tauepsilon, samples,
                        multicore, verbose, dic, ...))

.dcemri.bayes.single <- function(conc, time, nriters=3500, thin=3,
                                 burnin=500, tune=267, ab.gamma=c(0,1),
                                 ab.theta=c(0,1), ab.vp=c(1,19),
                                 ab.tauepsilon=c(1,1/1000), aif.model=0,
                                 aif.parameter=c(2.4,0.62,3,0.016), vp=1) {
  
  if (sum(is.na(conc)) > 0) {
    return(NA)
  } else {
    n <- floor((nriters - burnin) / thin)
    if (tune > nriters/2) {
      tune <- floor(nriters/2)
    }
    n0 <- rep(0, n)
    singlerun <- .C("dce_bayes_run_single",
                    as.integer(c(nriters, thin, burnin, tune)),
                    as.double(conc),
                    as.double(ab.gamma),
                    as.double(ab.theta),
                    as.double(ab.vp),
                    as.double(ab.tauepsilon),
                    as.double(c(aif.model, aif.parameter)),
                    as.integer(vp), # is this correct?
                    as.double(time),
                    as.integer(length(time)),
                    as.double(n0),
                    as.double(n0),
                    as.double(n0),
                    as.double(n0), 
                    as.double(n0), 
                    PACKAGE="dcemriS4")    
    list("ktrans"= singlerun[[11]],
         "kep"= singlerun[[12]],
         "vp"= singlerun[[13]],
         "sigma2"= 1/singlerun[[14]],
         "deviance" = singlerun[[15]])
  }
}

.dcemri.bayes <- function(conc, time, img.mask, model="extended",
                          aif=NULL, user=NULL, nriters=3500, thin=3,
                          burnin=500, tune=267, ab.ktrans=c(0,1),
                          ab.kep=ab.ktrans, ab.vp=c(1,19),
                          ab.tauepsilon=c(1,1/1000), samples=FALSE,
                          multicore=FALSE, verbose=FALSE, dic=FALSE,
                          ...) {
  
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
  
  extractSamples <- function(sample, img.mask, NRI) {
    out <- array(NA, c(NRI, dim(img.mask)))
    out[img.mask] <- sample
    aperm(out, c(2:length(dim(out)), 1)) # not too sure about drop()
  }
  
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
  
  if (J > 1 && K > 1) {
    if (sum(dim(img.mask) - dim(conc)[-length(dim(conc))]) != 0) {
      stop("Dimensions of \"conc\" do not agree with \"img.mask\"")
    }
  }
  if (nriters < thin) {
    stop("Please check settings for nriters")
  }
  if (burnin < 0) {
    stop("Please check settings for burnin")
  }
  if (thin < 1) {
    stop("Please check settings for thin")
  }
  if (tune < 50) {
    stop("Please check settings for tune")
  }
  if (burnin < tune) {
    burnin <- tune
    nriters <- nriters + tune
  } else {
    nriters <- nriters + burnin
  }
  
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
  am <- grep("^[Aa]|^[Mm][^Ee]", names(p))
  aif.parameter <- unlist(p[am])
  if (!is.null(p$D) && p$D != 1) {
    a <- grep("^[Aa]", names(p))
    aif.parameter[a] <- p$D * aif.parameter[a]
  }
  
  ## translate "model" to "aif.model" and "vp.do"
  switch(model,
         weinmann = {
           aif.model <- 0
           vp.do <- 0
         },
         extended = {
           aif.model <- 0
           vp.do <- 1
         },
         orton.exp = {
           aif.model <- 1
           vp.do <- 1
         },
         kety.orton.exp = {
           aif.model <- 1
           vp.do <- 0
         },
         stop("Model is not supported."))

  ## img.mask <- array(img.mask,c(I,J,K)) # why?
  nvoxels <- sum(img.mask)
  if (verbose) {
    cat("  Deconstructing data...", fill=TRUE)
  }
  img.mask <- ifelse(img.mask > 0, TRUE, FALSE)
  conc.mat <- matrix(conc[img.mask], nvoxels)
  conc.mat[is.na(conc.mat)] <- 0
  conc.list <- vector("list", nvoxels)
  for (i in 1:nvoxels) {
    conc.list[[i]] <- conc.mat[i,]
  }
  if (verbose) {
    cat("  Estimating the kinetic parameters...", fill=TRUE)
  }
  if (multicore && require("multicore")) {
    bayes.list <- mclapply(conc.list, FUN=.dcemri.bayes.single,
                           time=time, nriters=nriters, thin=thin,
                           burnin=burnin, tune=tune, ab.gamma=ab.ktrans,
                           ab.theta=ab.kep, ab.vp=ab.vp,
                           ab.tauepsilon=ab.tauepsilon, aif.model=aif.model,
                           aif.parameter=aif.parameter, vp=vp.do)
  } else {
    bayes.list <- lapply(conc.list, FUN=.dcemri.bayes.single, time=time,
                         nriters=nriters, thin=thin, burnin=burnin,
                         tune=tune, ab.gamma=ab.ktrans, ab.theta=ab.kep,
                         ab.vp=ab.vp, ab.tauepsilon=ab.tauepsilon,
                         aif.model=aif.model, aif.parameter=aif.parameter,
                         vp=vp.do)    
  }
  rm(conc.list) ; gc()
  
  if (verbose) {
    cat("  Extracting results...", fill=TRUE)
  }
  
  n <- (nriters - burnin) / thin # number of samples from posterior
  ktrans <- kep <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  sigma2 <- rep(NA, nvoxels)
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    Vp <- list(par=rep(NA, nvoxels), error=rep(NA, nvoxels))
  }
  if (dic) {
    med.deviance <- rep(NA, nvoxels)
  }
  if (samples) {
    sigma2.samples <- ktrans.samples <- kep.samples <- rep(NA, n*nvoxels)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      Vp.samples <- rep(NA, n*nvoxels)
    }
    if (dic) {
      deviance.samples <- rep(NA, n*nvoxels)
    }
  }

  for (k in 1:nvoxels) {
    index <- nvoxels * (k-1) + (1:n)
    ktrans$par[k] <- median(bayes.list[[k]]$ktrans)
    kep$par[k] <- median(bayes.list[[k]]$kep)
    ktrans$error[k] <- sd(bayes.list[[k]]$ktrans)
    kep$error[k] <- sd(bayes.list[[k]]$kep)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      Vp$par[k] <- median(bayes.list[[k]]$vp)
      Vp$error[k] <- sd(bayes.list[[k]]$vp)
      if (samples) {
        Vp.samples[index] <- bayes.list[[k]]$v
      }
    }
    sigma2[k] <- median(bayes.list[[k]]$sigma2)
    if (dic) {
      med.deviance[k] <- median(bayes.list[[k]]$deviance)
    }
    if (samples) {
      ktrans.samples[index] <- bayes.list[[k]]$ktrans
      kep.samples[index] <- bayes.list[[k]]$kep
      sigma2.samples[index] <- bayes.list[[k]]$sigma2
      if (dic) {
        deviance.samples[index] <- bayes.list[[k]]$deviance
      }
    }
  }
  rm(bayes.list) ; gc()

  if (verbose) {
    cat("  Reconstructing results...", fill=TRUE)
  }

  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- ktrans$par
  B[img.mask] <- ktrans$error
  ktrans.out <- list(par=A, error=B)
  rm(A,B,ktrans)
  A <- B <- array(NA, c(I,J,K))
  A[img.mask] <- kep$par
  B[img.mask] <- kep$error
  kep.out <- list(par=A, error=B)
  rm(A,B,kep)
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    A <- B <- array(NA, c(I,J,K))
    A[img.mask] <- Vp$par 
    B[img.mask] <- Vp$error
    Vp.out <- list(par=A, error=B)
    rm(A,B,Vp)
  }

  A <- array(NA, c(I,J,K))
  A[img.mask] <- sigma2
  sigma2.out <- A
  rm(A,sigma2)

  if (dic) {
    A <- array(NA, c(I,J,K))
    A[img.mask] <- med.deviance
    med.deviance <- A
    rm(A)
  }

  if (samples) {
    if (verbose) {
      cat("  Reconstructing samples...", fill=TRUE)
    }
    ktrans.out$samples <- extractSamples(ktrans.samples, img.mask, n)
    kep.out$samples <- extractSamples(kep.samples, img.mask, n)
    if (model %in% c("extended", "orton.exp", "orton.cos")) {
      Vp.out$samples <- extractSamples(Vp.samples, img.mask, n)
    }
    sigma2.samples <- extractSamples(sigma2.samples, img.mask, n)
    if (dic) {
      deviance.samples <- extractSamples(deviance.samples, img.mask, n)
    }
  }

  returnable <- list(ktrans=ktrans.out$par,
                     kep=kep.out$par,
                     ktranserror=ktrans.out$error,
                     keperror=kep.out$error,
                     ve=ktrans.out$par/kep.out$par,
                     time=time)
  if (model %in% c("extended", "orton.exp", "orton.cos")) {
    returnable[["vp"]] <- Vp.out$par
    returnable[["vperror"]] <- Vp.out$error
    if (samples) {
      returnable[["vp.samples"]] <- Vp.out$samples
    } 
  } 
  ## DIC
  if (dic) {
    if (verbose) {
      cat("  Computing DIC...", fill=TRUE)
    }
    fitted <- array(NA, c(I,J,K,length(time)))
    for (i in 1:I) {
      for (j in 1:J) {
        for (k in 1:K) {
          if (img.mask[i,j,k]) {
            par <- list("ktrans"=ktrans.out$par[i,j,k],
                        "kep"=kep.out$par[i,j,k])
            if (vp.do) {
              par["vp"] <- Vp.out$par[i,j,k]
            }
            fitted[i,j,k,] <- kineticModel(time, par, model=model, aif=aif)
          }
        }
      }
    }
    conc <- array(conc, c(I,J,K,length(time)))
    fitted <- fitted - conc
    fitted <- fitted * fitted
    fitted <- apply(fitted, 1:3, sum)
    deviance.med <- length(time) * log(sigma2.out) + fitted / sigma2.out
    pD <- med.deviance - deviance.med
    DIC <- med.deviance + pD
    returnable[["DIC"]] <- sum(DIC,na.rm=TRUE)
    returnable[["pD"]] <- sum(pD,na.rm=TRUE)
    returnable[["DIC.map"]] <- DIC
    returnable[["pD.map"]] <- pD 
    returnable[["deviance.med"]] <- deviance.med
    returnable[["med.deviance"]] <- med.deviance
    if (samples) {
      returnable[["deviance.samples"]] <- deviance.samples
    }
    rm(DIC)
    rm(pD)
    rm(deviance.med)
    rm(med.deviance)
    rm(deviance.samples)
    gc()
  }
  if (samples) {
    temp <- ktrans.out$samples
    rm(ktrans.out)
    returnable[["ktrans.samples"]] <- temp
    temp <- kep.out$samples
    rm(kep.out)
    returnable[["kep.samples"]] <- temp
    returnable[["sigma2.samples"]] <- sigma2.samples
  }
  ## rm(Vp.out) ; gc()
  
  return(returnable)
}

