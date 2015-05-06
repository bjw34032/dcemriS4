##
##
## Copyright (c) 2011 Volker Schmid, Julia Sommer
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

#############################################################################
## setGeneric("dcemri.space")
#############################################################################

setGeneric("dcemri.space",
           function(conc, ...) standardGeneric("dcemri.space"))
setMethod("dcemri.space", signature(conc="array"), 
          function(conc, time, img.mask, starting=0,
				  model="extended",
				  aif=NULL,  
          nriters=5000,
				  burnin=500, 
				  tuning=200,
				  thin=10,
				  ab.gamma=c(0.0001,0.0001), 
				  ab.gamma3=3*c(0.0001,0.0001),					  
				  ab.theta=c(0.0001,0.0001),
				  ab.theta3=3*c(0.0001,0.0001),
				  ab.vp=c(1,19),
				  ab.tauepsilon=c(1,1/1000), 
          spatial=0,
				  slice=0,
				  gemupdate=0,
				  uptauep=0,
				  retunecycles=3,
				  tunepct=5,
				  t0=0,
				  samples=FALSE,
				  multicore=FALSE, 
				  verbose=FALSE, 
				  dic=FALSE,
				  user=NULL)
          .dcemriWrapper("dcemri.space", conc, time, img.mask, starting,
				  model,
				  aif,  
				  nriters,
				  burnin, 
				  tuning,
				  thin,
				  ab.gamma, 
				  ab.gamma3,					  
				  ab.theta,
				  ab.theta3,  
				  ab.vp,
				  ab.tauepsilon, 
				  spatial,
				  slice,
				  gemupdate,
				  uptauep,
				  retunecycles,
				  tunepct,
				  t0,
				  samples,
				  multicore, 
				  verbose, 
				  dic,
				  user))

.dcemri.space <- function(conc, time, img.mask, starting=0,
				model="extended",
				aif=NULL,  
				nriters=5000,
				burnin=500, 
				tuning=200,
				thin=10,
				ab.gamma=c(1,0.0001), 
				ab.gamma3=3*c(1,0.0001),					  
				ab.theta=c(1,0.0001),
				ab.theta3=3*c(1,0.0001),
				ab.vp=c(1,19),
				ab.tauepsilon=c(1,1/1000), 
				spatial=0,
				slice=0,
				gemupdate=0,
				uptauep=0,
				retunecycles=3,
				tunepct=5,
				t0=0,
				samples=FALSE,
				multicore=FALSE, 
				verbose=FALSE, 
				dic=FALSE,
				user=NULL) {
  ## dcemri.space - a function for fitting 1-compartment PK models 
  ## with spatial prior to DCE-MRI images using Bayes inference
  ## 
  ## authors: Volker Schmid, Julia Sommer
  ##
  ## input:
  ##        conc: array of Gd concentration,
  ##        time: timepoints of aquisition,
  ##        img.mask: array of voxels to fit,
  ##        D(=0.1): Gd dose in mmol/kg,
  ##        model: AIF... "weinman" or "parker",
  ##  
  ##		  ab.taugamma
  ##		  ab.tautheta
  ##		  ab.gamma  : precision for voxels in same slice
  ##		  ab.gamma3 : precision for voxels in different slices
  ##		  ab.theta  : precision for voxels in same slice
  ##		  ab.theta3	: precision for voxels in different slices
  ##		  ab.vp
  ##		  ab.tauepsilon
  ##		  aif.model 
  ##		  aif.parameter
  ##		  settings: 
  ##			spatial 0: no spatial model, fit voxelwise
  ##					1: 2D spatial model, fit layers separately
  ##					3: 3D spatial model
  ##			slice	0: all slices 
  ##					k: restrict to kths slice
  ##			gemupdate: 0: default
  ##					   3: same weights in GMRF for ktrans and kep
  ##			uptauep: 0: fixed precisions 
  ##					 1: estimate precision of error
  ##	 				 2: 1+global smoothing
  ##					 3: 1+adaptive smoothing
  ##			tuning: tune after n iterations
  ##			tunepct: percentage of voxels which can fail tuning
  ##			retunecycles: number of times to do tuning
  ##		  
  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
  ##         (ktranserror and keperror), samples if samples=TRUE
  ##  

  if (spatial==2)spatial=1
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
  # a1 and a2 have to be multiplied with D
  if (!is.null(p$D) && p$D != 1) {
    a <- grep("^[Aa]", names(aif.parameter))
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
           vp.do <- 3 # !!
         },
         orton.exp = {
           aif.model <- 1
           vp.do <- 3 # !!
         },
         kety.orton.exp = {
           aif.model <- 1
           vp.do <- 0
         },
         stop("Model is not supported."))

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
  
  img.mask <- ifelse(img.mask > 0, 1, 0)
  img.mask <- array(img.mask,c(I,J,K)) ## convert array in arbitrary dimension to 3D array
  
  T<-length(time)

  x <- which(apply(img.mask,1,sum)>0)
  y <- which(apply(img.mask,2,sum)>0)
 
  img.mask <- img.mask[x,y,]
  conc <- conc[x,y,,]
  
  N <- prod(dim(img.mask))
  N1 <- sum(img.mask)
  
  if (verbose) {
    cat(paste(N1,"voxels in mask."))
    cat(" Tuning MCMC algorithm 0 %")
  }

  II <- nrow(conc)
  JJ <- ncol(conc)
  KK <- nsli(conc)
  
  if (!is.numeric(dim(conc))) {
    II <- JJ <- KK <- 1
  } else {
    if (length(dim(conc)) == 2) {
      JJ <- KK <- 1
    }
    if (length(dim(conc)) == 3) {
      KK <- 1
    }
  }

  conc[is.na(conc)]<-0
  conc<- array(conc,c(II,JJ,KK,T)) ## convert array in arbitrary dimension to 4D array
  img.mask<- array(img.mask,c(II,JJ,KK)) ## convert array in arbitrary dimension to 3D array
 
  # starting values
	if(is.list(starting)==FALSE) # default starting values
	{
		kep_1 <- 0.5  
		ve_1 <- 0.5
		ktrans_1 <- kep_1 * ve_1
	
		ktrans1_start <- rep(ktrans_1,N)
		kep1_start <- rep(kep_1,N)
		
	}
	else 
	{			
		ktrans1_start <- starting$ktrans
		kep1_start <- starting$kep
		
		ktrans1_start[is.na(ktrans1_start)] <- 0
		kep1_start[is.na(kep1_start)] <- 0
		
	}
   
  # tuning step
  singlerun <- .C("dce_space",
                  as.integer(c(tuning+2, tuning+1, tuning)),
                  as.double(conc),
                  as.integer(img.mask),
                  as.integer(dim(conc)),
                  as.double(ab.gamma), #5
                  as.double(ab.gamma3),
                  as.double(ab.theta),
                  as.double(ab.theta3),
                  as.double(ab.vp),
                  as.double(ab.tauepsilon),#10
                  as.double(c(aif.model, aif.parameter)),
                  as.integer(c(vp.do,spatial,slice,gemupdate,uptauep,retunecycles,tunepct)),
				  as.double(c(0,time)), # time shifted, s.t. it fits with conc
                  as.double(rep(0,N)),					
                  as.double(rep(ab.gamma[1]/ab.gamma[2],N)), # 15
                  as.double(rep(ab.theta[1]/ab.theta[2],N)), 
				  as.double(ktrans1_start), # ktrans
				  as.double(kep1_start), # kep
                  as.double(if(vp.do==3){rep(0.1,N)}else{rep(0,N)}), # vp
                  as.double(rep(ab.tauepsilon[1]/ab.tauepsilon[2],N)), # tau_epsilon 20
                  as.double(rep(.5,N)), 
                  as.double(rep(0.5,N)), 
                  as.double(rep(0.5,N)), 
                  as.double(rep(0.5,N)), #
                  as.integer(rep(0,N)), #25
                  as.integer(rep(0,N)),
                  as.integer(rep(0,N)),
                  as.integer(rep(0,N)),
                  as.double(rep(ab.gamma[1]/ab.gamma[2],N)), 
                  as.double(rep(ab.theta[1]/ab.theta[2],N)), # 30
                  as.double(rep(ab.gamma3[1]/ab.gamma3[2],N)), 
                  as.double(rep(ab.theta3[1]/ab.theta3[2],N)), 
                  as.double(rep(ab.vp[1]/ab.vp[2],N)), 
                  as.double(rep(ab.vp[1]/ab.vp[2],N)), 
                  as.double(rep(ab.vp[1]/ab.vp[2],N)), #35
                  as.double(rep(ab.vp[1]/ab.vp[2],N)),
                  as.double(rep(0,N)), #
                  PACKAGE="dcemriS4")

  cat("\b\b. Burnin phase")
  
  # burnin
    	    singlerun <- .C("dce_space",
                            as.integer(c(burnin+1-thin, burnin-thin, 0)),
                            as.double(conc),
                            as.integer(img.mask),
                            as.integer(dim(conc)),				
                            as.double(ab.gamma),
                            as.double(ab.gamma3),
                            as.double(ab.theta),
                            as.double(ab.theta3),
                            as.double(ab.vp),
                            as.double(ab.tauepsilon),
                            as.double(c(aif.model, aif.parameter)),
                            as.integer(c(vp.do,spatial,slice,gemupdate,uptauep,0,0)),
							as.double(c(0,time)), # time shifted, s.t. it fits with conc
                            as.double(singlerun[[14]]),
                            as.double(singlerun[[15]]),
                            as.double(singlerun[[16]]),
                            as.double(singlerun[[17]]),
                            as.double(singlerun[[18]]),
                            as.double(singlerun[[19]]),
                            as.double(singlerun[[20]]),
                            as.double(singlerun[[21]]),
                            as.double(singlerun[[22]]),
                            as.double(singlerun[[23]]),
                            as.double(singlerun[[24]]),
                            as.integer(rep(0,N)), #25
                            as.integer(rep(0,N)),
                            as.integer(rep(0,N)),
                            as.integer(rep(0,N)),
                            as.double(singlerun[[29]]),
                            as.double(singlerun[[30]]),
                            as.double(singlerun[[31]]),
                            as.double(singlerun[[32]]),
                            as.double(singlerun[[33]]),
                            as.double(singlerun[[34]]),
                            as.double(singlerun[[35]]),
                            as.double(singlerun[[36]]),
                            as.double(rep(0,N)), #
                            PACKAGE="dcemriS4")
  
  cat(" done. MCMC iteration 0")

  samplesize <- floor(nriters/thin)
  
  ktrans<-array(NA,c(II,JJ,KK,samplesize))
  kep<-array(NA,c(II,JJ,KK,samplesize))
  vp<-array(NA,c(II,JJ,KK,samplesize))
  sigma2<-array(NA,c(II,JJ,KK,samplesize))
  deviance<-array(NA,c(II,JJ,KK,samplesize))
  
  acc_ktrans<-array(NA,c(II,JJ,KK,samplesize))
  acc_kep<-array(NA,c(II,JJ,KK,samplesize))
  acc <- array(NA,c(II,JJ,KK,samplesize))
   
  iters<-0
  
  for (i in 1:samplesize)
    {
      #print(singlerun[[18]][1])
      singlerun <- .C("dce_space",
                      as.integer(c(thin, 0, 0)),
                      as.double(conc),
                      as.integer(img.mask),
                      as.integer(dim(conc)),
                      as.double(ab.gamma),
                      as.double(ab.gamma3),
                      as.double(ab.theta),
                      as.double(ab.theta3),
                      as.double(ab.vp),
                      as.double(ab.tauepsilon),
                      as.double(c(aif.model, aif.parameter)),
                      as.integer(c(vp.do,spatial,slice,gemupdate,uptauep,0,0)),
					  as.double(c(0,time)), # time shifted, s.t. it fits with conc
                      as.double(singlerun[[14]]),
                      as.double(singlerun[[15]]),
                      as.double(singlerun[[16]]),
                      as.double(singlerun[[17]]),
                      as.double(singlerun[[18]]),
                      as.double(singlerun[[19]]),
                      as.double(singlerun[[20]]),
                      as.double(singlerun[[21]]),
                      as.double(singlerun[[22]]),
                      as.double(singlerun[[23]]),
                      as.double(singlerun[[24]]),
                      as.integer(rep(0,N)), #25 acc_gamma
                      as.integer(rep(0,N)), # acc_theta
                      as.integer(rep(0,N)), # acc_eta
                      as.integer(rep(0,N)), # acc
                      as.double(singlerun[[29]]),
                      as.double(singlerun[[30]]),
                      as.double(singlerun[[31]]),
                      as.double(singlerun[[32]]),
                      as.double(singlerun[[33]]),
                      as.double(singlerun[[34]]),
                      as.double(singlerun[[35]]),
                      as.double(singlerun[[36]]),
                      as.double(rep(0,N)), #
                      PACKAGE="dcemriS4")
      
      ktrans[,,,i]<-array(singlerun[[17]],c(II,JJ,KK))
      kep[,,,i]<-array(singlerun[[18]],c(II,JJ,KK))
      vp[,,,i]<-array(singlerun[[19]],c(II,JJ,KK))
      sigma2[,,,i]<-1/array(singlerun[[20]],c(II,JJ,KK))
      deviance[,,,i]<-array(singlerun[[37]],c(II,JJ,KK))
	  	  
	  acc_ktrans[,,,i]<-array(singlerun[[25]],c(II,JJ,KK))
	  acc_kep[,,,i]<-array(singlerun[[26]],c(II,JJ,KK))	  
	  acc[,,,i]<-array(singlerun[[28]],c(II,JJ,KK))
	  
      if (iters==0)
        {
          cat("\b")
        }
      else
        {
          for (j in 1:floor(log(10*iters)/log(10)))cat("\b")
        }
      iters<-iters+thin
      cat (iters)
    }

  cat("\b s done. Preparing Results.\n")

  ktrans.med<-array(NA,c(I,J,K))
  kt<-apply(ktrans,1:3,median)
  kt[img.mask==0]<-NA
  ktrans.med[x,y,]<-kt

  kep.med<-array(NA,c(I,J,K))
  kp<-apply(kep,1:3,median)
  kp[img.mask==0]<-NA
  kep.med[x,y,]<-kp

  ve.med<-array(NA,c(I,J,K))
  temp<-apply(ktrans/kep,1:3,median)
  temp[img.mask==0]<-NA
  ve.med[x,y,]<-temp

  vp.med<-array(NA,c(I,J,K))
  v<-apply(vp,1:3,median)
  v[img.mask==0]<-NA
  vp.med[x,y,]<-v

  sigma2.med<-array(NA,c(I,J,K))
  s2<-apply(sigma2,1:3,median)
  s2[img.mask==0]<-NA
  sigma2.med[x,y,]<-s2

  returnable <- list(ktrans=ktrans.med,
                     kep=kep.med,
                     ve=ve.med,
                     vp=vp.med,
                     sigma2=sigma2.med
                     )

  if (dic) {
    if (verbose) {
      cat("  Computing DIC...", fill=TRUE)
    }
    fitted <- array(NA, c(II,JJ,KK,T))
    for (i in 1:II) {
      for (j in 1:JJ) {
        for (k in 1:KK) {
          if (img.mask[i,j,k]==1) {
            par <- NULL
			# Caution: order of assignment is important!!!
			if (vp.do) {
				par["vp"] <- v[i,j,k]
			}
            par["ktrans"]=kt[i,j,k]
            par["kep"]=kp[i,j,k]

            fitted[i,j,k,] <- kineticModel(time, par, model=model, aif=aif)
          }
        }
      }
    }
	
	returnable[["fitted"]] <- fitted
	
    conc <- array(conc, c(II,JJ,KK,length(time)))
    fitted <- fitted - conc
    fitted <- fitted * fitted
    fitted <- apply(fitted, 1:3, sum)
    deviance.med <- length(time) * log(s2) + fitted / s2
    med.deviance <- apply(deviance,1:3,median,na.rm=TRUE)
	    
    pD <- med.deviance - deviance.med
    DIC <- med.deviance + pD
    returnable[["DIC"]] <- sum(DIC,na.rm=TRUE)
    returnable[["pD"]] <- sum(pD,na.rm=TRUE)
    returnable[["DIC.map"]] <- DIC
    returnable[["pD.map"]] <- pD 
    returnable[["deviance.med"]] <- deviance.med
    returnable[["med.deviance"]] <- med.deviance
    if (samples) {
      returnable[["deviance.sample"]] <- deviance
    }	
  }

  if (samples)
    {
      returnable[["ktrans.sample"]] <- ktrans
      returnable[["kep.sample"]] <- kep
      returnable[["ve.sample"]] <- ktrans/kep
      
	  returnable[["kep.acc"]] <- acc_kep
	  returnable[["ktrans.acc"]] <- acc_ktrans
	  returnable[["acc"]] <- acc
	  	  
      if (vp.do==3)
        {
          returnable[["vp.sample"]] <- vp
        }

      returnable[["sigma2.sample"]] <- sigma2
    }	
  return(returnable)
}
  
  
