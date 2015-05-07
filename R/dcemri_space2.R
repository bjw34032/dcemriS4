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
## setGeneric("dcemri.space2")
#############################################################################
setGeneric("dcemri.space2",
           function(conc, ...) standardGeneric("dcemri.space2"))
setMethod("dcemri.space2", signature(conc="array"), 
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
				  ab.gamma2comp=c(0.0001,0.0001),
				  ab.theta2comp=c(0.0001,0.0001),
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
          .dcemriWrapper("dcemri.space2", conc, time, img.mask, starting,
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
				  ab.gamma2comp,
				  ab.theta2comp,
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

.dcemri.space2 <- function(conc, time, img.mask, starting=0,
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
				ab.gamma2comp=c(0.0001,0.0001),
				ab.theta2comp=c(0.0001,0.0001),
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
	  ## dcemri.space2 - a function for fitting 2-compartment PK models 
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
	  ##					 1: estimate precision of error (estimate tau_epsilon)
	  ##	 				 2: 1+global smoothing (estimate tau_epsilon and global tau_theta / tau_gamma)
	  ##					 3: 1+adaptive smoothing (estimate tau_epsilon and local tau_theta_ij / tau_gamma_ij)
	  ##			tuning: tune after n iterations
	  ##			tunepct: percentage of voxels which can fail tuning
	  ##			retunecycles: number of times to do tuning
	  ##		  
	  ## output: list with ktrans, kep, ve, std.error of ktrans and kep
	  ##         (ktranserror and keperror), samples if samples=TRUE
	  ##  
	
    if (spatial==2)spatial=1
    
	  # 2Comp model
#	  if(model=="extended"){model="extended2"}
#	  if(model=="weinmann"){model="weinmann2"}
	  	  
	  switch(model,
	         weinmann = ,
			 weinmann2 = ,
			 extended2 = {
				 aif <- ifelse(is.null(aif), "tofts.kermode", aif)
				 aif.names <- c("tofts.kermode","fritz.hansen","empirical")
				 if (! aif %in% aif.names) {
					 stop(sprintf("Only aif=\"%s\" or aif=\"%s\" or aif=\"%s\" are acceptable AIFs for model=\"weinmann\" or model=\"extended\"", aif.names[1], aif.names[2], aif.names[3]), call.=FALSE)
				 } 
			 },
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
			 weinmann2 = {
				 aif.model <- 0
				 vp.do <- 0
			 },
			 extended2 = {
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
	  } 
	  else {
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
			kep_1 <- 0.2  
			kep_2 <- 5
			ve_1 <- 0.3
			ve_2 <- 0.3
			ktrans_1 <- kep_1 * ve_1
			ktrans_2 <- kep_2 * ve_2	
			
			ktrans1_start <- rep(ktrans_1,N)
			ktrans2_start <- rep(ktrans_2,N)
			kep1_start <- rep(kep_1,N)
			kep2_start <- rep(kep_2,N)
			
		}
		else 
		{			
			ktrans1_start <- starting$ktrans
			ktrans2_start <- starting$ktrans2
			kep1_start <- starting$kep
			kep2_start <- starting$kep2
					
			ktrans1_start[is.na(ktrans1_start)] <- 0
			ktrans2_start[is.na(ktrans2_start)] <- 0
			kep1_start[is.na(kep1_start)] <- 0
			kep2_start[is.na(kep2_start)] <- 0
			
		}
		
	  	# tuning step
	  	singlerun <- .C("dce_space_2comp",
			  as.integer(c(tuning+2, tuning+1, tuning)),
			  as.double(conc),
			  as.integer(c(img.mask,0)),
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
			  as.double(rep(0,N+1)),	# t0						  
			  as.double(c(ktrans1_start,0)), # ktrans #15
			  as.double(c(kep1_start,0)), # kep
			  as.double(c(ktrans2_start,0)), # ktrans2 
			  as.double(c(kep2_start,0)), # kep2			  
			  as.double(if(vp.do==3){rep(0.1,N+1)}else{rep(0,N+1)}), # vp
			  as.double(rep(ab.tauepsilon[1]/ab.tauepsilon[2],N+1)), # tau_epsilon 20	  			  
			  as.double(rep(ab.gamma[1]/ab.gamma[2],N+1)), # tau_gamma
			  as.double(rep(ab.theta[1]/ab.theta[2],N+1)), # tau_theta
			  as.double(rep(ab.gamma2comp[1]/ab.gamma2comp[2],N+1)), # tau_gamma_2comp
			  as.double(rep(ab.theta2comp[1]/ab.theta2comp[2],N+1)), # tau_theta_2comp
			  as.double(rep(ab.vp[1]/ab.vp[2],N+1)), # tau_vp 25			  
			  as.double(rep(ab.gamma3[1]/ab.gamma3[2],N+1)), #tau_gamma3
			  as.double(rep(ab.theta3[1]/ab.theta3[2],N+1)), #tau_theta3
			  as.double(rep(ab.gamma3[1]/ab.gamma3[2],N+1)), #tau_gamma3_2comp
			  as.double(rep(ab.theta3[1]/ab.theta3[2],N+1)), #tau_theta3_2comp
			  as.double(rep(ab.vp[1]/ab.vp[2],N+1)),  # tau_vp3 #30			  
			  as.double(rep(0.5,N+1)), #sigmatheta
			  as.double(rep(0.5,N+1)), #sigmagamma
			  as.double(rep(0.5,N+1)), #sigmatheta_2comp
			  as.double(rep(0.5,N+1)), #sigmagamma_2comp
			  as.double(rep(1,N+1)), #sigmavp #35
			  as.integer(rep(0,N+1)), # acc_gamma
			  as.integer(rep(0,N+1)), # acc_theta
			  as.integer(rep(0,N+1)), # acc_gamma_2comp
			  as.integer(rep(0,N+1)), # acc_theta_2comp
			  as.integer(rep(0,N+1)), # acc_vp #40
			  as.integer(rep(0,N+1)), # acc
			  as.double(rep(0,N+1)), # devi
			  as.double(rep(0,(N+1)*(T+1))), # fit
			  PACKAGE="dcemriS4")
			  			  	  
	 	 cat("\b\b. Burnin phase 2Comp")
	  
	  	# burnin
	  	singlerun <- .C("dce_space_2comp",
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
			  as.double(singlerun[[15]]), #ktrans
			  as.double(singlerun[[16]]), #kep
			  as.double(singlerun[[17]]), #ktrans2
			  as.double(singlerun[[18]]), #kep2
			  as.double(singlerun[[19]]), #vp
			  as.double(singlerun[[20]]), #tau_epsilon
			  as.double(singlerun[[21]]),
			  as.double(singlerun[[22]]),
			  as.double(singlerun[[23]]),
			  as.double(singlerun[[24]]),
			  as.double(singlerun[[25]]), #tau_vp
			  as.double(singlerun[[26]]), 
			  as.double(singlerun[[27]]),
			  as.double(singlerun[[28]]),
			  as.double(singlerun[[29]]),
			  as.double(singlerun[[30]]), # tau_vp3			  
			  as.double(singlerun[[31]]), #sigmatheta
			  as.double(singlerun[[32]]), #sigmagamma
			  as.double(singlerun[[33]]), #sigmatheta_2comp
			  as.double(singlerun[[34]]), #sigmagamma_2comp			  
			  as.double(singlerun[[35]]), #sigmavp			  
			  as.integer(rep(0,N+1)), # acc_
			  as.integer(rep(0,N+1)),
			  as.integer(rep(0,N+1)),
			  as.integer(rep(0,N+1)),
			  as.integer(rep(0,N+1)),
			  as.integer(rep(0,N+1)),
			  as.double(rep(0,N+1)), # devi
			  as.double(rep(0,(N+1)*(T+1))), # fit
			  PACKAGE="dcemriS4")
	  	  
	  cat(" done. 2comp MCMC iteration 0")
		  
	  samplesize <- floor(nriters/thin)
	  
	  #print(samplesize)
	  #print(c(II,JJ,KK,samplesize))
	  
	  #print("dim(singlerun[[14]])")
	  #dim(singlerun[[14]])
	  #print("dim(singlerun[[15]])")
	  #dim(singlerun[[15]])
	  
	  ktrans<-array(NA,c(II,JJ,KK,samplesize))
	  
	  #print("ktrans initialized")
	  
	  kep<-array(NA,c(II,JJ,KK,samplesize))
	  ktrans2<-array(NA,c(II,JJ,KK,samplesize))
	  kep2<-array(NA,c(II,JJ,KK,samplesize))
	  vp<-array(NA,c(II,JJ,KK,samplesize))
	  sigma2<-array(NA,c(II,JJ,KK,samplesize))
	  deviance<-array(NA,c(II,JJ,KK,samplesize))
	  
	  acc_ktrans<-array(NA,c(II,JJ,KK,samplesize))
	  acc_kep<-array(NA,c(II,JJ,KK,samplesize))
	  acc_ktrans2<-array(NA,c(II,JJ,KK,samplesize))
	  acc_kep2<-array(NA,c(II,JJ,KK,samplesize))	  
	  acc <- array(NA,c(II,JJ,KK,samplesize))
	  
	  #print("parameter initialization finished")
	  	  
	  iters<-0
	  
	  for (i in 1:samplesize)
	  {
		  #cat("Before next call of dce_space_2comp...")	
		  singlerun <- .C("dce_space_2comp",
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
				  as.double(c(0,time)),
				  as.double(singlerun[[14]]),
				  as.double(singlerun[[15]]), #ktrans
				  as.double(singlerun[[16]]), #kep
				  as.double(singlerun[[17]]), #ktrans2
				  as.double(singlerun[[18]]), #kep2
				  as.double(singlerun[[19]]), #vp
				  as.double(singlerun[[20]]), #tau_epsilon
				  as.double(singlerun[[21]]),
				  as.double(singlerun[[22]]),
				  as.double(singlerun[[23]]),
				  as.double(singlerun[[24]]),
				  as.double(singlerun[[25]]), #tau_vp
				  as.double(singlerun[[26]]), 
				  as.double(singlerun[[27]]),
				  as.double(singlerun[[28]]),
				  as.double(singlerun[[29]]),
				  as.double(singlerun[[30]]), # tau_vp3			  
				  as.double(singlerun[[31]]), #sigmatheta
				  as.double(singlerun[[32]]), #sigmagamma
				  as.double(singlerun[[33]]), #sigmatheta_2comp
				  as.double(singlerun[[34]]), #sigmagamma_2comp			  
				  as.double(singlerun[[35]]), #sigmavp			  
				  as.integer(rep(0,N+1)), # acc_gamma
				  as.integer(rep(0,N+1)), # acc_theta
				  as.integer(rep(0,N+1)), # acc_gamma_2comp
				  as.integer(rep(0,N+1)), # acc_theta_2comp
				  as.integer(rep(0,N+1)), # acc_vp
				  as.integer(rep(0,N+1)), # acc
				  as.double(rep(0,N+1)), # devi
				  as.double(rep(0,(N+1)*(T+1))), # fit
				  PACKAGE="dcemriS4")				
		
			#cat("... after call of dce_space_2comp")	
				  
		   ktrans[,,,i]<-array(singlerun[[15]],c(II,JJ,KK))
		   kep[,,,i]<-array(singlerun[[16]],c(II,JJ,KK))
		   ktrans2[,,,i]<-array(singlerun[[17]],c(II,JJ,KK))
		   kep2[,,,i]<-array(singlerun[[18]],c(II,JJ,KK))
		   vp[,,,i]<-array(singlerun[[19]],c(II,JJ,KK))
		   sigma2[,,,i]<-1/array(singlerun[[20]],c(II,JJ,KK))
		   deviance[,,,i]<-array(singlerun[[42]],c(II,JJ,KK))
		   
		   acc_ktrans[,,,i]<-array(singlerun[[36]],c(II,JJ,KK))
		   acc_kep[,,,i]<-array(singlerun[[37]],c(II,JJ,KK))
		   acc_ktrans2[,,,i]<-array(singlerun[[38]],c(II,JJ,KK))
		   acc_kep2[,,,i]<-array(singlerun[[39]],c(II,JJ,KK))
		   
		   acc[,,,i]<-array(singlerun[[41]],c(II,JJ,KK))
		 		 
		  		  
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
     } # end for (i in 1:samplesize)
	
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
	 				 
	 ktrans2.med<-array(NA,c(I,J,K))
	 kt2<-apply(ktrans2,1:3,median)
	 kt2[img.mask==0]<-NA
	 ktrans2.med[x,y,]<-kt2
	 
	 kep2.med<-array(NA,c(I,J,K))
	 kp2<-apply(kep2,1:3,median)
	 kp2[img.mask==0]<-NA
	 kep2.med[x,y,]<-kp2
	 
	 returnable[["ktrans2"]] <- ktrans2.med
	 returnable[["kep2"]] <- kep2.med
	 		 					 				 	
	 if (dic)
	 {
	    if (verbose) 
		{
	      cat("  Computing DIC...", fill=TRUE)
	    }
		
		fitted <- array(NA, c(I,J,K,T))
	    for (i in 1:I) {
	      for (j in 1:J) {
	        for (k in 1:K) {
	          if (img.mask[i,j,k]==1) 
			  {
	            par <- NULL
				# Caution: order of assignment is important!!!					
				if (vp.do) {
					par["vp"] <- vp.med[i,j,k]
				}
				
				par["ktrans"]=ktrans.med[i,j,k]
				par["kep"]=kep.med[i,j,k]
				par["ktrans2"]=ktrans2.med[i,j,k]
				par["kep2"]=kep2.med[i,j,k]
							
				fitted[i,j,k,] <- kineticModel(time, par, model=model, aif=aif)
				
	          }
	        }
	      }
	    }
		
		returnable[["fitted"]] <- fitted
		
		conc <- array(conc, c(I,J,K,length(time)))
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
	    if (samples) 
		{
	      returnable[["deviance.sample"]] <- deviance
	    }
		
	  }
	
	  if(samples)
	  {
	      returnable[["ktrans.sample"]] <- ktrans
	      returnable[["kep.sample"]] <- kep
	      returnable[["ve.sample"]] <- ktrans/kep	  
		  returnable[["kep.acc"]] <- acc_kep
		  returnable[["ktrans.acc"]] <- acc_ktrans
		  returnable[["acc"]] <- acc
		  
		  returnable[["ktrans2.sample"]] <- ktrans2
		  returnable[["kep2.sample"]] <- kep2
		  returnable[["ve2.sample"]] <- ktrans2/kep2		  
		  returnable[["kep2.acc"]] <- acc_kep2
		  returnable[["ktrans2.acc"]] <- acc_ktrans2		  		  
	      
	      if (vp.do==3)
	        {
	          returnable[["vp.sample"]] <- vp
	        }
	
	      returnable[["sigma2.sample"]] <- sigma2
	  }
	  	
	  return(returnable)
}
  
