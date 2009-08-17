dcemri.lm.s4 <- function(nifticonc, time, mask, model="extended", aif=NULL,
    nprint=0, user=NULL, verbose=FALSE, ...) {
  ## input:
  ##	    nifticonc: a niftis4 object of Gd concentrations
  ##	    time: timepoints of aquisition
  ##	    img.mask: array of voxels to fit
  

  lm.fit <- dcemri.lm(nifticonc, time, mask, model=model, aif=aif, nprint=nprint, user=user, verbose=verbose, ...)
  
  time <- lm.fit$time
  lm.fit$time <- NULL

  s4.fit <- lapply(lm.fit, function(x) { as.nifti(x, nifticonc) })
  lapply(names(s4.fit), function(x) { 
	niftiobject <- s4.fit[[x]]
	niftiobject@"intent_code" <- nifit.intent.code[["Estimate"]]
	niftiobject@"intent_name" <- substr(x,1,15)
      }) 

  s4.fit[["time"]] <- time
  return(s4.fit)
}
