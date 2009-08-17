dcemri.lm.s4 <- function(nifti, time, mask, ...) {
  ## input:
  ##	    nifticonc: a niftis4 object of Gd concentrations
  ##	    time: timepoints of aquisition
  ##	    img.mask: array of voxels to fit
  

  lm.fit <- dcemri.lm(nifti, time, mask, ...)
  fit.nifti(lm.fit, nifti)
}

dcemri.bayes.s4 <- function(nifti, time, img.mask, ...) {
  bayes.fit <- dcemri.bayes(nifti, time, img.mask, ...) 
  fit.nifti(bayes.fit, nifti)
}

fit.nifti <- function(fit, nifti) {
  time <- fit$time
  fit$time <- NULL

  s4.fit <- lapply(fit, function(x) { as.nifti(x, nifti) })
  
  lapply(names(s4.fit), function(x) {
	niftiobject <- s4.fit[[x]]
	niftiobject@"intent_code" <- nifit.intent.code[["Estimate"]]
	niftiobject@"intent_name" <- substr(x,1,15)
      })

  s4.fit[["time"]] <- time
  return(s4.fit)
}
