dcemri.lm.s4 <- function(conc, time, mask, ...) {
  as.nifti(dcemri::dcemri.lm(conc, time, mask, ...), conc)
}

dcemri.bayes.s4 <- function(conc, time, img.mask, ...) {
  as.nifti(dcemri::dcemri.bayes(conc, time, img.mask, ...), conc)
}

dcemri.spline.s4 <- function(conc, time, img.mask, ...) {
  as.nifti(dcemri::dcemri.spline(conc, time, img.mask, ...), conc)
}

setGeneric("dcemri.lm")
setMethod("dcemri.lm", signature(conc="nifti"), dcemri.lm.s4)
setMethod("dcemri.lm", signature(conc="array"), function(conc,...) { 
      dcemri.lm.s4(as("nifti", conc), ...)
    })
setMethod("dcemri.lm", signature(conc="anlz"), function(conc,...) { 
      dcemri.lm.s4(as("nifti", conc), ...)
    })

setGeneric("dcemri.bayes")
setMethod("dcemri.bayes", signature(conc="nifti"), dcemri.bayes.s4)
setMethod("dcemri.bayes", signature(conc="array"), function(conc,...) { 
      dcemri.bayes.s4(as("nifti", conc), ...)
    })
setMethod("dcemri.bayes", signature(conc="anlz"), function(conc,...) { 
      dcemri.bayes.s4(as("nifti", conc), ...)
    })

setGeneric("dcemri.spline")
setMethod("dcemri.spline", signature(conc="nifti"), dcemri.spline.s4)
setMethod("dcemri.spline", signature(conc="array"), function(conc,...) { 
      dcemri.spline.s4(as("nifti", conc), ...)
    })
setMethod("dcemri.spline", signature(conc="anlz"), function(conc,...) { 
      dcemri.spline.s4(as("nifti", conc), ...)
    })


