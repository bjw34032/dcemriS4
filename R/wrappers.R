dcemri.lm.s4 <- function(nim, time, mask, ...) {
  as.nifti(dcemri.lm(nim, time, mask, ...), nim)
}

dcemri.bayes.s4 <- function(nim, time, img.mask, ...) {
  as.nifti(dcemri.bayes(nim, time, img.mask, ...), nim)
}

dcemri.spline.s4 <- function(nim, time, img.mask, ...) {
  as.nifti(dcemri.spline(nim, time, img.mask, ...), nim)
}
