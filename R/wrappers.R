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
## $Id: $
##

#############################################################################
## Convert output from dcemri functions into S4 objects
#############################################################################

dcemri.lm.s4 <- function(conc, time, mask, ...) {
  as.nifti(dcemri::dcemri.lm(conc, time, mask, ...), conc)
}

dcemri.bayes.s4 <- function(conc, time, img.mask, ...) {
  as.nifti(dcemri::dcemri.bayes(conc, time, img.mask, ...), conc)
}

dcemri.spline.s4 <- function(conc, time, img.mask, ...) {
  as.nifti(dcemri::dcemri.spline(conc, time, img.mask, ...), conc)
}

#############################################################################
## setGeneric("dcemri.lm")
#############################################################################

setGeneric("dcemri.lm",
           function(conc, ...) standardGeneric("dcemri.lm"))
setMethod("dcemri.lm", signature(conc="nifti"), dcemri.lm.s4)
setMethod("dcemri.lm", signature(conc="array"),
          function(conc, ...) { 
            dcemri.lm.s4(as(conc, "nifti"), ...)
          })
setMethod("dcemri.lm", signature(conc="anlz"),
          function(conc, ...) { 
            dcemri.lm.s4(as(conc, "nifti"), ...)
          })

#############################################################################
## setGeneric("dcemri.bayes")
#############################################################################

setGeneric("dcemri.bayes",
           function(conc, ...) standardGeneric("dcemri.bayes"))
setMethod("dcemri.bayes", signature(conc="nifti"), dcemri.bayes.s4)
setMethod("dcemri.bayes", signature(conc="array"), function(conc,...) { 
      dcemri.bayes.s4(as("nifti", conc), ...)
    })
setMethod("dcemri.bayes", signature(conc="anlz"), function(conc,...) { 
      dcemri.bayes.s4(as("nifti", conc), ...)
    })

#############################################################################
## setGeneric("dcemri.spline")
#############################################################################

setGeneric("dcemri.spline",
           function(conc, ...) standardGeneric("dcemri.spline"))
setMethod("dcemri.spline", signature(conc="nifti"), dcemri.spline.s4)
setMethod("dcemri.spline", signature(conc="array"), function(conc,...) { 
      dcemri.spline.s4(as("nifti", conc), ...)
    })
setMethod("dcemri.spline", signature(conc="anlz"), function(conc,...) { 
      dcemri.spline.s4(as("nifti", conc), ...)
    })
