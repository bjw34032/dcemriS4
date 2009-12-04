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
## $Id$
##

#############################################################################
## setClass("nifti")
#############################################################################

setClass("nifti",
         representation("sizeof_hdr"="numeric",
                        "data_type"="character",
                        "db_name"="character",
                        "extents"="numeric",
                        "session_error"="numeric",
                        "regular"="character",
                        "dim_info"="numeric",
                        "dim_"="vector",
                        "intent_p1"="numeric",
                        "intent_p2"="numeric",
                        "intent_p3"="numeric",
                        "intent_code"="numeric",
                        "datatype"="numeric",
                        "bitpix"="numeric",
                        "slice_start"="numeric",
                        "pixdim"="vector",
                        "vox_offset"="numeric",
                        "scl_slope"="numeric",
                        "scl_inter"="numeric",
                        "slice_end"="numeric",
                        "slice_code"="numeric",
                        "xyzt_units"="numeric",
                        "cal_max"="numeric",
                        "cal_min"="numeric",
                        "slice_duration"="numeric",
                        "toffset"="numeric",
                        "glmax"="numeric",
                        "glmin"="numeric",
                        "descrip"="character",
                        "aux_file"="character",
                        "qform_code"="numeric",
                        "sform_code"="numeric",
                        "quatern_b"="numeric",
                        "quatern_c"="numeric",
                        "quatern_d"="numeric",
                        "qoffset_x"="numeric",
                        "qoffset_y"="numeric",
                        "qoffset_z"="numeric",
                        "srow_x"="vector",
                        "srow_y"="vector",
                        "srow_z"="vector",
                        "intent_name"="character",
                        "magic"="character",
                        "extender"="vector",
			"reoriented"="logical"),
         prototype("sizeof_hdr"=348,
                   "data_type"="",
                   "db_name"="",
                   "extents"=numeric(1),
                   "session_error"=numeric(1),
                   "regular"="",
                   "dim_info"=numeric(1),
                   "dim_"=numeric(8),
                   "intent_p1"=numeric(1),
                   "intent_p2"=numeric(1),
                   "intent_p3"=numeric(1),
                   "intent_code"=numeric(1),
                   "datatype"=2,
                   "bitpix"=8,
                   "slice_start"=numeric(1),
                   "pixdim"=numeric(8),
                   "vox_offset"=352,
                   "scl_slope"=numeric(1),
                   "scl_inter"=numeric(1),
                   "slice_end"=numeric(1),
                   "slice_code"=numeric(1),
                   "xyzt_units"=numeric(1),
                   "cal_max"=numeric(1),
                   "cal_min"=numeric(1),
                   "slice_duration"=numeric(1),
                   "toffset"=numeric(1),
                   "glmax"=numeric(1),
                   "glmin"=numeric(1),
                   "descrip"="",
                   "aux_file"="",
                   "qform_code"=numeric(1),
                   "sform_code"=numeric(1),
                   "quatern_b"=numeric(1),
                   "quatern_c"=numeric(1),
                   "quatern_d"=numeric(1),
                   "qoffset_x"=numeric(1),
                   "qoffset_y"=numeric(1),
                   "qoffset_z"=numeric(1),
                   "srow_x"=numeric(4),
                   "srow_y"=numeric(4),
                   "srow_z"=numeric(4),
                   "intent_name"="",
                   "magic"="n+1",
		   "extender"=numeric(4),
		   "reoriented"=FALSE),
         contains="array")

#############################################################################
## setClass("niftiExtension")
#############################################################################

setClass("niftiExtension",
         representation(extensions="list"),
         prototype(extensions=list()),
         contains="nifti")

#############################################################################
## setClass("niftiExtensionSection")
#############################################################################

setClass("niftiExtensionSection",
         representation(esize="numeric",
                        ecode="numeric",
                        edata="character"),
         prototype(esize=numeric(1),
                   ecode=numeric(1),
                   edata=""))

#############################################################################
## setMethod("show", "nifti")
#############################################################################

setMethod("show", "nifti", function(object) {
  cat("NIfTI-1 format", fill=TRUE)
  cat("  Type            :", class(object), fill=TRUE)
  cat("  Data Type       : ", object@"datatype", " (", convert.datatype(object@datatype), ")", sep="", fill=TRUE)
  cat("  Bits per Pixel  :", object@bitpix, fill=TRUE)
  cat("  Slice Code      : ", object@"slice_code", " (", convert.slice(object@"slice_code"), ")", sep="", fill=TRUE)
  cat("  Intent Code     : ", object@"intent_code", " (", convert.intent(object@"intent_code"), ")", sep="", fill=TRUE)
  cat("  Qform Code      : ", object@"qform_code", " (", convert.form(object@"qform_code"), ")", sep="", fill=TRUE)
  cat("  Sform Code      : ", object@"sform_code", " (", convert.form(object@"sform_code"), ")", sep="", fill=TRUE)
  cat("  Dimension       :", paste(object@"dim_"[2:(1+object@"dim_"[1])], collapse=" x "), fill=TRUE)
  cat("  Pixel Dimension :", paste(round(object@pixdim[2:(1+object@"dim_"[1])], 2), collapse=" x "), fill=TRUE)
  cat("  Voxel Units     :", convert.units(xyzt2space(object@"xyzt_units")), fill=TRUE)
  cat("  Time Units      :", convert.units(xyzt2time(object@"xyzt_units")), fill=TRUE)
})

#############################################################################
## setValidity("nifti")
#############################################################################

setValidity("nifti", function(object) {
  retval <- NULL
  indices <- 2:(1+object@"dim_"[1])
  ## sizeof_hdr must be 348
  if (object@"sizeof_hdr" != 348)
    retval <- c(retval, "sizeof_hdr != 348")
  ## datatype needed to specify type of image data
  if (!object@datatype %in% convert.datatype())
    retval <- c(retval, "datatype not recognized")
  ## bitpix should correspond correctly to datatype
  if (!object@bitpix == convert.bitpix()[[convert.datatype(object@datatype)]]) 
    retval <- c(retval, "bitpix does not match the datatype")
  ## dim should be non-zero for dim[1] dimensions
  if (!all(as.logical(object@"dim_"[indices])))
    retval <- c(retval, "dim[1]/dim mismatch")
  ## number of data dimensions should match dim[1]
  if (length(indices) != length(dim(object@.Data)))
    retval <- c(retval, "dim[1]/img mismatch")
  ## 
  if (object@"cal_min" != min(object@.Data) ||
      object@"cal_max" != max(object@.Data))
    retval <- c(retval, "range(img) != c(cal_min,cal_max)")
  ## pixdim[0] is required when qform_code != 0
  if (object@"qform_code" != 0 && object@pixdim[1] == 0)
    retval <- c(retval, "pixdim[1] is required")              
  ## pixdim[n] required when dim[n] is required
  if (!all(as.logical(object@"dim_"[indices]) &&
           as.logical(object@"pixdim"[indices])))
    retval <- c(retval, "dim/pixdim mismatch")
  ## data dimensions should match dim 
  if (!all(object@"dim_"[indices] == dim(object@.Data)))
    retval <- c(retval, "dim/img mismatch")
  ## vox_offset required for an "n+1" header
  if (match(object@"magic", "n+1") == 1 && object@"vox_offset" == 0)
    retval <- c(retval, "vox_offset required when magic=\"n+1\"")
  ## magic must be "ni1\0" or "n+1\0"
  if (!(match(object@"magic", "n+1") == 1 || match(object@"magic", "ni1") == 1))
    retval <- c(retval, "magic != \"n+1\" and magic != \"ni1\"")
  if (is.null(retval)) return(TRUE)
  else return(retval)
})

#############################################################################
## setValidity("niftiExtension")
#############################################################################

setValidity("niftiExtension", function(object) {
  ## Allegedly setValidity will always check for superclasses
  ## So we need only check that the list is empty or only contains niftiExtensionSections and check the validity of each of those
  retval <- NULL
  validSection <- getValidity(getClassDef("niftiExtensionSection"))
  lapply(object@"extensions",
         function(x) { 
           if (!is(x, "niftiExtensionSection")) {
             retval <<- c(retval, paste("@extensions list contains non-niftiExtensionSection element:", class(x)))
           } else {
             if (!(validSection(x) == TRUE)) {
               retval <<- c(retval, validSection(x))
             }
           }
         })
  if (is.null(retval))
    return(TRUE)
  else
    return(retval)
})

#############################################################################
## setValidity("niftiExtensionSection")
#############################################################################

setValidity("niftiExtensionSection", function(object) {
  retval <- NULL
  if (object@esize %% 16 != 0)
    retval <- c(retval, "esize is not a multiple of 16")
  if ((object@esize - 8) < nchar(object@edata, type="bytes"))
    retval <- c(retval, "esize is too small for the data contained within the section")
  if (is.null(retval))
    return(TRUE)
  else
    return(retval)
})

## setGeneric("img", function(object) { standardGeneric("img") })
## setMethod("img", "nifti", function(object) { object@.Data })
## setGeneric("img<-", function(x, value) { standardGeneric("img<-") })
## setReplaceMethod("img", signature(x="nifti", value="array"),
##           function(x, value) {
##             x@.Data <- value
##             x
##           })

## setGeneric("hdr", function(object) { standardGeneric("hdr") })
## setMethod("hdr", "nifti", # signature(object="nifti", name="ANY"),
##           function(object) { object@"descrip" })
## setGeneric("hdr<-", function(x, name, value) { standardGeneric("hdr<-") })
## setReplaceMethod("hdr", signature(x="nifti", name="character", value="ANY"),
##           function(x, name, value) {
##             x@name <- value
##             validObject(x)
##             x
##           })

#############################################################################
## nifti()
#############################################################################

nifti <- function(img=array(0, dim=rep(1,4)), dim, datatype=2, cal.min=NULL,
                  cal.max=NULL, ...) {
  ## Set dimensions
  if (missing(dim)) {
    if (is.array(img))
      dim <- base::dim(img)
    else
      dim <- c(1, length(img))
  }
  ld <- length(dim)
  if (ld < 3)
    stop(sprintf("length(dim) must be at least 3 and is %d.", ld))

  ## Create "dim" and "pixdim" slots
  x <- rep(1, 8)
  x[1] <- length(dim)
  y <- rep(0.0, 8)
  for (i in 2:length(x)) {
    x[i] <- ifelse(is.na(dim(img)[i-1]), 1, dim(img)[i-1])
    y[i] <- ifelse(is.na(dim(img)[i-1]), 1.0, 1.0)
  }
  ## min/max values for visualization
  if (is.null(cal.max))
    cal.max <- as.numeric(max(img, na.rm=TRUE))
  if (is.null(cal.min))
    cal.min <- as.numeric(min(img, na.rm=TRUE))
  ## Set datatype
  switch(as.character(datatype),
         "2" = bitpix <- 8,
         "4" = bitpix <- 16,
         "8" = bitpix <- 32,
         "16" = bitpix <- 32,
         "64" = bitpix <- 64,
         "512" = bitpix <- 16,
         stop(paste("Data type", datatype, "unsupported."))
         )
  ## Create the object
  niftiClass <- "nifti"
  if (getOption("NIfTI.audit.trail")) {
    niftiClass <- "niftiAuditTrail"
  }
  obj <- new(niftiClass, .Data=array(img, dim=dim), "dim_"=x, "pixdim"=y,
             "cal_max"=cal.max, "cal_min"=cal.min, "datatype"=datatype,
             "bitpix"=bitpix, ...)
  if (getOption("NIfTI.audit.trail")) {
    audit.trail(obj) <- niftiAuditTrailCreated(call=match.call())
  }
  validObject(obj)
  return(obj)
}

#############################################################################
## is.nifti()
#############################################################################

is.nifti <- function(x) {
  if (!is(x, "nifti"))
    return(FALSE)
  else
    return (TRUE)
}

#############################################################################
## descrip() accessor function to @descrip
#############################################################################

setGeneric("descrip", function(object) { standardGeneric("descrip") })
setMethod("descrip", "nifti", function(object) { object@descrip })
setGeneric("descrip<-", function(x, value) { standardGeneric("descrip<-") })
setReplaceMethod("descrip", "nifti",
                 function(x, value) { 
		   x@descrip <- value 
		   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", match.call(),
                                          paste("descrip <-", value))
		   x 
		 })

#############################################################################
## aux.file() accessor function to @"aux_file"
#############################################################################

setGeneric("aux.file", function(object) { standardGeneric("aux.file") })
setMethod("aux.file", "nifti", function(object) { object@"aux_file" })
setGeneric("aux.file<-", function(x, value) { standardGeneric("aux.file<-") })
setReplaceMethod("aux.file", "nifti",
                 function(x, value) {
		   x@"aux_file" <- value
		   audit.trail(x) <- niftiAuditTrailEvent(x, "modification", match.call(), paste("aux.file <-", value))
		   x 
		 })

#############################################################################
## audit.trail() accessor function to @"trail"
#############################################################################
## These functions will work even if the audit trail functionality is not
## activated. They should help reduce the difference in code paths.
#############################################################################

setGeneric("audit.trail", function(object) { standardGeneric("audit.trail") })
setMethod("audit.trail", "nifti",
          function(object) { 
            if (getOption("NIfTI.audit.trail") &&
                is(object, "niftiAuditTrail")) {
	      return(object@"trail")
            } else {
              NULL
            }
          })

setGeneric("audit.trail<-",
           function(x, value) { standardGeneric("audit.trail<-") })
setReplaceMethod("audit.trail", signature(x="nifti", value="character"),
                 function(x, value) {
                   if (getOption("NIfTI.audit.trail")) {
                     if (!is(x, "niftiAuditTrail")) {
                       x <- as(x, "niftiAuditTrail")
                     }
                     x@"trail" <- value
                   } 
                   x
                 })
setReplaceMethod("audit.trail", signature(x="nifti", value="XMLAbstractNode"),
                 function(x, value) {
                   if (getOption("NIfTI.audit.trail")) {
                     if (!is(x, "niftiAuditTrail")) {
                       x <- as(x, "niftiAuditTrail")
                     }
                     x@"trail" <- saveXML(value)
                   } 
                   x
                 })

setReplaceMethod("[", signature(x="nifti", i="missing", j="missing",
                                value="array"),
                 function(x, value) {
                   x <- as.nifti(value, x)
                   validObject(x)
                   x
                 })

setReplaceMethod("[", signature(x="nifti", i="ANY", j="missing", value="ANY"), 
                 function(x, i, value) {
                   # For some reason this line is slow; I don't understand it
                   x@.Data[i] <- value
                   #if (any(is.na(value))) {
                   if (any(value < x@"cal_min", na.rm=TRUE))
                     x@"cal_min" <- min(value, na.rm=TRUE)
                   if (any(value > x@"cal_max", na.rm=TRUE))
                     x@"cal_max" <- max(value, na.rm=TRUE)
                   #} else {
                   #  xr <- range(x, na.rm=TRUE)
                   #  x@"cal_min" <- xr[1]
                   #  x@"cal_max" <- xr[2]
                   #}
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", 
                                          call=sys.call(-3),
                                          comment=paste("Non-numeric replace ["))
                   x
                 })
setReplaceMethod("[", signature(x="nifti", i="numeric", j="missing", value="ANY"), 
                 function(x, i, value) {
                   # For some reason this line is slow; I don't understand it
                   x@.Data[i] <- value
                   #if (any(is.na(value))) {
                   if (any(value < x@"cal_min", na.rm=TRUE))
                     x@"cal_min" <- min(value, na.rm=TRUE)
                   if (any(value > x@"cal_max", na.rm=TRUE))
                     x@"cal_max" <- max(value, na.rm=TRUE)
                   #} else {
                   #  xr <- range(x, na.rm=TRUE)
                   #  x@"cal_min" <- xr[1]
                   #  x@"cal_max" <- xr[2]
                   #}
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", 
                                          call=sys.call(-3))
                   x
                 })

setReplaceMethod("[", signature(x="nifti", i="ANY", j="ANY", value="ANY"),
                 function(x, i, j, ..., value) {
                   # For some reason this line is slow; I don't understand it
                   x@.Data[i,j,...] <- value
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", call=sys.call(-3),
                                          comment=paste("Non-numeric replace ["))
					  
                   x
                 })
setReplaceMethod("[", signature(x="nifti", i="numeric", j="numeric", value="ANY"),
                 function(x, i, j, ..., value) {
                   # For some reason this line is slow; I don't understand it
                   x@.Data[i,j,...] <- value
                   #if (all(is.na(value))) {
                   #  xr <- range(x, na.rm=TRUE)
                   #  x@"cal_min" <- xr[1]
                   #  x@"cal_max" <- xr[2]
                   #} else {
                   if (any(value < x@"cal_min", na.rm=TRUE))
                     x@"cal_min" <- min(value, na.rm=TRUE)
                   if (any(value > x@"cal_max", na.rm=TRUE))
                     x@"cal_max" <- max(value, na.rm=TRUE)
                   #}
                   audit.trail(x) <-
                     niftiAuditTrailEvent(x, "modification", 
                                          comment=paste("[", paste(i, j, ..., sep=", "), "] <- ", value, sep=""))
                   x
                 })



#############################################################################
## quaternion2rotation()
#############################################################################

quaternion2rotation <- function(b, c, d) {
  a <- sqrt(1 - (b*b+c*c+d*d))
  R <- matrix(c(a*a+b*b-c*c-d*d, 2*b*c+2*a*d, 2*b*d-2*a*c,  # column 1
                2*b*c-2*a*d, a*a+c*c-b*b-d*d, 2*c*d+2*a*b,  # column 2
                2*b*d+2*a*c, 2*c*d-2*a*b, a*a+d*d-c*c-b*b), # column 3
              3, 3)
  return(R)
}

