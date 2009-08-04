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
## Time-stamp: <2009-07-22 15:15:11 (bjw34032)>
## $ Id: $
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
                        "extender"="vector"),
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
                   "extender"=numeric(4)),
         contains="array")

#############################################################################
## setClass("extension")
#############################################################################

setClass("niftiExtension",
         representation(esize="numeric",
                        ecode="numeric",
                        edata="character"),
         prototype(esize=numeric(1),
                   ecode=numeric(1),
                   edata=""),
         contains="nifti")

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
  if (!object@datatype %in% c(2^(1:11),768,1280,1536,1792))
    retval <- c(retval, "datatype not recognized")
  ## bitpix should correspond correctly to datatype
  
  ## dim should be non-zero for dim[1] dimensions
  if (!all(as.logical(object@"dim_"[indices])))
    retval <- c(retval, "dim[1]/dim mismatch")
  ## number of data dimensions should match dim[1]
  if (length(indices) != length(dim(object@.Data)))
    retval <- c(retval, "dim[1]/img mismatch")
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

nifti <- function(img=array(0, dim=rep(1,4)), dim, ...) {
  if (missing(dim)) {
    if (is.array(img))
      dim <- base::dim(img)
    else
      dim <- c(1, length(img))
  }
  ld <- length(dim)
  if (ld < 3)
    stop(sprintf("length(dim) must be at least 3 and is %d.", ld))
  
  x <- c(length(dim), dim[1], dim[2], dim[3],
         ifelse(is.na(dim[4]), 1, dim[4]), rep(1,3))
  y <- c(0.0, rep(1.0,length(dim)), rep(0.0,3))
  cal.max <- quantile(img, probs=0.95, na.rm=TRUE)
  cal.min <- quantile(img, probs=0.05, na.rm=TRUE)
  obj <- new("nifti", .Data=array(img, dim=dim), "dim_"=x, "pixdim"=y,
             "cal_max"=cal.max, "cal_min"=cal.min, ...)
  validObject(obj)
  return(obj)
}

is.nifti <- function(x) {
  if (!is(x, "nifti"))
    return(FALSE)
  else
    return (TRUE)
}

if (!isGeneric("descrip")) {
  if (is.function("descrip"))
    fun <- descrip
  else
    fun <- function(object) { standardGeneric("descrip") }
  setGeneric("descrip", fun)
}
setMethod("descrip", "nifti", function(object) { object@descrip })
setGeneric("descrip<-", function(x, value) { standardGeneric("descrip<-") })
setReplaceMethod("descrip", "nifti",
                 function(x, value) { x@descrip <- value ; x })

if (!isGeneric("aux.file")) {
  if (is.function("aux.file"))
    fun <- aux.file
  else
    fun <- function(object) { standardGeneric("aux.file") }
  setGeneric("aux.file", fun)
}
setMethod("aux.file", "nifti", function(object) { object@"aux_file" })
setGeneric("aux.file<-", function(x, value) { standardGeneric("aux.file<-") })
setReplaceMethod("aux.file", "nifti",
                 function(x, value) { x@"aux.file" <- value ; x })

convert.datatype <- function(dt) {
  ## defgroup NIFTI1_DATATYPE_ALIASES
  ## brief aliases for the nifti1 datatype codes
  switch(as.character(dt),
         "2" = "UINT8",
         "4" = "INT16",
         "8" = "INT32",
         "16" = "FLOAT32",
         "32" = "COMPLEX64",
         "64" = "FLOAT64",
         "128" = "RGB24",
         "256" = "INT8",
         "512" = "UINT16",
         "768" = "UINT32",
         "1024" = "INT64",
         "1280" = "UINT64",
         "1536" = "FLOAT128",
         "1792" = "COMPLEX128",
         "2048" = "COMPLEX256")
}

convert.intent <- function(ic) {
  ## defgroup NIFTI1_INTENT_CODES
  ##-------- These codes are for probability distributions -----------
  ## Most distributions have a number of parameters, below denoted
  ## by p1, p2, and p3, and stored in
  ##  - intent_p1, intent_p2, intent_p3 if dataset doesn't have 5th
  ##    dimension
  ##  - image data array                if dataset does have 5th
  ##    dimension
  ##
  ## Functions to compute with many of the distributions below can
  ## be found in the CDF library from U Texas.
  ##
  ## Formulas for and discussions of these distributions can be
  ## found in the following books:
  ##
  ##  [U] Univariate Discrete Distributions,
  ##      NL Johnson, S Kotz, AW Kemp.
  ##
  ##  [C1] Continuous Univariate Distributions, vol. 1,
  ##       NL Johnson, S Kotz, N Balakrishnan.
  ##
  ##  [C2] Continuous Univariate Distributions, vol. 2,
  ##       NL Johnson, S Kotz, N Balakrishnan.
  ##
  switch(as.character(ic),
         "0" = "None",
         "2" = "Correl",
         "3" = "Ttest",
         "4" = "Ftest",
         "5" = "Zscore",
         "6" = "Chisq",
         "7" = "Beta",
         "8" = "Binom",
         "9" = "Gamma",
         "10" = "Poisson",
         "11" = "Normal",
         "12" = "Ftest_Nonc",
         "13" = "Chisq_Nonc",
         "14" = "Logistic",
         "15" = "Laplace",
         "16" = "Uniform",
         "17" = "Ttest_Nonc",
         "18" = "Weibull",
         "19" = "Chi",
         "20" = "Invgauss",
         "21" = "Extval",
         "22" = "Pval",
         "23" = "Logpval",
         "24" = "Log10pval",
         "1001" = "Estimate",      # estimate of some parameter
         "1002" = "Label",         # index into some set of labels
         "1003" = "Neuroname",     # index into the NeuroNames labels set
         "1004" = "Genmatrix",     # M x N matrix at each voxel
         "1005" = "Symmatrix",     # N x N symmetric matrix at each voxel
         "1006" = "Dispvect",      # a displacement field
         "1007" = "Vector",        # a displacement vector
         "1008" = "Pointset",      # a spatial coordinate
         "1009" = "Triangle",      # triple of indexes 
         "1010" = "Quaternion",    # a quaternion
         "1011" = "Dimless")       # Dimensionless value - no params
}

convert.form <- function(fc) {
  ## defgroup NIFTI1_XFORM_CODES
  ## brief nifti1 xform codes to describe the "standard" 
  ## coordinate system
  switch(as.character(fc),
         "0" = "Unknown",       # Arbitrary coordinates (Method 1)
         "1" = "Scanner_Anat",  # Scanner-based anatomical coordinates
         "2" = "Aligned_Anat",  # Coordinates aligned to another file's,
                                # or to anatomical "truth"
         "3" = "Talairach",     # Coordinates aligned to Talairach-
                                # Tournoux Atlas; (0,0,0)=AC, etc.
         "4" = "MNI_152")       # MNI 152 normalized coordinates
}

convert.units <- function(units) {
  ## defgroup NIFTI1_UNITS
  ## brief nifti1 units codes to describe the unit of measurement for
  ## each dimension of the dataset
  switch(as.character(units),
         "0" = "Unknown",       # unspecified units
         "1" = "meter",         # meters
         "2" = "mm",            # millimeters
         "3" = "micron",        # micrometers
         "8" = "sec",           # seconds
         "16" = "msec",         # milliseconds
         "24" = "usec",         # microseconds
         "32" = "Hz",           # Hertz
         "40" = "ppm",          # parts per million
         "48" = "rads")         # radians per second
}

convert.slice <- function(sc) {
  switch(as.character(sc),
         "0" = "Unknown",
         "1" = "Seq_Inc",       # sequential increasing
         "2" = "Seq_Dec",       # sequential decreasing
         "3" = "Alt_Inc",       # alternating increasing
         "4" = "Alt_Dec",       # alternating decreasing
         "5" = "Alt_Inc2",      # alternating increasing #2
         "6" = "Alt_Dec2")      # alternating decreasing #2
}

## Bitwise conversion subroutines
xyzt2space <- function (xyzt) {
  ## define XYZT_TO_SPACE(xyzt)       ( (xyzt) & 0x07 )
  require("bitops")
  bitAnd(xyzt, 7)
}

xyzt2time <- function (xyzt) {
  ## define XYZT_TO_TIME(xyzt)        ( (xyzt) & 0x38 )
  require("bitops")
  bitAnd(xyzt, 56)
}

dim2freq <- function(di) {
  ## define DIM_INFO_TO_FREQ_DIM(di)   ( ((di)     ) & 0x03 )
  require("bitops")
  bitAnd(di, 3)
}

dim2phase <- function(di) {
  ## define DIM_INFO_TO_PHASE_DIM(di)  ( ((di) >> 2) & 0x03 )
  require("bitops")
  bitAnd(bitShiftR(di, 2), 3)
}

dim2slice <- function(di) {
  ## define DIM_INFO_TO_SLICE_DIM(di)  ( ((di) >> 4) & 0x03 )
  require("bitops")
  bitAnd(bitShiftR(di, 4), 3)
}

quaternion2rotation <- function(b, c, d) {
  a <- sqrt(1 - (b*b+c*c+d*d))
  R <- matrix(c(a*a+b*b-c*c-d*d, 2*b*c+2*a*d, 2*b*d-2*a*c,  # column 1
                2*b*c-2*a*d, a*a+c*c-b*b-d*d, 2*c*d+2*a*b,  # column 2
                2*b*d+2*a*c, 2*c*d-2*a*b, a*a+d*d-c*c-b*b), # column 3
              3, 3)
  return(R)
}

#############################################################################
## image() for class="nifti"
#############################################################################

#if (!isGeneric("image")) {
#  if (is.function("image"))
#     fun <- image
#  else
#    fun <- function(object) { standardGeneric("image") }
#  setGeneric("image", fun)
#}
setMethod("image", signature(x="nifti"),
          function(x, z=1, w=1, col=gray(0:64/64),
                   plot.type=c("multiple","single"), zlim=NULL,
                   xlab="", ylab="", axes=FALSE, ...) {
            X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
            if (X == 0 || Y == 0 || Z == 0)
              stop("size of NIfTI volume is zero, nothing to plot")
            if (is.null(zlim))
              zlim <- c(x@"cal_min", x@"cal_max")
            if (is.na(W)) {
              ## three-dimensional array
              if (z < 1 || z > Z)
                stop("slice \"z\" out of range")
              oldpar <- par()
              if (plot.type[1] == "multiple")
                index <- 1:Z
              else
                index <- z
              lz <- length(index)
              par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
              for (z in index)
                graphics:::image(1:X, 1:Y, x[,,z], col=col, zlim=zlim,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
            } else {
              ## four-dimensional array
              if (w < 1 || w > W)
                stop("volume \"w\" out of range")
              if (z < 1 || z > Z)
                stop("slice \"z\" out of range")
              if (plot.type[1] == "multiple")
                index <- 1:Z
              else
                index <- z
              lz <- length(index)
              par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
              for (z in index)
                graphics:::image(1:X, 1:Y, x[,,z,w], col=col, zlim=zlim,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
            }
            
          })

#############################################################################
## overlay() for class="nifti"
#############################################################################

setGeneric("overlay", function(x, y, ...) standardGeneric("overlay"))
setMethod("overlay", signature(x="nifti", y="nifti"),
          function(x, y, z=1, w=1, col.x=gray(0:64/64), col.y=hotmetal(),
                   zlim.x=NULL, zlim.y=NULL, plot.type=c("multiple","single"),
                   xlab="", ylab="", axes=FALSE, ...) {
            ## should do plot.type=c("multiple","single") like plot.ts
            if (!all(dim(x)[1:3] == dim(y)[1:3]))
              stop("dimensions of \"x\" and \"y\" must be equal")
            X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
            if (X == 0 || Y == 0 || Z == 0)
              stop("size of NIfTI volume is zero, nothing to plot")
            if (is.null(zlim.x))
              zlim.x <- c(x@"cal_min", x@"cal_max")
            if (is.null(zlim.y))
              zlim.y <- c(y@"cal_min", y@"cal_max")
            if (is.na(W)) {
              ## three-dimensional array
              oldpar <- par()
              if (plot.type[1] == "single") {
                if (z < 1 || z > Z)
                  stop("slice \"z\" out of range")
                par(mfrow=c(1,1), mar=rep(0,4))
                graphics:::image(1:X, 1:Y, x[,,z], col=col.x, zlim=zlim.x,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
                graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                 add=TRUE)
              } else {
                par(mfrow=c(ceiling(sqrt(Z)), ceiling(sqrt(Z))), mar=rep(0,4))
                for (z in 1:Z) {
                  graphics:::image(1:X, 1:Y, x[,,z], col=col.x, zlim=zlim.x,
                                   axes=axes, xlab=xlab, ylab=ylab, ...)
                  graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                   add=TRUE)
                }
              }
            } else {
              ## four-dimensional array
              oldpar <- par()
              if (w < 1 || w > W)
                stop("volume \"w\" out of range")
              if (z < 1 || z > Z)
                stop("slice \"z\" out of range")
              if (plot.type == "single") {
                par(mfrow=c(1,1), mar=rep(0,4))
                graphics:::image(1:X, 1:Y, x[,,z,w], col=col.x, zlim=zlim.x,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
                graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                 add=TRUE)
              } else {
                par(mfrow=c(ceiling(sqrt(Z)), ceiling(sqrt(Z))), mar=rep(0,4))
                for (z in 1:Z) {
                  graphics:::image(1:X, 1:Y, x[,,z,w], col=col.x, zlim=zlim.x,
                                   axes=axes, xlab=xlab, ylab=ylab, ...)
                  graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                   add=TRUE)
                }
              }
              graphics:::image(1:X, 1:Y, x[,,z,w], col=colx, zlim=zlimx,
                               axes=axes, xlab=xlab, ylab=ylab, ...)
              graphics:::image(1:X, 1:Y, y[,,z], col=coly, zlim=zlimy,
                               add=TRUE)
            }
            par(oldpar)
            invisible()
          })

#############################################################################
## orthographic() for class="nifti"
#############################################################################

setGeneric("orthographic", function(x, ...) standardGeneric("orthographic"))
setMethod("orthographic", signature(x="nifti"),
          function(x, xyz=NULL, crosshairs=TRUE, col.crosshairs="red", w=1,
                   zlim=NULL, col=gray(0:64/64), xlab="", ylab="",
                   axes=FALSE, ...) {
            X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
            if (is.null(xyz))
              xyz <- ceiling(c(X,Y,Z)/2)
            if (X == 0 || Y == 0 || Z == 0)
              stop("size of NIfTI volume is zero, nothing to plot")
            if (is.null(zlim))
              zlim <- c(x@"cal_min", x@"cal_max")
            if (is.na(W)) {
              ## three-dimensional array
              oldpar <- par()
              par(mfrow=c(2,2), mar=rep(0,4))
              graphics::image(1:X, 1:Z, x[,xyz[2],], col=col,
                              asp=x@pixdim[4]/x@pixdim[2],
                              xlab=ylab, ylab=xlab, axes=axes, ...)
              if (crosshairs)
                abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
              graphics::image(1:Y, 1:Z, x[xyz[1],,], col=col,
                              asp=x@pixdim[4]/x@pixdim[3],
                              xlab=xlab, ylab=ylab, axes=axes, ...)
              if (crosshairs)
                abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
              graphics::image(1:X, 1:Y, x[,,xyz[3]], col=col,
                              asp=x@pixdim[3]/x@pixdim[2],
                              xlab=xlab, ylab=ylab, axes=axes, ...)
              if (crosshairs)
                abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
            } else {
              ## four-dimensional array
              if (w < 1 || w > W)
                stop("volume \"w\" out of range")
              oldpar <- par()
              par(mfrow=c(2,2), mar=rep(0,4))
              graphics::image(1:X, 1:Z, x[,xyz[2],,w], col=col,
                              asp=x@pixdim[4]/x@pixdim[2],
                              xlab=ylab, ylab=xlab, axes=axes, ...)
              if (crosshairs)
                abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
              graphics::image(1:Y, 1:Z, x[xyz[1],,,w], col=col,
                              asp=x@pixdim[4]/x@pixdim[3],
                              xlab=xlab, ylab=ylab, axes=axes, ...)
              if (crosshairs)
                abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
              graphics::image(1:X, 1:Y, x[,,xyz[3],w], col=col,
                              asp=x@pixdim[3]/x@pixdim[2],
                              xlab=xlab, ylab=ylab, axes=axes, ...)
              if (crosshairs)
                abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
            }
            par(oldpar)
            invisible()
          })
