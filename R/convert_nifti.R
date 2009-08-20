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

############################################################################
## NIfTI list structures
############################################################################

nifti.datatypes <- list(UINT8 = 2,
                        INT16 = 4,
                        INT32 = 8,
                        FLOAT32 = 16,
                        COMPLEX64 = 32,
                        FLOAT64 = 64,
                        RGB24 = 128,
                        INT8 = 256,
                        UINT16 = 512,
                        UINT32 = 768,
                        INT64 = 1024,
                        UINT64 = 1280,
                        FLOAT128 = 1536,
                        COMPLEX128 = 1792,
                        COMPLEX256 = 2048,
                        RGBA32 = 2304)

nifti.bitpix <- list(UINT8 = 8,
                     INT8 = 8,
                     UINT16 = 16,
                     INT16 = 16,
                     UINT32 = 32,
                     INT32 = 32,
                     UINT64 = 64,
                     INT64 = 64,
                     FLOAT32 = 32,
                     FLOAT64 = 64,
                     FLOAT128 = 128,
                     COMPLEX64 = 64,
                     COMPLEX128 = 128,
                     COMPLEX256 = 256,
                     RGB24 = 24,
                     RGBA32 = 32)

nifti.intent.code <- list("None" = 0,
                          "Correl" = 2,
                          "Ttest" = 3,
                          "Ftest" = 4,
                          "Zscore" = 5,
                          "Chisq" = 6,
                          "Beta" = 7,
                          "Binom" = 8,
                          "Gamma" = 9,
                          "Poisson" = 10,
                          "Normal" = 11,
                          "Ftest_Nonc" = 12,
                          "Chisq_Nonc" = 13,
                          "Logistic" = 14,
                          "Laplace" = 15,
                          "Uniform" = 16,
                          "Ttest_Nonc" = 17,
                          "Weibull" = 18,
                          "Chi" = 19,
                          "Invgauss" = 20,
                          "Extval" = 21,
                          "Pval" = 22,
                          "Logpval" = 23,
                          "Log10pval" = 24,
                          "Estimate" = 1001,  # estimate of some parameter
                          "Label" = 1002,     # index into some set of labels
                          "Neuroname" = 1003, # index into the NeuroNames
                                              # labels set
                          "Genmatrix" = 1004, # M x N matrix at each voxel
                          "Symmatrix" = 1005, # N x N symmetric matrix at
                                              # each voxel
                          "Dispvect" = 1006,  # a displacement field
                          "Vector" = 1007,    # a displacement vector
                          "Pointset" = 1008,  # a spatial coordinate
                          "Triangle" = 1009,  # triple of indexes
                          "Quaternion" = 1010, # a quaternion
                          "Dimless" = 1011)   # Dimensionless value - no params

############################################################################
## Conversion subroutines
############################################################################

convert.datatype <- function(datatype)
  names(which(nifti.datatypes == datatype))

convert.intent <- function(intent.code)
  names(which(nifti.intent.code == intent.code))

convert.form <- function(form.code) {
  switch(as.character(form.code),
         "0" = "Unknown",       # Arbitrary coordinates (Method 1)
         "1" = "Scanner_Anat",  # Scanner-based anatomical coordinates
         "2" = "Aligned_Anat",  # Coordinates aligned to another file's,
                                # or to anatomical "truth"
         "3" = "Talairach",     # Coordinates aligned to Talairach-
                                # Tournoux Atlas; (0,0,0)=AC, etc.
         "4" = "MNI_152")       # MNI 152 normalized coordinates
}

convert.units <- function(units) {
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

convert.slice <- function(slice.code) {
  switch(as.character(slice.code),
         "0" = "Unknown",
         "1" = "Seq_Inc",       # sequential increasing
         "2" = "Seq_Dec",       # sequential decreasing
         "3" = "Alt_Inc",       # alternating increasing
         "4" = "Alt_Dec",       # alternating decreasing
         "5" = "Alt_Inc2",      # alternating increasing #2
         "6" = "Alt_Dec2")      # alternating decreasing #2
}

############################################################################
## Bitwise conversion subroutines
############################################################################

xyzt2space <- function (xyzt) {
  ## define XYZT_TO_SPACE(xyzt) ( (xyzt) & 0x07 )
  require("bitops")
  bitAnd(xyzt, 7)
}

xyzt2time <- function (xyzt) {
  ## define XYZT_TO_TIME(xyzt) ( (xyzt) & 0x38 )
  require("bitops")
  bitAnd(xyzt, 56)
}

dim2freq <- function(diminfo) {
  ## define DIM_INFO_TO_FREQ_DIM(di) ( ((di)     ) & 0x03 )
  require("bitops")
  bitAnd(diminfo, 3)
}

dim2phase <- function(diminfo) {
  ## define DIM_INFO_TO_PHASE_DIM(di) ( ((di) >> 2) & 0x03 )
  require("bitops")
  bitAnd(bitShiftR(diminfo, 2), 3)
}

dim2slice <- function(diminfo) {
  ## define DIM_INFO_TO_SLICE_DIM(di) ( ((di) >> 4) & 0x03 )
  require("bitops")
  bitAnd(bitShiftR(diminfo, 4), 3)
}

############################################################################
## as.nifti()
############################################################################

as.nifti <- function(from, value=NULL, verbose=FALSE) {
  anlz.as.nifti <- function(from, value=nifti()) {
    ## So what kind of thing do we keep?
    slots <- c("dim_", "datatype", "bitpix", "pixdim", "descrip",
               "aux_file", ".Data")
    sapply(slots, function(x) { slot(value, x) <<- slot(from, x); NULL })
    #value@"data_type" <- convert.datatype(value@datatype)
    calset <- !(from@"cal_max" == 0 && from@"cal_min" == 0)
    value@"cal_max" <- ifelse(calset, from@"cal_min", from@"glmax")
    value@"cal_min" <- ifelse(calset, from@"cal_max", from@"glmin")
    return(value)
    ## The below code should apply the orient code as per 
    ## http://eeg.sourceforge.net/ANALYZE75.pdf and 
    ## http://eeg.sourceforge.net/AnalyzeDirect_ANALYZEformat.pdf
    ## However the NIfTI website says that this field is not often
    ##   set properly so I am unsure whether to apply this transform
    ##
    ##R<-diag(1, 4)
    ##R[,1]<- -R[,1] # as i is by default Left in the ANALYZE standard
    ##switch(from@"orient",
    ## orient:   slice orientation for this dataset.
    ##      0         transverse unflipped (i,j,k)=(L,A,S)
    #0=R,
    ##      1         coronal unflipped    (i,j,k)=(L,S,A)
    #1=R[,2:3] <- R[,3:2],
    ##      2         sagittal unflipped   (i,j,k)=(A,S,L)
    #2=R[,c(1,3)] <- R[,c(3,1)],
    ##      3         transverse flipped   (i,j,k)=(R,A,S)
    #3=R[,1]<--R[,1],
    ##      4         coronal flipped	    (i,j,k)=(R,S,A)
    #4=R[,1]<--R[,1];R[,2:3] <- R[,3:2],
    ##      5         sagittal flipped	    (i,j,k)=(A,S,R)
    #5=R[,1]<--R[,1];R[,c(1,3)] <- R[,c(3,1)],
    #stop("Unknown orientation: ",from@orient))
    #value@"xyzt_units" <- from@"vox_units"
    #value@"sform_code" <- 1
    #value@"srow_x" <- R[1,]
    #value@"srow_y" <- R[2,]
    #value@"srow_z" <- R[3,]
  }
  integertype <- function(from) {
    integer.ranges <- list("BINARY" = c(0, 1),
                           "UINT8" = c(0, 255),
                           "INT16" = c(-2^15, 2^15-1),
                           "INT32" = c(-2^31, 2^31-1))
    fromRange <- range(from)
    for (i in 1:length(integer.ranges)) {
      if (fromRange[1] >= integer.ranges[[i]][1] &&
	  fromRange[2] <= integer.ranges[[i]][2]) {
	return(names(integer.ranges)[i])
      }
    }
    warning("Range too large to be kept as integer forcing float")
    floattype(from)
  }

  floattype <- function(from)
    return("FLOAT32")

  if (is.null(value)) {
    nim <- nifti()
  } else {
    nim <- value
  }
  if (is(from, "anlz")) {
    nim <- anlz.as.nifti(from,value)
  } else if (is.array(from)) {
    ## Determine a sensible datatype
    dataClass <- class(from[1])
    datatypeString <- switch(dataClass,
                             logical = "BINARY",
                             integer = integertype(from),
                             numeric = floattype(from),
                             stop("Can't transform data in from: ",
                                  class(from[1])))
    nim@"data_type" <- datatypeString
    nim@"datatype" <- nifti.datatypes[[datatypeString]]
    nim@"bitpix" <- nifti.bitpix[[datatypeString]]
    nim@"cal_min" <- min(from)
    nim@"cal_max" <- max(from)
    nim@"dim_" <- c(length(dim(from)), dim(from))
    if (length(nim@"dim_") < 8)
      nim@"dim_" <- c(nim@"dim_", rep(1, 8 - length(nim@"dim_")))
    nim@.Data <- from
  } else {
    if (is.list(from)) {
      nim <- lapply(from, function(x) as.nifti(x, value))
      lapply(names(from),
             function(x) {
               if (is.nifti(nim[[x]])) {
                 nim[[x]]@"intent_code" <<- nifti.intent.code[["Estimate"]]
                 nim[[x]]@"intent_name" <<- substr(x, 1, 15)
               }
             })
    } else {
      if (verbose) 
        cat("Warning cannot convert class =", class(from),
            "to nifti object", fill=TRUE)
      nim <- from
    }
  }
  return(nim)
}

############################################################################
## as("anlz", "nifti")
############################################################################

setAs("anlz", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )

############################################################################
## as("array", "nifti")
############################################################################

setAs("array", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )

############################################################################
## as("list", "nifti")
############################################################################

setAs("list", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )

