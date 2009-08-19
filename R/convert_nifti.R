
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

nifti.datatypes.bitsperpix <- list(UINT8 = 8,
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

convert.datatype <- function(datatype)
  names(which(nifti.datatypes == datatype))

nifti.intent.code <- list(
         "None" = 0,
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
         "Estimate" = 1001,      # estimate of some parameter
         "Label" = 1002,         # index into some set of labels
         "Neuroname" = 1003,     # index into the NeuroNames labels set
         "Genmatrix" = 1004,     # M x N matrix at each voxel
         "Symmatrix" = 1005,     # N x N symmetric matrix at each voxel
         "Dispvect" = 1006,      # a displacement field
         "Vector" = 1007,        # a displacement vector
         "Pointset" = 1008,      # a spatial coordinate
         "Triangle" = 1009,      # triple of indexes
         "Quaternion" = 1010,    # a quaternion
         "Dimless" = 1011)       # Dimensionless value - no params

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
  integertype <- function(from) {
    intranges <- list("BINARY" = c(0,1),
                      "UINT8" = c(0,255),
                      "INT16" = c(-2^15,2^15-1),
                      "INT32" = c(-2^31,2^31-1))
    fromRange <- range(from)
    for (i in 1:length(intranges)) {
      if (fromRange[1] >= intranges[[i]][1] &&
	  fromRange[2]<= intranges[[i]][2]) {
	return(names(intranges)[i])
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
  ## Determine a sensible datatype
  if (is.array(from)) {
    dataClass <- class(from[1])
    datatypeString <- switch(dataClass,
                             logical = "BINARY",
                             integer = integertype(from),
                             numeric = floattype(from),
                             stop("Can't transform data in from: ",
                                  class(from[1])))
    nim@"data_type" <- datatypeString
    nim@"datatype" <- nifti.datatypes[[datatypeString]]
    nim@"bitpix" <- nifti.datatypes.bitsperpix[[datatypeString]]
    nim@"cal_min" <- min(from)
    nim@"cal_max" <- max(from)
    nim@"dim_" <- c(length(dim(from)), dim(from))
    if (length(nim@"dim_") < 8 )
      nim@"dim_" <- c(nim@"dim_", rep(1, 8 - length(nim@"dim_")))
    nim@.Data <- from
  } else {
    if (is.list(from)) {
      nim <- lapply(from, function(x) as.nifti(x, value))
      lapply(names(from), function(x) {
	if (is.nifti(nim[[x]])) {
	  nim[[x]]@"intent_code" <<- nifti.intent.code[["Estimate"]]
	  nim[[x]]@"intent_name" <<- substr(x,1,15)
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

setAs("array", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )
setAs("list", "nifti",
      function(from) { as.nifti(from) },
      function(from, value) { as.nifti(from, value) } )
