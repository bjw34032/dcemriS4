
nifti.datatypes <- list(
    UINT8= 2,
    INT16= 4,
    INT32= 8,
    FLOAT32= 16,
    COMPLEX64= 32,
    FLOAT64= 64,
    RGB24= 128,
    INT8= 256,
    UINT16= 512,
    UINT32= 768,
    INT64= 1024,
    UINT64= 1280,
    FLOAT128= 1536,
    COMPLEX128= 1792,
    COMPLEX256= 2048,
    RGBA32= 2304)

nifti.datatypes.bitsperpix <- list(
	UINT8=8,
	INT8=8,
	UINT16=16,
	INT16=16,
	UINT32=32,
	INT32=32,
	UINT64=64,
	INT64=64,
	FLOAT32=32,
	FLOAT64=64,
	FLOAT128=128,
	COMPLEX64=64,
	COMPLEX128=128,
	COMPLEX256=256,
	RGB24=24,
	RGBA32=32)

convert.datatype <- function(datatype) {
  names(which(nifti.datatypes == datatype))
}

convert.intent <- function(intent.code) {
  switch(as.character(intent.code),
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

as.nifti <- function(from, value=NULL) {
  integertype <- function(from) {
    intranges <- list("BINARY" = c(0,1),
                      "UINT8" = c(0,255),
                      "INT16" = c(-32768,32767),
                      "INT32" = c(-2147483648,2147483647))
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

  floattype <- function(from) {
    return("FLOAT32")
  }


  if (is.null(value)) {
    nim <- nifti()
  } else {
    nim <- value
  }
  # Determine a sensible datatype
  if (is.array(from)) {
    dataClass <- class(from[1])
    datatypeString <- switch(dataClass,
	logical="BINARY",
	integer=integertype(from),
	numeric=floattype(from),
	stop("Can't transform data in from: ",class(from[1]))
	)

    nim@"data_type" <- datatypeString
    nim@"datatype" <- nifti.datatypes[[datatypeString]]
    nim@"bitpix" <- nifti.datatypes.bitsperpix[[datatypeString]]
    nim@"cal_min" <- min(from)
    nim@"cal_max" <- max(from)
    nim@"dim_" <- c(length(dim(from)),dim(from))
    if (length(nim@"dim_") < 8 ) {
      nim@"dim_" <- c(nim@"dim_", rep(1, 8 - length(nim@"dim_")))
    }
    
    nim@.Data<-from
  } else {
    warning("not an array")
  }
  return(nim)
}

setAs("array", "nifti", function(from) { as.nifti(from) }, as.nifti)
