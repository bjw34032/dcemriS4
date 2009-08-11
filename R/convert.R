convert.datatype <- function(datatype) {
  switch(as.character(datatype),
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
