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
## Time-stamp: <2009-07-18 09:11:25 (bjw34032)>
## $ Id: $
##

readNIfTI <- function(fname, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Check if any file extensions are present
  NIFTI <- ifelse(length(grep("nii", fname)) != 0, TRUE, FALSE)
  GZ <- ifelse(length(grep("gz", fname)) != 0, TRUE, FALSE)

  if (GZ) {
    if (NIFTI) {
      ls.nii.gz <- system(paste("ls", fname), intern=TRUE, ignore.stderr=TRUE)
      if (length(ls.nii.gz) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "\n  file =", ls.nii.gz), fill=TRUE)
        nim <- read.nifti.content(sub(".nii.gz", "", fname), gzipped=TRUE,
                                  verbose=verbose, warn=warn)
        options(warn=oldwarn)
        return(nim)
      }
    } else {
      options(warn=oldwarn)
      stop(paste(fname, "is not recognized."))
    }
  } else {
    if (NIFTI) {
      ls.nii <- system(paste("ls", fname), intern=TRUE, ignore.stderr=TRUE)
      if (length(ls.nii) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "\n  file =", ls.nii), fill=TRUE)
        nim <- read.nifti.content(sub(".nii", "", fname), gzipped=FALSE,
                                  verbose=verbose, warn=warn)
        options(warn=oldwarn)
        return(nim)
      } else {
        options(warn=oldwarn)
        stop(paste(fname, "is not recognized."))
      }
    }
  }
}

############################################################################
############################################################################
############################################################################

read.nifti.content <- function(fname, onefile=TRUE, gzipped=TRUE,
                               verbose=FALSE, warn=-1) {
  ## Open appropriate file
  if (gzipped) {
    suffix <- ifelse(onefile, "nii.gz", "hdr.gz")
    fname <- paste(fname, suffix, sep=".")
    fid <- gzfile(fname, "rb")
  } else {
    suffix <- ifelse(onefile, "nii", "hdr")
    fname <- paste(fname, suffix, sep=".")
    fid <- file(fname, "rb")
  }
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Test for endian properties
  endian <- .Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  if (sizeof.hdr != 348) {
    close(fid)
      endian <- "swap"
    if (gzipped)
      fid <- gzfile(fname, "rb")
    else
      fid <- file(fname, "rb")
    sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
    if (verbose) cat("  ENDIAN = swap", fill=TRUE)
  }
  ## Construct S4 object
  nim <- new("nifti")
  nim@"sizeof_hdr" <- sizeof.hdr
  nim@"data_type" <- rawToChar(readBin(fid, "raw", n=10))
  nim@"db_name" <- rawToChar(readBin(fid, "raw", n=18))
  nim@"extents" <- readBin(fid, integer(), size=4, endian=endian)
  nim@"session_error" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"regular" <- rawToChar(readBin(fid, "raw", n=1))
  nim@"dim_info" <- readBin(fid, integer(), size=1, signed=FALSE,
                            endian=endian)
  nim@"dim_" <- readBin(fid, integer(), 8, size=2, endian=endian)
  nim@"intent_p1" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"intent_p2" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"intent_p3" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"intent_code" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"datatype" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"bitpix" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"slice_start" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"pixdim" <- readBin(fid, numeric(), 8, size=4, endian=endian)
  nim@"vox_offset" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"scl_slope" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"scl_inter" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"slice_end" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"slice_code" <- readBin(fid, integer(), size=1, signed=FALSE,
                              endian=endian)
  nim@"xyzt_units" <- readBin(fid, integer(), size=1, signed=FALSE,
                              endian=endian)
  nim@"cal_max" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"cal_min" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"slice_duration" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"toffset" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"glmax" <- readBin(fid, integer(), size=4, endian=endian)
  nim@"glmin" <- readBin(fid, integer(), size=4, endian=endian)
  ## was data_history substruct in Analyze 7.5
  nim@"descrip" <- rawToChar(readBin(fid, "raw", n=80))
  nim@"aux_file" <- rawToChar(readBin(fid, "raw", n=24))
  nim@"qform_code" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"sform_code" <- readBin(fid, integer(), size=2, endian=endian)
  nim@"quatern_b" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"quatern_c" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"quatern_d" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"qoffset_x" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"qoffset_y" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"qoffset_z" <- readBin(fid, numeric(), size=4, endian=endian)
  nim@"srow_x" <- readBin(fid, numeric(), 4, size=4, endian=endian)
  nim@"srow_y" <- readBin(fid, numeric(), 4, size=4, endian=endian)
  nim@"srow_z" <- readBin(fid, numeric(), 4, size=4, endian=endian)
  nim@"intent_name" <- rawToChar(readBin(fid, "raw", n=16))
  nim@"magic" <- rawToChar(readBin(fid, "raw", n=4))
  ## To flag such a struct as being conformant to the NIFTI-1 spec,
  ## the last 4 bytes of the header must be either the C String "ni1"
  ## or "n+1"; in hexadecimal, the 4 bytes [6E 69 31 00] or [6E 2B 31
  ## 00] (in any future version of this format, the 1 will be upgraded
  ## to 2, etc.).  Normally, such a "magic number" or flag goes at the
  ## start of the file, but trying to avoid clobbering widely-used
  ## ANALYZE 7.5 fields led to putting this marker last.  However,
  ## recall that "the last shall be first" (Matthew 20:16).
  if (!(nim@"magic" %in% c("n+1","ni1")))
    stop(" -- Unrecognized \"magic\" field! --")
  nim@"extender" <- readBin(fid, integer(), 4, size=1, signed=FALSE, endian=endian)
  ## If extension[0] is nonzero, it indicates that extended header
  ## information is present in the bytes following the extension
  ## array.  In a .nii file, this extended header data is before the
  ## image data (and vox_offset must be set correctly to allow for
  ## this).  In a .hdr file, this extended data follows extension and
  ## proceeds (potentially) to the end of the file.
  if (nim@"extender"[1] > 0 || nim@"vox_offset" > 352) {
    nimext <- as(nim, "niftiExtension")
    i <- 1
    while (seek(fid) < nim@"vox_offset" && i < 2) {
      nimext@esize <- readBin(fid, integer(), size=4, endian=endian)
      nimext@ecode <- readBin(fid, integer(), size=4, endian=endian)
      nimext@edata <- rawToChar(readBin(fid, "raw", n=nimext@esize-8,
                                        endian=endian))
      i <- i + 1
    }
    if (seek(fid) < nim@"vox_offset")
      stop("-- multiple extensions have been detected, not supported --")
    if (seek(fid) > nim@"vox_offset")
      stop("-- extension size (esize) has overshot voxel offset --")
    nim <- nimext
  }

  n <- prod(nim@"dim_"[2:5])
  data <-
    switch(as.character(nim@"datatype"),
           "2" = readBin(fid, integer(), n, nim@"bitpix"/8, signed=FALSE,
             endian=endian),
           "4" = readBin(fid, integer(), n, nim@"bitpix"/8, endian=endian),
           "8" = readBin(fid, integer(), n, nim@"bitpix"/8, endian=endian),
           "16" = readBin(fid, numeric(), n, nim@"bitpix"/8, endian=endian),
           "64" = readBin(fid, double(), n, nim@"bitpix"/8, endian=endian),
           "512" = readBin(fid, numeric(), n, nim@"bitpix"/8, endian=endian),
           stop(paste("Data type ", nim@"datatype", " unsupported in ",
                      fname, ".img", sep=""))
           )
  close(fid)

  ## 3D IMAGE (VOLUME) ORIENTATION AND LOCATION IN SPACE:
  ## There are 3 different methods by which continuous coordinates can
  ## attached to voxels.  The discussion below emphasizes 3D volumes,
  ## and the continuous coordinates are referred to as (x,y,z).  The
  ## voxel index coordinates (i.e., the array indexes) are referred to
  ## as (i,j,k), with valid ranges:
  ##   i = 0 .. dim[1]-1
  ##   j = 0 .. dim[2]-1  (if dim[0] >= 2)
  ##   k = 0 .. dim[3]-1  (if dim[0] >= 3)
  ## The (x,y,z) coordinates refer to the CENTER of a voxel.  In
  ## methods 2 and 3, the (x,y,z) axes refer to a subject-based
  ## coordinate system, with
  ##   +x = Right  +y = Anterior  +z = Superior.
  ## This is a right-handed coordinate system.  However, the exact
  ## direction these axes point with respect to the subject depends on
  ## qform_code (Method 2) and sform_code (Method 3).

  dims <- 2:(1+nim@"dim_"[1])

  if (nim@"qform_code" <= 0 && nim@"sform_code" <= 0) {
    if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
    nim@.Data <- array(data, nim@"dim_"[dims])
  } else {
    i <- 0:(nim@"dim_"[2]-1)
    j <- 0:(nim@"dim_"[3]-1)
    k <- 0:(nim@"dim_"[4]-1)
    ijk <- cbind(rep(i, nim@"dim_"[3] * nim@"dim_"[4]),
                 rep(rep(j, each=nim@"dim_"[2]), nim@"dim_"[4]),
                 rep(k, each=nim@"dim_"[2] * nim@"dim_"[3]))
    index.ijk <- (ijk[,1] +
                  ijk[,2] * nim@"dim_"[2] +
                  ijk[,3] * nim@"dim_"[2] * nim@"dim_"[3])
    ## check for qform codes
    if (nim@"qform_code" > 0) {
      if (verbose) cat("  NIfTI-1: qform_code > 0", fill=TRUE)
      qfac <- nim@"pixdim"[1]
      R <- quaternion2rotation(nim@"quatern_b",
                               nim@"quatern_c",
                               nim@"quatern_d")
      ## HACK!!! To ensure matrix is integer-valued
      R <- round(R)
      qoffset <- c(nim@"qoffset_x", nim@"qoffset_y", nim@"qoffset_z")
      if (qfac < 0)
        R[3,3] <- -R[3,3]
      if (all(abs(R) == diag(3))) {
        ## HACK!!! Multiply x-dimension for proper orientation in R
        R[1,] <- -R[1,]
        xyz <-
          t(sweep(R %*% t(sweep(ijk, 2, as.array(nim@"pixdim"[2:4]), "*")),
                  1, as.array(qoffset), "+"))
        index.xyz <- (xyz[,1] +
                      xyz[,2] * nim@"dim_"[2] +
                      xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])      
        if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
        nim@.Data <- array(data[order(index.xyz)], nim@"dim_"[dims])
      } else {
        stop("-- rotation matrix is NOT (approximately) diagonal with +/- 1s --")
      }
    }
    ## check for sform codes
    if (nim@"sform_code" > 0) {
      if (verbose) cat("  NIfTI-1: sform_code > 0", fill=TRUE)
      xyz <- matrix(0, length(data), 3)
      xyz[,1] <- (nim@"srow_x"[1] * ijk[,1] + nim@"srow_x"[2] * ijk[,2] +
                  nim@"srow_x"[3] * ijk[,3] + nim@"srow_x"[4])
      ## HACK!!! Multiply x-dimension for proper orientation in R
      xyz[,1] <- -xyz[,1]
      xyz[,2] <- (nim@"srow_y"[1] * ijk[,1] + nim@"srow_y"[2] * ijk[,2] +
                  nim@"srow_y"[3] * ijk[,3] + nim@"srow_y"[4])
      xyz[,3] <- (nim@"srow_z"[1] * ijk[,1] + nim@"srow_z"[2] * ijk[,2] +
                  nim@"srow_z"[3] * ijk[,3] + nim@"srow_z"[4])
      index.xyz <- (xyz[,1] +
                    xyz[,2] * nim@"dim_"[2] +
                    xyz[,3] * nim@"dim_"[2] * nim@"dim_"[3])
      if (verbose) cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      nim@.Data <- array(data[order(index.xyz)], nim@"dim_"[dims])
    }
  }

  ## Warnings?
  options(warn=oldwarn)
  ## Check validity
  validObject(nim)
  return(nim)
}

############################################################################
############################################################################
############################################################################

readANALYZE <- function(fname, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Check if any file extensions are present
  ANLZ <- ifelse(length(grep("hdr|img", fname)) != 0, TRUE, FALSE)
  GZ <- ifelse(length(grep("gz", fname)) != 0, TRUE, FALSE)

  if (GZ) {
    ls.hdr.gz <- system(paste("ls", fname), intern=TRUE, ignore.stderr=TRUE)
    if (length(ls.hdr.gz) != 0) {
      if (verbose)
        cat(paste("  fname =", fname, "\n  file =", ls.hdr.gz), fill=TRUE)
      aim <- read.analyze.content(sub(".hdr.gz", "", fname), gzipped=TRUE,
                                  verbose=verbose, warn=warn)
      options(warn=oldwarn)
      return(aim)
    } else {
      options(warn=oldwarn)
      stop(paste(fname, "is not recognized."))
    }
  } else {
    ls.hdr <- system(paste("ls", fname), intern=TRUE, ignore.stderr=TRUE)
    if (length(ls.hdr) != 0) {
      if (verbose)
        cat(paste("  fname =", fname, "\n  file =", ls.hdr), fill=TRUE)
      aim <- read.analyze.content(sub(".hdr", "", fname), gzipped=FALSE,
                                  verbose=verbose, warn=warn)
      options(warn=oldwarn)
      return(aim)
    } else {
      ls.hdr.gz <- system(paste("ls ", fname, ".hdr.gz", sep=""),
                          intern=TRUE, ignore.stderr=TRUE)
      if (length(ls.hdr.gz) != 0) {
        if (verbose)
          cat(paste("  fname =", fname, "\n  file  =", ls.hdr.gz), fill=TRUE)
        aim <- read.analyze.content(fname, gzipped=TRUE, verbose=verbose,
                                    warn=warn)
        options(warn=oldwarn)
        return(aim)
      } else {
        options(warn=oldwarn)
        stop(paste(fname, "is not recognized."))
      }
    }
  }
}

############################################################################
############################################################################
############################################################################

read.analyze.content <- function(fname, gzipped=TRUE, verbose=FALSE,
                                 warn=-1) {
  ## Open header file
  if (gzipped) {
    fname <- paste(fname, "hdr.gz", sep=".")
    if (verbose) cat("  hdr =", fname, fill=TRUE)
    fid <- gzfile(fname, "rb")
  } else {
    fname <- paste(fname, "hdr", sep=".")
    if (verbose) cat("  hdr =", fname, fill=TRUE)
    fid <- file(fname, "rb")
  }
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Test for endian properties
  endian <- .Platform$endian
  sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  if (sizeof.hdr != 348) {
    close(fid)
    endian <- "swap"
    if (gzipped)
      fid <- gzfile(fname, "rb")
    else
      fid <- file(fname, "rb")
    sizeof.hdr <- readBin(fid, integer(), size=4, endian=endian)
  }
  ## Construct S4 object
  aim <- new("anlz")
  aim@"sizeof_hdr" <- sizeof.hdr
  aim@"db_type" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"db_name" <- rawToChar(readBin(fid, "raw", n=18))
  aim@"extents" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"session_error" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"regular" <- rawToChar(readBin(fid, "raw", n=1))
  aim@"hkey_un0" <- rawToChar(readBin(fid, "raw", n=1))
  aim@"dim_" <- readBin(fid, integer(), 8, size=2, endian=endian)
  aim@"vox_units" <- rawToChar(readBin(fid, "raw", n=4))
  aim@"cal_units" <- rawToChar(readBin(fid, "raw", n=8))
  aim@"unused1" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"datatype" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"bitpix" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"dim_un0" <- readBin(fid, integer(), size=2, endian=endian)
  aim@"pixdim" <- readBin(fid, numeric(), 8, size=4, endian=endian)
  aim@"vox_offset" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"funused1" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"funused2" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"funused3" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"cal_max" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"cal_min" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"compressed" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"verified" <- readBin(fid, numeric(), size=4, endian=endian)
  aim@"glmax" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"glmin" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"descrip" <- rawToChar(readBin(fid, "raw", n=80))
  aim@"aux_file" <- rawToChar(readBin(fid, "raw", n=24))
  aim@"orient" <- rawToChar(readBin(fid, "raw", n=1))
  aim@"originator" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"generated" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"scannum" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"patient_id" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"exp_date" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"exp_time" <- rawToChar(readBin(fid, "raw", n=10))
  aim@"hist_un0" <- rawToChar(readBin(fid, "raw", n=3))
  aim@"views" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"vols_added" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"start_field" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"field_skip" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"omax" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"omin" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"smax" <- readBin(fid, integer(), size=4, endian=endian)
  aim@"smin" <- readBin(fid, integer(), size=4, endian=endian)
  close(fid)
  ## Open image file
  if (gzipped) {
    fname <- paste(fname, ".img.gz", sep=".")
    fid <- gzfile(fname, "rb")
  } else {
    fname <- paste(fname, "img", sep=".")
    fid <- file(fname, "rb")
  }
  n <- prod(aim@"dim_"[2:5])
  data <- switch(as.character(aim@"datatype"),
                 "1" = readBin(fid, integer(), n, aim@"bitpix"/8, signed=FALSE, endian=endian),
                 "2" = readBin(fid, integer(), n, aim@"bitpix"/8, signed=FALSE, endian=endian),
                 "4" = readBin(fid, integer(), n, aim@"bitpix"/8, endian=endian),
                 "8" = readBin(fid, integer(), n, aim@"bitpix"/8, endian=endian),
                 "16" = readBin(fid, numeric(), n, aim@"bitpix"/8, endian=endian),
                 "64" = readBin(fid, double(), n, aim@"bitpix"/8, endian=endian),
                 stop(paste("Data type ", aim@"datatype", " (",
                            convert.datatype.anlz(aim@"datatype"), 
                            ") unsupported in", fname, sep="")))
  close(fid)

  dims <- 2:(1+aim@"dim_"[1])
  aim@.Data <- array(data, aim@"dim_"[dims])
  
  ## Warnings?
  options(warn=oldwarn)
  ## Check validity
  validObject(aim)
  return(aim)
}
