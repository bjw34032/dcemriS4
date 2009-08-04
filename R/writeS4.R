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
## Time-stamp: <2009-07-18 09:09:06 (bjw34032)>
## $ Id: $
##

writeNIfTI <- function(nim, filename, gzipped=TRUE, verbose=FALSE, warn=-1) {
  ## Extensions?
  if (nim@"extender"[1] > 0 || nim@"vox_offset" > 352) {
    stop("-- extensions are not currently supported! --")
  }
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Basic error checking
  if (!validObject(nim))
    stop("-- aborting writeNIfTI --")
  ## Write header file...
  if (gzipped) {
    fid <- gzfile(paste(filename, "nii.gz", sep="."), "wb")
  } else {
    fid <- file(paste(filename, "nii", sep="."), "wb")
  }

  writeBin(as.integer(nim@"sizeof_hdr"), fid, size=4)
  writeChar(nim@"data_type", fid, nchar=10, eos=NULL)
  writeChar(nim@"db_name", fid, nchar=18, eos=NULL)
  writeBin(as.integer(nim@"extents"), fid, size=4)
  writeBin(as.integer(nim@"session_error"), fid, size=2)
  writeChar(nim@"regular", fid, nchar=1, eos=NULL)
  writeBin(as.integer(nim@"dim_info"), fid, size=1)
  writeBin(as.integer(nim@"dim_"), fid, size=2)
  writeBin(nim@"intent_p1", fid, size=4)
  writeBin(nim@"intent_p2", fid, size=4)
  writeBin(nim@"intent_p3", fid, size=4)
  writeBin(as.integer(nim@"intent_code"), fid, size=2)
  writeBin(as.integer(nim@"datatype"), fid, size=2)
  writeBin(as.integer(nim@"bitpix"), fid, size=2)
  writeBin(as.integer(nim@"slice_start"), fid, size=2)
  writeBin(nim@"pixdim", fid, size=4)
  writeBin(nim@"vox_offset", fid, size=4)
  writeBin(nim@"scl_slope", fid, size=4)
  writeBin(nim@"scl_inter", fid, size=4)
  writeBin(as.integer(nim@"slice_end"), fid, size=2)
  writeBin(as.integer(nim@"slice_code"), fid, size=1)
  writeBin(as.integer(nim@"xyzt_units"), fid, size=1)
  writeBin(nim@"cal_max", fid, size=4)
  writeBin(nim@"cal_min", fid, size=4)
  writeBin(nim@"slice_duration", fid, size=4)
  writeBin(nim@"toffset", fid, size=4)
  writeBin(as.integer(nim@"glmax"), fid, size=4)
  writeBin(as.integer(nim@"glmin"), fid, size=4)
  writeChar(nim@"descrip", fid, nchar=80, eos=NULL)
  writeChar(nim@"aux_file", fid, nchar=24, eos=NULL)
  writeBin(as.integer(nim@"qform_code"), fid, size=2)
  writeBin(as.integer(nim@"sform_code"), fid, size=2)
  writeBin(nim@"quatern_b", fid, size=4)
  writeBin(nim@"quatern_c", fid, size=4)
  writeBin(nim@"quatern_d", fid, size=4)
  writeBin(nim@"qoffset_x", fid, size=4)
  writeBin(nim@"qoffset_y", fid, size=4)
  writeBin(nim@"qoffset_z", fid, size=4)
  writeBin(nim@"srow_x", fid, size=4)
  writeBin(nim@"srow_y", fid, size=4)
  writeBin(nim@"srow_z", fid, size=4)
  writeChar(nim@"intent_name", fid, nchar=16, eos=NULL)
  writeChar(nim@"magic", fid, nchar=4, eos=NULL)
  writeBin(as.integer(nim@"extender"), fid, size=1)
  ## Write image file...
  dims <- 2:(1+nim@"dim_"[1])
  if (nim@"qform_code" <= 0 && nim@"sform_code" <= 0) {
    if (verbose)
      cat("  dims =", nim@"dim_"[dims], fill=TRUE)
    writeBin(as.vector(nim), fid, size=nim@"bitpix"/8)
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
      if (verbose) {
        cat("  NIfTI-1: qform_code > 0", fill=TRUE)
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
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
        writeBin(as.vector(nim@.Data[order(index.xyz)]), fid,
                 size=nim@"bitpix"/8)
        ## nim@.Data <- array(data[order(index.xyz)], nim@"dim_"[dims])
      } else {
        stop("-- rotation matrix is NOT diagonal with +/- 1s --")
      }
      ## stop("-- qform_code > 0 not implemented --")
    }
    ## check for sform codes
    if (nim@"sform_code" > 0) {
      if (verbose) {
        cat("  NIfTI-1: sform_code > 0", fill=TRUE)
        cat("  dims =", nim@"dim_"[dims], fill=TRUE)
      }
      xyz <- matrix(0, length(nim@.Data), 3)
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
      writeBin(as.vector(nim@.Data[order(index.xyz)]), fid,
               size=nim@"bitpix"/8)
    }
  }

  close(fid)
  ## Warnings?
  options(warn=oldwarn)
  invisible()
}

############################################################################
############################################################################
############################################################################

writeANALYZE <- function(aim, filename, gzipped=TRUE, verbose=FALSE, warn=-1) {
  ## Warnings?
  oldwarn <- options()$warn
  options(warn=warn)
  ## Basic error checking
  if (!validObject(aim))
    stop("Ouch!")
  ## Write header file...
  if (gzipped)
    fid <- gzfile(paste(fname, ".hdr.gz", sep=""), "wb")
  else
    fid <- file(paste(fname, ".hdr", sep=""), "wb")
  
  writeBin(as.integer(aim@"sizeof_hdr"), fid, size=4)
  writeChar(aim@"db_type", fid, nchar=10, eos=NULL)
  writeChar(aim@"db_name", fid, nchar=18, eos=NULL)
  writeBin(as.integer(aim@"extents"), fid, size=4)
  writeBin(as.integer(aim@"session_error"), fid, size=2)
  writeChar(aim@"regular", fid, nchar=1, eos=NULL)
  writeChar(aim@"hkey_un0", fid, nchar=1, eos=NULL)
  writeBin(as.integer(aim@"dim_"), fid, size=2)
  writeChar(aim@"vox_units", fid, nchar=4, eos=NULL)
  writeChar(aim@"cal_units", fid, nchar=8, eos=NULL)
  writeBin(as.integer(aim@"unused1"), fid, size=2)
  writeBin(as.integer(aim@"datatype"), fid, size=2)
  writeBin(as.integer(aim@"bitpix"), fid, size=2)
  writeBin(as.integer(aim@"dim_un0"), fid, size=2)
  writeBin(aim@"pixdim", fid, size=4)
  writeBin(aim@"vox_offset", fid, size=4)
  writeChar(aim@"vox_units", fid, nchar=4, eos=NULL)
  writeChar(aim@"cal_units", fid, nchar=8, eos=NULL)
  writeBin(aim@"cal_max", fid, size=4)
  writeBin(aim@"cal_min", fid, size=4)
  writeBin(aim@"compressed", fid, size=4)
  writeBin(aim@"verified", fid, size=4)
  writeBin(as.integer(aim@"glmax"), fid, size=4)
  writeBin(as.integer(aim@"glmin"), fid, size=4)
  writeChar(aim@"descrip", fid, nchar=80, eos=NULL)
  writeChar(aim@"aux_file", fid, nchar=24, eos=NULL)
  writeChar(aim@"orient", fid, nchar=1, eos=NULL)
  writeChar(aim@"originator", fid, nchar=10, eos=NULL)
  writeChar(aim@"generated", fid, nchar=10, eos=NULL)
  writeChar(aim@"scannum", fid, nchar=10, eos=NULL)
  writeChar(aim@"patient_id", fid, nchar=10, eos=NULL)
  writeChar(aim@"exp_date", fid, nchar=10, eos=NULL)
  writeChar(aim@"exp_time", fid, nchar=10, eos=NULL)
  writeChar(aim@"hist_un0", fid, nchar=3, eos=NULL)
  writeBin(as.integer(aim@"views"), fid, size=4)
  writeBin(as.integer(aim@"vols_added"), fid, size=4)
  writeBin(as.integer(aim@"start_field"), fid, size=4)
  writeBin(as.integer(aim@"field_skip"), fid, size=4)
  writeBin(as.integer(aim@"omax"), fid, size=4)
  writeBin(as.integer(aim@"omin"), fid, size=4)
  writeBin(as.integer(aim@"smax"), fid, size=4)
  writeBin(as.integer(aim@"smin"), fid, size=4)
  close(fid)
  
  ## Write image file...
  if (gzipped) {
    fid <- gzfile(paste(filename, "img.gz", sep="."), "wb")
  } else {
    fid <- file(paste(filename, "img", sep="."), "wb")
  }
  dims <- 2:(1+aim@"dim_"[1])
  if (verbose)
    cat("  dims =", aim@"dim_"[dims], fill=TRUE)
  writeBin(as.vector(aim), fid, size=aim@"bitpix"/8)
  close(fid)
  invisible()
}
