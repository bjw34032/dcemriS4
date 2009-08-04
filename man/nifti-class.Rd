% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-03 16:29:09 (bjw34032)>

\name{nifti-class}
\docType{class}
\alias{nifti-class}
\title{Class "nifti" for Medical Imaging Data}
\description{
  This class encapsulates both header and image information from the
  NIfTI data standard.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nifti",
    \dots)}, but most often as the result of a call to
  \code{\link{nifti}}.
}
\section{Slots}{
  \describe{
    \item{\code{sizeof_hdr}:}{Object of class \code{"numeric"}.  MUST be
      348},
    \item{\code{data_type}:}{Object of class \code{"character"}.
      ++UNUSED++}, 
    \item{\code{db_name}:}{Object of class \code{"character"}.
      ++UNUSED++}, 
    \item{\code{extents}:}{Object of class \code{"numeric"}.
      ++UNUSED++}, 
    \item{\code{session_error}:}{Object of class \code{"numeric"}.
      ++UNUSED++}, 
    \item{\code{regular}:}{Object of class \code{"character"}.
      ++UNUSED++}, 
    \item{\code{dim_info}:}{Object of class \code{"numeric"}.  MRI slice
      ordering.}, 
    \item{\code{dim_}:}{Object of class \code{"vector"}.  Data array
      dimensions.}, 
    \item{\code{intent_p1}:}{Object of class \code{"numeric"}.  1st
      intent parameter.}, 
    \item{\code{intent_p2}:}{Object of class \code{"numeric"}.  2nd
      intent parameter.},
    \item{\code{intent_p3}:}{Object of class \code{"numeric"}.  3rd
      intent parameter.}, 
    \item{\code{intent_code}:}{Object of class \code{"numeric"}.
      NIFTIINTENT code.}, 
    \item{\code{datatype}:}{Object of class \code{"numeric"}.  Defines
      data type!}, 
    \item{\code{bitpix}:}{Object of class \code{"numeric"}.  Number
      bits/voxel.}, 
    \item{\code{slice_start}:}{Object of class \code{"numeric"}.  First
      slice index.}, 
    \item{\code{pixdim}:}{Object of class \code{"vector"}.  Grid
      spacings.}, 
    \item{\code{vox_offset}:}{Object of class \code{"numeric"}.  Offset
      into .nii file}, 
    \item{\code{scl_slope}:}{Object of class \code{"numeric"}.  Data
      scaling: slope.}, 
    \item{\code{scl_inter}:}{Object of class \code{"numeric"}.  Data
      scaling: offset.}, 
    \item{\code{slice_end}:}{Object of class \code{"numeric"}.  Last
      slice index.}, 
    \item{\code{slice_code}:}{Object of class \code{"numeric"}.  Slice
      timing order.}, 
    \item{\code{xyzt_units}:}{Object of class \code{"numeric"}.  Units
      of \code{pixdim[1..4]}}, 
    \item{\code{cal_max}:}{Object of class \code{"numeric"}.  Max
      display intensity},
    \item{\code{cal_min}:}{Object of class \code{"numeric"}.  Min
      display intensity}, 
    \item{\code{slice_duration}:}{Object of class \code{"numeric"}.
      Time for 1 slice.}, 
    \item{\code{toffset}:}{Object of class \code{"numeric"}.  Time axis
      shift.}, 
    \item{\code{glmax}:}{Object of class \code{"numeric"}.  ++UNUSED++},
    \item{\code{glmin}:}{Object of class \code{"numeric"}.  ++UNUSED++},
    \item{\code{descrip}:}{Object of class \code{"character"}.  Any text
      you like.}, 
    \item{\code{aux_file}:}{Object of class \code{"character"}.
      Auxiliary filename.},
    \item{\code{qform_code}:}{Object of class \code{"numeric"}.
      NIFTIXFORM code.}, 
    \item{\code{sform_code}:}{Object of class \code{"numeric"}.
      NIFTIXFORM code.}, 
    \item{\code{quatern_b}:}{Object of class \code{"numeric"}.
      Quaternion b param.}, 
    \item{\code{quatern_c}:}{Object of class \code{"numeric"}.
      Quaternion c param.}, 
    \item{\code{quatern_d}:}{Object of class \code{"numeric"}.
      Quaternion d param.}, 
    \item{\code{qoffset_x}:}{Object of class \code{"numeric"}.
      Quaternion x shift.}, 
    \item{\code{qoffset_y}:}{Object of class \code{"numeric"}.
      Quaternion y shift.}, 
    \item{\code{qoffset_z}:}{Object of class \code{"numeric"}.
      Quaternion z shift.}, 
    \item{\code{srow_x}:}{Object of class \code{"vector"}.  1st row
      affine transform.}, 
    \item{\code{srow_y}:}{Object of class \code{"vector"}.  2nd row
      affine transform.}, 
    \item{\code{srow_z}:}{Object of class \code{"vector"}.  3rd row
      affine transform.}, 
    \item{\code{intent_name}:}{Object of class \code{"character"}.  Name
      or meaning of data.}, 
    \item{\code{magic}:}{Object of class \code{"character"}.  MUST be
      "ni1\0" or "n+1\0".}, 
    \item{\code{extender}:}{Object of class \code{"vector"}.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "nifti")}: Display object briefly.}
  }
}
\keyword{classes}
