% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{Bitwise Conversion Subroutines}
\alias{xyzt2space}
\alias{xyzt2time}
\alias{dim2freq}
\alias{dim2phase}
\alias{dim2slice}
\title{Bitwise Conversion Subroutines}
\description{
  Units of spatial and temporal dimensions, and MRI-specific spatial and
  temporal information.
}
\usage{
xyzt2space(xyzt)
xyzt2time(xyzt)
dim2freq(diminfo)
dim2phase(diminfo)
dim2slice(diminfo)
}
\arguments{
  \item{xyzt}{represents the units of pixdim[1..4] in the NIfTI header.}
  \item{diminfo}{represents MRI slice ordering in the NIfTI header.}
}
\details{
  The functions \code{xyzt2space} and \code{xyzt2time} can be used to
  mask off the undesired bits from the \kbd{xyzt_units} fields, leaving
  \dQuote{pure} space and time codes.
  \url{http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/xyzt_units.html}

  The functions \code{dim2freq}, \code{dim2phase}, and \code{dim2slice}
  can be used to extract values from the \kbd{dim_info} byte.
  \url{http://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/dim_info.html}
}
\value{
  For \kbd{diminfo}: the frequency, phase and slice dimensions encode
  which spatial dimension (1,2, or 3) corresponds to which acquisition
  dimension for MRI data.  For \kbd{xyzt_units}: the codes are used to
  indicate the units of pixdim.  Dimensions 1,2,3 are for x,y,z;
  dimension 4 is for time (t).
}
\references{
  Neuroimaging Informatics Technology Initiative (NIfTI)
  \url{http://nifti.nimh.nih.gov/}
}
\seealso{
  \code{\link{convert.units}}, \code{\link{convert.slice}}
}
%\examples{}
\author{B. Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
