\name{Bitwise Conversion Subroutines}
\alias{xyzt2space}
\alias{xyzt2time}
\alias{dim2freq}
\alias{dim2phase}
\alias{dim2slice}
\title{Bitwise Conversion Subroutines}
\description{
  Units of spatial and temporal dimensions and MRI-specific spatial and
  temporal informations.
}
\usage{
xyzt2space(xyzt)
xyzt2time(xyzt)
dim2freq(di)
dim2phase(di)
dim2slice(di)
}
\arguments{
  \item{xyzt}{}
  \item{di}{}
}
\value{
  
}
\details{
  The codes below can be used in \kbd{xyzt_units} to indicate the units
  of pixdim.  As noted earlier, dimensions 1,2,3 are for x,y,z;
  dimension 4 is for time (t).
  \itemize{
    \item{If dim[4]=1 or dim[0] < 4, there is no time axis.}
    \item{A single time series (no space) would be specified with
      \itemize{
	\item{dim[0] = 4 (for scalar data) or dim[0] = 5 (for vector
	  data)}
	\item{dim[1] = dim[2] = dim[3] = 1}
	\item{dim[4] = number of time points}
	\item{pixdim[4] = time step}
	\item{\kbd{xyzt_units} indicates units of pixdim[4]}
	\item{dim[5] = number of values stored at each time point}
      }
    }
  }
    
  Bits 0..2 of xyzt_units specify the units of pixdim[1..3] (e.g.,
  spatial units are values 1..7).  Bits 3..5 of \kbd{xyzt_units} specify
  the units of pixdim[4] (e.g., temporal units are multiples of 8).

  This compression of 2 distinct concepts into 1 byte is due to the
  limited space available in the 348 byte ANALYZE 7.5 header.  The
  macros \code{xyzt2space} and \code{xyzt2time} can be used to mask off
  the undesired bits from the \kbd{xyzt_units} fields, leaving
  \dQuote{pure} space and time codes.  Inversely, the macro
  \code{spacetime2xyzt} can be used to assemble a space code
  (0,1,2,...,7) with a time code (0,8,16,32,...,56) into the combined
  value for \kbd{xyzt_units}. 

  A few fields are provided to store some extra information that is
  sometimes important when storing the image data from an FMRI time
  series experiment.  (After processing such data into statistical
  images, these fields are not likely to be useful.)

  \itemize{
    \item{freq_dim = These fields encode which spatial dimension (1,2,
      or 3)}
    \item{phase_dim = corresponds to which acquisition dimension for MRI
      data.}
    \item{slice_dim = Examples:
      \itemize{
	\item{Rectangular scan multi-slice EPI: freq_dim = 1, phase_dim
	  = 2, slice_dim = 3 (or some permutation)}
	\item{Spiral scan multi-slice EPI: freq_dim = phase_dim = 0,
	  slice_dim = 3 (since the concepts of frequency- and
	  phase-encoding directions don't apply to spiral scan)}
      }
    }
  }

  The macros \code{dim2freq}, \code{dim2phase}, and \code{dim2slice} can
  be used to extract these values from the \kbd{dim_info} byte.
}
\references{
  Neuroimaging Informatices Technology Initiative
  \url{http://nifti.nimh.nih.gov/}
}
\seealso{
  \code{\link{convert.units}}, \code{\link{convert.slice}}
}
%\examples{
%}
\author{B. Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
