% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: ftm.Rd 332 2010-01-29 16:54:07Z bjw34032 $

\name{ftm}
\alias{ftm}
\alias{ftm,array-method}
\title{Fast Template Matching via Cross-Correlation}
\description{
  Motion correction and/or co-registration of three-dimensional arrays
  (medical imaging data) are performed by applying a user-defined mask
  of voxels.  Normalized cross-correlations (in 3D) are computed using
  the FFT.
}
\usage{
\S4method{ftm}{array}(input, ...)
}
\arguments{
  \item{input}{is a four-dimensional array of signal intensities.}
%  \item{mask}{is a three-dimensional array of logical values (voxels
%    to be included are \code{TRUE}, and voxels to be excluded are
%    \code{FALSE}).}
%  \item{template}{is a three dimensional array (the same structure as
%    \code{mask}) that contains the signal intensities associated with
%    the reference data.}
%  \item{plot}{is a flag that determines if a plot of the motion
%    correction as a function of the fourth dimension (default =
%    \code{FALSE}).}
  \item{...}{Additional variables passed to the \code{plot} function.}
}
\details{
  An extremely basic method of motion correction/co-registration is
  implemented by estimating \dQuote{local} cross-correlations based on a
  binary mask that is a subset of the original three-dimensional volume.
  All convolutions are preformed via the FFT (\code{\link{fft}}) and
  repetitive calculations are minimized where possible.

  Only whole-voxel translations are considered.  This does not begin to
  capture the true effects of motion in soft tissue, but we assume that
  the object of interest (e.g., tumor) is a fairly rigid structure.
  Potential extensions include rigid-body, affine and nonlinear
  registration techniques along with interploation schemes in order to
  capture intra-voxel manipulations of the data.
}
\value{
  A list of objects are returned:
  \item{out}{Motion-corrected version of the four-dimensional array.}
  \item{offset}{Translations (in 3D) for each volume in the 4D array.}
  \item{t.center}{Estimated center of the binary mask.}
}
\references{
  Lewis, J.P. (2003) Fast normalized cross-correlation.\cr
  \url{www.idiom.com/~zilla/}
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\seealso{
  \code{\link{conv.fft}}, \code{\link{find.center}}, \code{\link{shift3D}}
}
%\examples{}
%\keyword{}
