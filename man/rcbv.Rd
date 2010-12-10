% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id:$

\name{rCBV}
\alias{rCBV}
\alias{rCBV.fast}
\alias{rCBV.fast,array-method}
\alias{rCBV.fast,anlz-method}
\alias{rCBV.fast,nifti-method}
\title{Regional Cerebral Blood Volume}
\description{
  ...
}
\usage{
rCBV(Ct, Ca, time, Hf=1, rho=1)
\S4method{rCBV.fast}{array}(signal, mask, aif, time, multicore=FALSE, verbose=FALSE)
}
\arguments{
  \item{Ct}{...}
  \item{Ca}{...}
  \item{time}{is the vector of acquisition times associated with the
    dynamic data.}
  \item{Hf}{...}
  \item{rho}{...}
  \item{signal}{is a multidimensional array of signal intensities (or
    concentrations).  The last dimension is assumed to be a function of
    the acquisition times, while the previous dimenions are assued to be
    spatial.}
  \item{mask}{is a (logical) multidimensional array that identifies the
    voxels to be analyzed.}
  \item{aif}{Arterial Input Function.}
  \item{multicore}{is a logical variable (default = \code{FALSE}) that
    allows parallel processing via \pkg{multicore}.}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
}
\details{
  ...
}
\value{
  A \code{nifti} object containing the estimates of regional cerebral
  blood volume (rCBV).
}
%\references{}
%\seealso{}
%\examples{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{misc}
