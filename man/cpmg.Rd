% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id:$

\name{CPMG}
\alias{T2.lm}
\alias{T2.fast}
\alias{T2.fast,array-method}
\alias{T2.fast,anlz-method}
\alias{T2.fast,nifti-method}
\title{Quantitative T2}
\description{
  ...
}
\usage{
T2.lm(signal, TE, guess, control=nls.lm.control())
\S4method{T2.fast}{array}(cpmg, cpmg.mask, TE, control=nls.lm.control(maxiter=150),
         multicore=FALSE, verbose=FALSE)
}
\arguments{
  \item{signal}{is the vector of signal intensities as a function of
    echo times.}
  \item{TE}{is the vector of echo times (in seconds).}
  \item{guess}{is the vector of initial values for the parameters of
    interest: \eqn{\rho}{rho} and \eqn{T2}{T2}.}
  \item{control}{An optional list of control settings for \code{nls.lm}.
    See \code{nls.lm.control} for the names of the settable control
    values and their effect.}
  \item{cpmg}{is a multidimensional array of signal intensities.  The
    last dimension is assumed to be a function of the echo times, while
    the previous dimenions are assued to be spatial.}
  \item{cpmg.mask}{is a (logical) multidimensional array that
    identifies the voxels to be analyzed.}
  \item{multicore}{is a logical variable (default = \code{FALSE}) that
    allows parallel processing via \pkg{multicore}.}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
}
\details{
  ...
}
\value{
  A list structure is produced with (all or some of the) parameter
  estimates
  \item{\eqn{\rho}{rho}}{Scaling factor between signal intensity and T2
    (proton density).}
  \item{T2}{T2 relaxation time.}
}
%\references{}
\seealso{
  \code{\link{R1.fast}}, \code{\link{R10.lm}}
}
\examples{
## Example?
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{misc}
