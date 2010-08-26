% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id:$

\name{CPMG}
\alias{T2.lm}
\alias{T2.fast}
\title{Quantitative T2}
\description{
  ...
}
\usage{
T2.lm(signal, TE, guess, nprint=0)
\S4method{T2.fast}{nifti}(cpmg, cpmg.mask, TE, verbose=FALSE)
\S4method{T2.fast}{anlz}(cpmg, cpmg.mask, TE, verbose=FALSE)
\S4method{T2.fast}{array}(cpmg, cpmg.mask, TE, verbose=FALSE)
}
\arguments{
  \item{signal}{is the vector of signal intensities as a function of
    echo times.}
  \item{TE}{is the vector of echo times (in seconds).}
  \item{guess}{is the vector of initial values for the parameters of
    interest: \eqn{\rho}{rho} and \eqn{T2}{T2}.}
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{cpmg}{is a multidimensional array of signal intensities.  The
    last dimension is assumed to be a function of the echo times, while
    the previous dimenions are assued to be spatial.}
  \item{cpmg.mask}{is a (logical) multidimensional array that
    identifies the voxels to be analyzed.}
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
%\examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\keyword{misc}
