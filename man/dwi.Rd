% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dwi.Rd 332 2010-01-29 16:54:07Z bjw34032 $

\name{Apparent Diffusion Coefficient}
\alias{adc.lm}
\alias{ADC.fast,array-method}
\alias{ADC.fast}
\title{Estimate the Apparent Diffusion Coefficient (ADC)}
\description{
  Estimation of apparent diffusion coefficient (ADC) values, using a
  single exponential function, is achieved through nonlinear
  optimization.
}
\usage{
adc.lm(signal, b, guess, control=nls.lm.control())
\S4method{ADC.fast}{array}(dwi, bvalues, dwi.mask,
         control=nls.lm.control(maxiter=150), multicore=FALSE, verbose=FALSE)
}
\arguments{
  \item{signal}{Signal intensity vector as a function of b-values.}
  \item{b,bvalues}{Diffusion weightings (b-values).}
  \item{guess}{Initial values of \eqn{S_0}{S0} and \eqn{D}.}
  \item{control}{An optional list of control settings for \code{nls.lm}.
    See \code{nls.lm.control} for the names of the settable control
    values and their effect.}
  \item{dwi}{Multidimensional array of diffusion-weighted images.}
  \item{dwi.mask}{Logical array that defines the voxels to be analyzed.}
  \item{multicore}{is a logical variable (default = \code{FALSE}) that
    allows parallel processing via \pkg{multicore}.}
  \item{verbose}{Additional information will be printed when
    \code{verbose=TRUE}.}
}
\details{
  The \code{adc.lm} function estimates parameters for a vector of
  observed MR signal intensities using the following relationship
  \deqn{S(b) = S_0 \exp(-bD),}
  where \eqn{S_0}{S0} is the baseline signal intensity and \eqn{D} is
  the apparent diffusion coefficient (ADC).  It requires the routine
  \code{nls.lm} that applies the Levenberg-Marquardt algorithm.  Note,
  low b-values (\eqn{<50} or \eqn{<100} depending on who you read)
  should be avoided in the parameter estimation because they do not
  represent information about the diffusion of water in tissue.

  The \code{ADC.fast} function rearranges the assumed multidimensional
  (2D or 3D) structure of the DWI data into a single matrix to take
  advantage of internal R functions instead of loops, and called
  \code{adc.lm}.
}
\value{
  A list structure is produced with estimates of \eqn{S_0}, \eqn{D} and
  information about the convergence of the nonlinear optimization
  routine.
}
\references{
  Buxton, R.B. (2002) \emph{Introduction to Functional Magnetic
    Resonance Imaging: Principles & Techniques}, Cambridge University
  Press: Cambridge, UK.

  Callahan, P.T. (2006) \emph{Principles of Nuclear Magnetic Resonance
    Microscopy}, Clarendon Press: Oxford, UK.
  
  Koh, D.-M. and Collins, D.J. (2007) Diffusion-Weighted MRI in the
  Body: Applications and Challenges in Oncology, \emph{American Journal
    of Roentgenology}, \bold{188}, 1622-1635.
}
\seealso{
  \code{\link[minpack.lm]{nls.lm}}
}
\examples{
S0 <- 10
b <- c(0, 50, 400, 800)  # units?
D <- 0.7e-3              # mm^2 / sec (normal white matter)

## Signal intensities based on the (simplified) Bloch-Torry equation
dwi <- function(S0, b, D) {
  S0 * exp(-b*D)
}

set.seed(1234)
signal <- array(dwi(S0, b, D) + rnorm(length(b), sd=0.15),
                c(rep(1,3), length(b)))
ADC <- ADC.fast(signal, b, array(TRUE, rep(1,3)))
unlist(ADC) # text output

par(mfrow=c(1,1)) # graphical output
plot(b, signal, xlab="b-value", ylab="Signal intensity")
lines(seq(0,800,10), dwi(S0, seq(0,800,10), D), lwd=2, col=1)
lines(seq(0,800,10), dwi(ADC$S0, seq(0,800,10), ADC$D), lwd=2, col=2)
legend("topright", c("True","Estimated"), lwd=2, col=1:2)
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\keyword{models}
