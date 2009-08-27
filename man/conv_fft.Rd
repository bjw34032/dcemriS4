% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-08-05 15:47:57 (bjw34032)>
% $Id: conv_fft.Rd 191 2009-08-25 15:12:31Z bjw34032 $

\name{conv.fft}
\alias{conv.fft}
\title{Convolution of 3D Arrays using the Fourier Transform}
\description{
  Convolve a 
}
\usage{
conv.fft(A, B, C, FFTA=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{is a three-dimensional array (\dQuote{the template}).}
  \item{B}{is a three-dimensional array (\dQuote{the target}).}
  \item{C}{is a vector of length three (the center of \dQuote{the template}).}
  \item{FFTA}{is the three-dimensional Fourier transform of \code{A},
    this may save time when looping over multiple \dQuote{targets}.}
}
\details{
  The arrays \eqn{A} and \eqn{B} are transformed into the Fourier domain
  and multiplied together (equivalent to a convolution in the image
  domain across all spatial locations simultaneously).  
}
\value{
  A three-dimensional array, the same dimension as the input arrays,
  that is the convolution of the \dQuote{target} to the
  \dQuote{template} at all spatial locations.
}
\references{
  Briggs, W.L. and Henson, V.E. (1995) \emph{The DFT: An Owner's Manual
    for the Discrete Fourier Transform}, SIAM: Philadelphia.
}
\author{Brandon Whitcher}
\seealso{
  \code{\link{fft}}, \code{\link{ftm}}, \code{\link{shift3D}}
}
\examples{
cube <- array(0, c(20,20,1))
cube[9:12,9:12,1] <- 1
tkernel <- array(0, c(20,20,1))
tkernel[,,1] <- c(.5, 1, .5, rep(0,20-3)) %o% c(.5, 1, .5, rep(0,20-3))
tcenter <- find.center(ifelse(tkernel > 0, TRUE, FALSE))
out <- conv.fft(tkernel, cube, tcenter)
out[8:13,8:13,1]  ## text output

par(mfrow=c(2,2))  ## graphic output
image(drop(tkernel), col=tim.colors(), main="Template")
image(drop(cube), col=tim.colors(), main="Target")
image(drop(out), col=tim.colors(), main="Output")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
