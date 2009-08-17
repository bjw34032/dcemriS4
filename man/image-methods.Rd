% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: $

\name{image-methods}
\docType{methods}
\alias{image.nifti}
\alias{image-methods}
\alias{image,ANY-method}
\alias{image,anlz-method}
\alias{image,nifti-method}
\title{Methods for Function `image' in Package `dcemri'}
\description{
  Produce images for \code{nifti} objects.
}
\usage{
\S4method{image}{nifti}(x, z=1, w=1, col=gray(0:64/64), plot.type=c("multiple","single"),
      xlab="", ylab="", axes=FALSE, ...)
\S4method{image}{anlz}(x, z=1, w=1, col=gray(0:64/64), plot.type=c("multiple","single"),
      xlab="", ylab="", axes=FALSE, ...)
}
\arguments{
 \item{x}{is an object of class \code{"nifti"}.}
 \item{z}{is the slice to be displayed (ignored when \code{plot.type =
     "multiple"}).}
 \item{w}{is the time point to be displayed (4D arrays only).}
 \item{col}{is grayscale (by default).}
 \item{plot.type}{allows the choice between all slices being displayed,
   in a matrix (left-to-right, top-to-bottom), or a single slice.}
 \item{xlab}{is set to "" since all margins are set to zero.}
 \item{ylab}{is set to "" since all margins are set to zero.}
 \item{axes}{is set to \code{FALSE} since all margins are set to zero.}
 \item{\dots}{other arguments to the \code{image} function may be
   provided here.}
}
\section{Methods}{
  \describe{
    \item{x = "ANY"}{Generic function: see \code{\link[graphics]{image}}.}
    \item{x = "nifti"}{Produce images for \code{x}.}
    \item{x = "anlz"}{Produce images for \code{x}.}
  }
}
\details{
}
\value{
}
\references{
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\seealso{
  \code{\link{orthographic-methods}}, \code{\link{overlay-methods}}
}
\examples{
}
\keyword{methods}
