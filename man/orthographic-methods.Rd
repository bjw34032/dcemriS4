% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{orthographic-methods}
\docType{methods}
\alias{orthographic}
\alias{orthographic-methods}
\alias{orthographic,nifti-method}
\alias{orthographic,anlz-method}
\alias{orthographic,array-method}
\alias{orthographic.nifti}
\title{Methods for Function `orthographic' in Package `dcemri'}
\description{
  Produce orthographic display for \code{nifti} objects.
}
\usage{
\S4method{orthographic}{nifti}(x, xyz = NULL, crosshairs = TRUE,
                       col.crosshairs = "red", w = 1, zlim = NULL,
                       col = gray(0:64/64), xlab = "", ylab = "",
                       axes = FALSE, ...)
}
\arguments{
  \item{x}{is an object of class \code{"nifti"} or similar.}
  \item{xyz}{is the coordinate for the center of the crosshairs.}
  \item{crosshairs}{is a logical value for the presence of crosshairs in
    all three orthogonal planes (default = TRUE).}
  \item{col.crosshairs}{is the color of the crosshairs (default = red).}
  \item{w}{is the time point to be displayed (4D arrays only).}
  \item{col}{is grayscale (by default).}
  \item{xlab}{is set to "" since all margins are set to zero.}
  \item{ylab}{is set to "" since all margins are set to zero.}
  \item{axes}{is set to \code{FALSE} since all margins are set to zero.}
  \item{\dots}{other arguments to the \code{image} function may be
    provided here.}
}
\section{Methods}{
  \describe{
    \item{x = "anlz"}{Produce orthographic display for \code{x}.}
    \item{x = "array"}{Produce orthographic display for \code{x}.}
    \item{x = "nifti"}{Produce orthographic display for \code{x}.}
  }
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\seealso{
  \code{\link{image-methods}}, \code{\link{overlay-methods}}
}
%\examples{}
\\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
