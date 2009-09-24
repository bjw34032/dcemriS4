% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{overlay-methods}
\docType{methods}
\alias{overlay}
\alias{overlay-methods}
\alias{overlay,nifti,nifti-method}
\alias{overlay,nifti,anlz-method}
\alias{overlay,nifti,array-method}
\alias{overlay,anlz,nifti-method}
\alias{overlay,anlz,anlz-method}
\alias{overlay,anlz,array-method}
\alias{overlay,array,nifti-method}
\alias{overlay,array,anlz-method}
\alias{overlay,array,array-method}
\alias{overlay.nifti}
\title{Methods for Function overlay}
\description{Methods for function \code{overlay}}
\usage{
\S4method{overlay}{nifti,nifti}(x, y, z=1, w=1, col.x=gray(0:64/64),
      col.y=hotmetal(), zlim.x=NULL, zlim.y=NULL,
      plot.type=c("multiple","single"), xlab="", ylab="",
      axes=FALSE, ...)
\S4method{overlay}{anlz,anlz}(x, y, z=1, w=1, col.x=gray(0:64/64),
      col.y=hotmetal(), zlim.x=NULL, zlim.y=NULL,
      plot.type=c("multiple","single"), xlab="", ylab="",
      axes=FALSE, ...)
}
\arguments{
  \item{x,y}{is an object of class \code{"nifti"} or \code{"anlz"}.}
  \item{z}{is the slice to be displayed (ignored when \code{plot.type =
      "multiple"}).}
  \item{w}{is the time point to be displayed (4D arrays only).}
  \item{col.x}{is grayscale (by default).}
  \item{col.y}{is hotmetal (by default).}
  \item{zlim.x,zlim.y}{are set to \code{NULL} (by default) and taken
    from the header information.}
  \item{plot.type}{allows the choice between all slices being displayed,
    in a matrix (left-to-right, top-to-bottom), or a single slice.}
  \item{xlab}{is set to "" since all margins are set to zero.}
  \item{ylab}{is set to "" since all margins are set to zero.}
  \item{axes}{is set to \code{FALSE} since all margins are set to zero.}
  \item{\dots}{other arguments to the \code{image} function may be
    provided here.}
}
\details{
  
}
\section{Methods}{
\describe{
  \item{x = "nifti", y = "nifti"}{Produce overlay of \code{y} on \code{x}.}
  \item{x = "anlz", y = "anlz"}{Produce overlay of \code{y} on \code{x}.}
 }
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\seealso{
  \code{\link{image-methods}}, \code{\link{overlay-methods}}
}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
