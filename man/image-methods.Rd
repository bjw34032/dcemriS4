% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-03 16:42:44 (bjw34032)>

\name{image-methods}
\docType{methods}
\alias{image-methods}
\alias{image,ANY,ANY-method}
\alias{image,nifti,nifti-method}
\title{Methods for Function `image' in Package `dcemri'}
\description{
  Produce images for \code{"nifti"} objects.
}
\usage{
\S4method{image}{nifti,nifti}(x, z=1, w=1, col=gray(0:128/128), plot.type=c("multiple","single"),
     xlab="", ylab="", axes=FALSE, ...)
}
\arguments{
 \item{x}{is an object of class \code{"nifti"}.}
 \item{z}{is the slice to be displayed (ignored when \code{plot.type =
     "multiple"}).}
 \item{w}{is the time point to be displayed (for 4D arrays only).}
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
  }
}
\keyword{methods}
