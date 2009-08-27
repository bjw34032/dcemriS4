% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{aux.file-methods}
\docType{methods}
\alias{aux.file-methods}
\alias{aux.file,anlz-method}
\alias{aux.file,nifti-method}
\alias{aux.file}
\alias{aux.file<-,nifti-method}
\alias{aux.file<-,anlz-method}
\alias{aux.file<-}
\title{Extract or Replace NIfTI/Analyze Auxiliary File}
\description{
  Operators that act on character strings in the NIfTI or Analyze
  header.
}
\usage{
\S4method{aux.file}{nifti}(x)
\S4method{aux.file}{anlz}(x)
}
\arguments{
  \item{x}{is an object of class \code{nifti} or \code{anlz}.}
}
\section{Methods}{
  \describe{
    \item{object = "anlz"}{Extract or replace Analyze auxiliary file.}
    \item{object = "nifti"}{Extract or replace NIfTI auxiliary file.}
  }
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\examples{}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
