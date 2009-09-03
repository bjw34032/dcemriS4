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
  Methods that act on the \dQuote{auxiliary file} character string in
  the NIfTI or Analyze header.
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
\examples{
mniRL <- readNIfTI(system.file("nifti/avg152T1_RL_nifti.nii.gz",
                   package="dcemriS4"))
aux.file(mniRL)
aux.file(mniRL) <- "avg152T1_RL_nifti"
aux.file(mniRL)
}
\keyword{methods}
