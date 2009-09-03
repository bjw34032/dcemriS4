% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{descrip-methods}
\docType{methods}
\alias{descrip-methods}
\alias{descrip,anlz-method}
\alias{descrip,nifti-method}
\alias{descrip<-,anlz-method}
\alias{descrip<-,nifti-method}
\alias{descrip}
\alias{descrip<-}
\title{Extract or Replace NIfTI/Analyze Description}
\description{
  Methods that act on the \dQuote{description} character string in the
  NIfTI or Analyze header.
}
\usage{
\S4method{descrip}{nifti}(x)
\S4method{descrip}{anlz}(x)
}
\arguments{
  \item{x}{is an object of class \code{nifti} or \code{anlz}.}
}
\section{Methods}{
  \describe{
  \item{object = "anlz"}{Extract or replace Analyze description.}
  \item{object = "nifti"}{Extract or replace NIfTI description.}
  }
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\examples{
mniRL <- readNIfTI(system.file("nifti/avg152T1_RL_nifti.nii.gz",
                   package="dcemriS4"))
descrip(mniRL)
descrip(mniRL) <- paste(descrip(mniRL), version$version.string, sep="; ")
descrip(mniRL)
}
\keyword{methods}
