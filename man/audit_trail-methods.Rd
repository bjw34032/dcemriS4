% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: $

\name{audit.trail-methods}
\docType{methods}
\alias{audit.trail-methods}
\alias{audit.trail,nifti-method}
\alias{audit.trail}
\alias{audit.trail<-,nifti-method}
\alias{audit.trail<-}
\title{Extract or Replace NIfTI Audit Trail}
\description{
  Operators that act on the audit trail (XML) in the NIfTI header.
}
\usage{
\S4method{audit.trail}{nifti}(x)
}
\arguments{
  \item{x}{is an object of class \code{nifti}.}
}
\section{Methods}{
  \describe{
    \item{object = "nifti"}{Extract or replace NIfTI audit trail.}
  }
}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\examples{}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
