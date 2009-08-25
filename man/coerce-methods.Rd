% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{coerce-methods}
\docType{methods}
%\alias{coerce-methods}
\alias{coerce,array,nifti-method}
\alias{coerce,list,nifti-method}
\alias{coerce<-,array,nifti-method}
\alias{coerce<-,list,nifti-method}
\title{Force an Object to Belong to the NIfTI Class}
\description{
 ~~ Methods for function \code{coerce} in Package `methods' ~~
}
\usage{
as(object, Class)
as(object, Class) <- value
}
\arguments{
  \item{object}{}
  \item{value}{}
}
\section{Methods}{
\describe{

\item{from = "array", to = "nifti"}{An object of class \code{array} is
  coerced into a NIfTI object.}

\item{from = "list", to = "nifti"}{All objects of class \code{array} in
  the list are coerced into NIfTI objects.  All other objects are left
  alone.  The original list structure is retained.}
}}
\references{
}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
\seealso{
}
\examples{
}
\keyword{methods}
%\keyword{ ~~ other possible keyword(s)}
