% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{as.nifti}
\alias{as.nifti}
\title{as.nifti}
\description{
  Internal function that converts multidimensional arrays to NIfTI class
  objects.
}
\usage{
as.nifti(from, value = NULL, verbose = FALSE)
}
\arguments{
  \item{from}{is the object to be converted.}
  \item{value}{is the \code{nifti} class object to use as a template for
    various NIfTI header information.}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
}
%\details{}
\value{
  An object of class \code{nifti}.
}
%\references{}
%\seealso{}
%\examples{}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\keyword{}
