% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{reorient}
\alias{reorient}
\alias{inverseReorient}
\title{reorient}
\description{
  Transforms in the NIfTI header are parsed and normalized versions of
  these transforms are applied.
}
\usage{
reorient(nim, data, verbose=FALSE, invert=FALSE)
inverseReorient(nim, verbose=FALSE)
}
\arguments{
  \item{nim}{is an object of class \code{nifti}.}
  \item{data}{is an array associated with \code{nim}.}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
    \code{FALSE}).}
  \item{invert}{stores the inverse transform.}
}
\details{
  This function utilizes the \code{performPermutation} function
  internally.
}
%\value{}
%\references{}
\seealso{
  \code{\link{performPermutation}}
}
%\examples{}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\keyword{}
