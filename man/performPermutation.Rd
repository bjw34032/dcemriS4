% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{performPermutation}
\alias{performPermutation}
\title{performPermutation}
\description{
  performPermutation takes an orthogonal permutation matrix T, an array of
dimensions and a single dimensional representation of data with those
dimensions. It will return a transformed array with the transformed dimensions.
}
\usage{
performPermutation(T, real.dimensions, data, verbose=FALSE)
}
\arguments{
  \item{T}{An orthogonal matrix}
  \item{real.dimensions}{An array of single dimension, representing the length of dimensions in data}
  \item{data}{A single dimensional array representation of the data to be transformed}
  \item{verbose}{An option telling the function to be verbose}
}
\details{
  This function is mainly used by the \code{\link{reorient}} function to
transform nifti data into neuroradiological convention.
}
%\value{}
%\references{}
\seealso{\code{\link{reorient}},\code{\link{inverseReorient}}}
%\examples{}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\keyword{}
