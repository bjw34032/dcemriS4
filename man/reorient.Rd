% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{reorient}
\alias{reorient}
\alias{inverseReorient}
\title{reorient}
\description{
  Reorient takes nim and data; a nifti header and a single dimensional array representing the nifti data; and some other options and will attempt to return a nifti object with the data reoriented into neuroradiological convention. It will parse the transforms in the nifti header and apply a normalized version of these transforms. It will not update the transforms in the header.
}
\usage{
reorient(nim, data, verbose=FALSE, invert=FALSE)
inverseReorient(nim, verbose=FALSE)
}
\arguments{
  \item{nim}{A nifti header}
  \item{data}{The data array associated with nim}
  \item{verbose}{Be verbose}
  \item{invert}{Invert our reorientation for saving}
}
\details{This function will use the performPermutation function internally.}
%\value{}
%\references{}
\seealso{\code{\link{performPermutation}}}
%\examples{}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\keyword{}
