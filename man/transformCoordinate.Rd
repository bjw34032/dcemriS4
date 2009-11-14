% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{translateCoordinate}
\alias{translateCoordinate}
\title{translateCoordinate}
\description{
  ...
}
\usage{
translateCoordinate(i, nim, verbose=FALSE)
}
\arguments{
  \item{i}{An index vector in nim}
  \item{nim}{A NIfTI object}
  \item{verbose}{Be verbose}
}
\details{
  translateCoordinate takes a NIfTI object and an index vector in to the voxel space of the object and translates that voxel index into the continuous co-ordinate space defined by the objects qform and sform.
}
%\value{
%  ...
%}
%\references{}
%\seealso{}
%\examples{}
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\keyword{}
