% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% 
% $Id: extract.Rd 191 2009-08-25 15:12:31Z bjw34032 $

\name{Seed Growing for a 4D Array}
\alias{extract.aif}
\title{Seed Growing for a 4D Array}
\description{
  Seed growing algorithm to find voxels in a three-dimensional array
  according to their correlation to a seed voxel.  The correlation is
  measured according to the fourth dimension of the array.
}
\usage{
extract.aif(img, x, y, z, thresh = 0.9)
}
\arguments{
  \item{img}{is the four-dimensional array of medical imaging data.}
  \item{x,y,z}{are the coordinates of the seed voxel.}
  \item{thresh}{is the minimum correlation for inclusion in the region.}
}
\details{
  Correlation coefficients are computed for every voxel in the input
  array.  A recursive algorithm is then used to grow the region of
  interest (\acronym{ROI}) from the seed voxel in three dimensions.  All
  adjacent voxels, where the correlation exceeds the threshold, are
  included.
}
\value{
  \item{coord}{is a matrix of the three-dimesional coordinates
    \eqn{(x,y,z)} for all voxels found by the algorithm.}
  \item{conc}{is a matrix whose rows correspond to the voxels found by
    the algorithm and whose columns are the fourth dimension from the
    input array (e.g., contrast agent concentration time curve).}
  \item{mask}{is an array of boolean values, where only voxels included
    by the algorithm are given a value greater than zero.}
  \item{cor}{is an array that mimics the \code{mask}, but contains the
    estimated correlation coefficients for all voxels included by the
    algorithm.}
}
%\references{}
%\seealso{}
%\examples{}
\author{Volker Schmid <\email{volker.schmid@users.sourceforge.net}>}
\keyword{misc}
