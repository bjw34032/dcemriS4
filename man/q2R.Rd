% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{quaternion2rotation}
\alias{quaternion2rotation}
\title{Convert Quaternion into a Rotation Matrix}
\description{
  The rotation matrix \eqn{R} is calculated from the quaternion
  parameters.
}
\usage{
quaternion2rotation(b, c, d)
}
\arguments{
  \item{b}{is the quaternion \eqn{b} parameter.}
  \item{c}{is the quaternion \eqn{c} parameter.}
  \item{d}{is the quaternion \eqn{d} parameter.}
}
\details{
  The quaternion representation is chosen for its compactness in
  representing rotations.  The orientation of the (x,y,z) axes relative
  to the (i,j,k) axes in 3D space is specified using a unit quaternion
  [a,b,c,d], where a*a+b*b+c*c+d*d=1.  The (b,c,d) values are all that
  is needed, since we require that a = sqrt(1.0-(b*b+c*c+d*d)) be
  nonnegative.  The (b,c,d) values are stored in the
  (quatern_b,quatern_c,quatern_d) fields.
}
\value{
  The (proper) 3x3 rotation matrix.
}
\references{
  NIfTI-1\cr
  \url{http://nifti.nimh.nih.gov/}
}
\examples{
R <- quaternion2rotation(1, 0, 0)
## This R matrix is represented by quaternion [a,b,c,d] = [0,1,0,0]
## (which encodes a 180 degree rotation about the x-axis).
R
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
