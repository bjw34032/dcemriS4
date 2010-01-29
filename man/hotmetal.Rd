% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-14 17:34:49 (bjw34032)>
% $Id$

\name{hotmetal}
\alias{hotmetal}
\title{Hot Metal Color Table}
\description{
  The hotmetal color table patterned after the one used in Matlab. 
}
\usage{
hotmetal(n=64)
}
\arguments{
  \item{n}{is the number of color levels (default = 64).}
}
\details{
  Based on the \code{\link[fields]{tim.colors}} function in the
  \pkg{fields} package.  The \code{hotmetal} function has been modified
  to break any dependence on code in the \pkg{fields} package.  Spline
  interpolation (\code{interpSpline}) is used when the number of
  requested colors is not the default.
}
\value{
  A vector of character strings giving the colors in hexadecimal
  format.
}
\seealso{
  \code{\link{topo.colors}}, \code{\link{terrain.colors}}
}
\examples{
hotmetal(10) 
image(outer(1:20,1:20,"+"), col=hotmetal(75), main="hotmetal")
}
\keyword{aplot}
