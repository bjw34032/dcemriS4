% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% Time-stamp: <2009-07-04 12:39:48 (bjw34032)>

\name{find.center}
\alias{find.center}
\title{Find the Center of a Binary Mask}
\description{
  The center of a binary mask is determined.
}
\usage{
find.center(M)
}
\arguments{
  \item{M}{is a binary mask (multidimensional array of logical values).}
}
\details{
  This method most likely only works with convex three-dimensional
  shapes (e.g., a hyper-rectangle).  Further testing is required to
  know the limits of the current implementation.
}
\value{
  A vector of values the same length as the input array.
}
%\references{}
\author{Brandon Whitcher}
\seealso{
  \code{\link{ftm}}
}
\examples{
M <- array(FALSE, rep(10,3))
M[6:10,6:10,6:10] <- TRUE
Mc <- find.center(M)
print(Mc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
