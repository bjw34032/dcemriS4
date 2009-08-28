% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: miscellaneous.Rd 191 2009-08-25 15:12:31Z bjw34032 $

\name{Dimension Accessor Functions}
\alias{nsli}
\alias{NSLI}
\alias{ntim}
\alias{NTIM}
\title{Dimension Accessor Functions}
\description{
  Functions to extract the higher dimensions from Analyze/NIfTI data.
}
\usage{
nsli(x)
ntim(x)
NSLI(x)
NTIM(x)
}
\arguments{
  \item{x}{is a three- or four-dimensional array (e.g., read in from an
    ANALYZE/NIfTI file).}
}
\details{
  Simple calls to \code{dim} to replicate the functionality of
  \code{nrow} and \code{ncol} for higher dimensions of an array that are
  commonly required when manipulating medical imaging data.
}
\value{
  Third (slice) or fourth (time) dimension of the array.
}
% \references{}
\seealso{
  \code{\link{readNIfTI}}, \code{\link{readANALYZE}}
}
% \examples{}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>}
\keyword{misc}
