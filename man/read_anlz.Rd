% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{readANALYZE}
\alias{readANALYZE}
\alias{read.analyze.content}
\title{readANALYZE}
\description{
  These functions read in the header information and multi-dimensional
  array from a binary file in Analyze 7.5 format.
}
\usage{
readANALYZE(fname, verbose = FALSE, warn = -1)
read.analyze.content(fname, gzipped = TRUE, verbose = FALSE, warn = -1)
}
\arguments{
  \item{fname}{Pathname of the Analyze pair of files .img and .hdr
    without the suffix.}
  \item{gzipped}{Allows the import/export of compressed (.gz) files
    (default = \code{TRUE}).}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
  \item{warn}{is a number to regulate the display of warnings (default =
    -1).  See \code{options} for more details.}
}
\details{
  The internal functions \code{readBin} and \code{rawToChar} are utilized
  in order to efficiently extract information from a binary file.  The
  types of data are limited to 1- and 2-byte integers, 4-byte floats and
  8-byte doubles.
}
\value{
  An object of class \code{anlz} is produced.
}
\references{
  ANALYZE 7.5\cr
  \url{http://www.mayo.edu/bir/PDF/ANALYZE75.pdf}
}
\seealso{
  \code{\link{readNIfTI}}
}
\examples{
avg152T1 <- readANALYZE(system.file("nifti/avg152T1.img.gz",
                                    package="dcemriS4"))
image(avg152T1)
\dontrun{
orthographic(avg152T1)
}
}
\author{Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>,
  Volker Schmid <\email{volkerschmid@users.sourceforge.net}>}
\keyword{file}
