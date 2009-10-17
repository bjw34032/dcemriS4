% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{anlz}
\Rdversion{1.1}
\alias{anlz}
\alias{is.anlz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructor for Analyze}
\description{
  Constructor for Analyze class objects.
}
\usage{
anlz(img = array(0, dim = rep(1, 4)), dim, ...)
is.anlz(x)
}
\arguments{
  \item{img}{is a multidimensional array of data.}
  \item{dim}{is the dimension of the data (default = \code{missing}).}
  \item{\dots}{allows for additional \sQuote{slots} to be specified.}
  \item{x}{is an object to be checked.}
}
%\details{}
\value{
  An object of class \code{anlz}.
}
\references{
  ANALYZE 7.5\cr
  \url{http://www.mayo.edu/bir/PDF/ANALYZE75.pdf}
}
\author{
  Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>
}
\seealso{
  \code{\linkS4class{anlz}}, \code{\link{nifti}}
}
\examples{
aim <- anlz() # default
aim
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
