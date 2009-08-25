% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{Convert ANALYZE Codes}
\alias{convert.datatype.anlz}
\alias{convert.orient.anlz}
\title{Convert NIfTI Codes}
\description{
  Codes that appear in the ANALYZE header are mapped to meaningful
  chartacter strings.
}
\usage{
convert.datatype.anlz(datatype)
convert.orient.anlz(orientation)
}
\arguments{
  \item{datatype}{defines data type.}
  \item{orientation}{defines the orientation.}
}
\details{
  \code{switch} statements are used to map a numeric code to the
  appropriate string.
}
\value{
  A character string.
}
\references{
  ANALYZE 7.5 \url{}
}
%\seealso{}
%\examples{}
\author{B. Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
