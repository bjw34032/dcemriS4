% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{niftiExtensionSection-class}
\Rdversion{1.1}
\docType{class}
\alias{niftiExtensionSection-class}
\title{Class "niftiExtensionSection"}
\description{
  ...
}
\section{Objects from the Class}{
  Objects can be created by calls of the form
  \code{new("niftiExtensionSection", data, dim, dimnames, ...)}. 
}
\section{Slots}{
  \describe{
    \item{\code{esize}:}{Object of class \code{"numeric"}; the number of
      bytes that form the extended header data.}
    \item{\code{ecode}:}{Object of class \code{"numeric"}; a
      non-negative integer that indicates the format of the extended
      header data that follows (default = 1002).} 
    \item{\code{edata}:}{Object of class \code{"character"}; Note that
      the other contents of the extended header data section are totally
      unspecified by the NIfTI-1 standard.}
  }
}
\section{Methods}{
  No methods defined with class "niftiExtensionSection" in the signature.
}
\references{
  NIfTI-1\cr
  \url{http://nifti.nimh.nih.gov/}
}
\author{
  Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>,
  Andrew Thornton <\email{zeripath@users.sourcefore.net}>
}
\seealso{
  \code{\linkS4class{niftiExtension}}, \code{\linkS4class{nifti}}
}
\examples{
showClass("niftiExtensionSection")
}
\keyword{classes}
