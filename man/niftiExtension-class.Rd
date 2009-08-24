% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: $

\name{niftiExtension-class}
\Rdversion{1.1}
\docType{class}
\alias{niftiExtension-class}
\title{Class "niftiExtension"}
\description{
  ~~ A concise (1-5 lines) description of what the class is.  ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("niftiExtension", data, dim, dimnames, ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"array"} ~~ }
    \item{\code{extensions}:}{Object of class \code{"list"} ~~ }
    \item{\code{sizeof_hdr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{data_type}:}{Object of class \code{"character"} ~~ }
    \item{\code{db_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{extents}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{session_error}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{regular}:}{Object of class \code{"character"} ~~ }
    \item{\code{dim_info}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dim_}:}{Object of class \code{"vector"} ~~ }
    \item{\code{intent_p1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{intent_p2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{intent_p3}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{intent_code}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{datatype}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bitpix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{slice_start}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pixdim}:}{Object of class \code{"vector"} ~~ }
    \item{\code{vox_offset}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{scl_slope}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{scl_inter}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{slice_end}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{slice_code}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{xyzt_units}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cal_max}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cal_min}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{slice_duration}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{toffset}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{glmax}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{glmin}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{descrip}:}{Object of class \code{"character"} ~~ }
    \item{\code{aux_file}:}{Object of class \code{"character"} ~~ }
    \item{\code{qform_code}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sform_code}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{quatern_b}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{quatern_c}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{quatern_d}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qoffset_x}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qoffset_y}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{qoffset_z}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{srow_x}:}{Object of class \code{"vector"} ~~ }
    \item{\code{srow_y}:}{Object of class \code{"vector"} ~~ }
    \item{\code{srow_z}:}{Object of class \code{"vector"} ~~ }
    \item{\code{intent_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{magic}:}{Object of class \code{"character"} ~~ }
    \item{\code{extender}:}{Object of class \code{"vector"} ~~ }
    \item{\code{reoriented}:}{Object of class \code{"logical"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{nifti}"}, directly.
Class \code{"\linkS4class{array}"}, by class "nifti", distance 2.
Class \code{"\linkS4class{matrix}"}, by class "nifti", distance 3, with explicit test and coerce.
Class \code{"\linkS4class{structure}"}, by class "nifti", distance 3.
Class \code{"\linkS4class{vector}"}, by class "nifti", distance 4, with explicit coerce.
Class \code{"\linkS4class{vector}"}, by class "nifti", distance 6, with explicit test and coerce.
}
\section{Methods}{
  No methods defined with class "niftiExtension" in the signature.
}
\references{ ~put references to the literature/web site here ~ }
\author{Andrew Thornton <\email{zeripath@users.sourceforge.net}>}
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\linkS4class{nifti}}, \code{\linkS4class{niftiAuditTrail}}
}
\examples{
showClass("niftiExtension")
}
\keyword{classes}
