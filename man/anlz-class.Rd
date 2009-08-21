% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: $

\name{anlz-class}
\Rdversion{1.1}
\docType{class}
\alias{anlz-class}
\alias{descrip<-,anlz-method}
\alias{descrip,anlz-method}
\alias{show,anlz-method}

\title{Class "anlz"}
\description{	 ~~ A concise (1-5 lines) description of what the class is.  ~~}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("anlz", data, dim, dimnames, ...)}.
	 ~~ describe objects here ~~ 
}
\section{Slots}{
	 \describe{
    \item{\code{.Data}:}{Object of class \code{"array"} ~~ }
    \item{\code{sizeof_hdr}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{db_type}:}{Object of class \code{"character"} ~~ }
    \item{\code{db_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{extents}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{session_error}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{regular}:}{Object of class \code{"character"} ~~ }
    \item{\code{hkey_un0}:}{Object of class \code{"character"} ~~ }
    \item{\code{dim_}:}{Object of class \code{"vector"} ~~ }
    \item{\code{vox_units}:}{Object of class \code{"character"} ~~ }
    \item{\code{cal_units}:}{Object of class \code{"character"} ~~ }
    \item{\code{unused1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{datatype}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{bitpix}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{dim_un0}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{pixdim}:}{Object of class \code{"vector"} ~~ }
    \item{\code{vox_offset}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{funused1}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{funused2}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{funused3}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cal_max}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{cal_min}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{compressed}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{verified}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{glmax}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{glmin}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{descrip}:}{Object of class \code{"character"} ~~ }
    \item{\code{aux_file}:}{Object of class \code{"character"} ~~ }
    \item{\code{orient}:}{Object of class \code{"character"} ~~ }
    \item{\code{originator}:}{Object of class \code{"character"} ~~ }
    \item{\code{generated}:}{Object of class \code{"character"} ~~ }
    \item{\code{scannum}:}{Object of class \code{"character"} ~~ }
    \item{\code{patient_id}:}{Object of class \code{"character"} ~~ }
    \item{\code{exp_date}:}{Object of class \code{"character"} ~~ }
    \item{\code{exp_time}:}{Object of class \code{"character"} ~~ }
    \item{\code{hist_un0}:}{Object of class \code{"character"} ~~ }
    \item{\code{views}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{vols_added}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{start_field}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{field_skip}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{omax}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{omin}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{smax}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{smin}:}{Object of class \code{"numeric"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{array}"}, from data part.
Class \code{"\linkS4class{matrix}"}, by class "array", distance 2, with explicit test and coerce.
Class \code{"\linkS4class{structure}"}, by class "array", distance 2.
Class \code{"\linkS4class{vector}"}, by class "array", distance 3, with explicit coerce.
Class \code{"\linkS4class{vector}"}, by class "array", distance 5, with explicit test and coerce.
}
\section{Methods}{
  \describe{
    \item{descrip<-}{\code{signature(x = "anlz")}: ... }
    \item{descrip}{\code{signature(object = "anlz")}: ... }
    \item{image}{\code{signature(x = "anlz")}: ... }
    \item{show}{\code{signature(object = "anlz")}: ... }
	 }
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\linkS4class{nifti}}, \code{\linkS4class{niftiExtension}}
}
\examples{
showClass("anlz")
}
\keyword{classes}
