\name{Convert NIfTI Codes}
\alias{convert.datatype}
\alias{convert.intent}
\alias{convert.form}
\alias{convert.units}
\alias{convert.slice}
\title{Convert NIfTI Codes}
\description{
  Codes that appear in the NIfTI-1 header are mapped to meaningful
  chartacter strings.
}
\usage{
convert.datatype(datatype)
convert.intent(intent.code)
convert.form(form.code)
convert.units(units)
convert.slice(slice.code)
}
\arguments{
  \item{datatype}{defines data type.}
  \item{intent.code}{is the NIfTI intent code.}
  \item{form.code}{is the x,y,z coordinate system.}
  \item{units}{is the units of pixdim[1..4].}
  \item{slice.code}{is the slice timing order.}
}
\details{
  \code{switch} statements are used to map a numeric code to the
  appropriate string.
}
\value{
  A character string.
}
\references{
  Neuroimaging Informatics Technology Initiative (NIfTI)
  \url{http://nifti.nimh.nih.gov/}
}
%\seealso{}
%\examples{}
\author{B. Whitcher <\email{bjw34032@users.sourceforge.net}>}
%\keyword{}
