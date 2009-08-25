% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: nifti-class.Rd 191 2009-08-25 15:12:31Z bjw34032 $

\name{niftiExtensionToAuditTrail}
\Rdversion{1.1}
\alias{enableAuditTrail}
\alias{newAuditTrail}
\alias{niftiExtensionToAuditTrail}
\alias{niftiAuditTrailToExtension}
\alias{niftiAuditTrailSystemNode}
\alias{niftiAuditTrailSystemNodeEvent}
\alias{niftiAuditTrailCreated}
\alias{niftiAuditTrailEvent}
\title{Facilitate the Creation and Modification of Audit Trails}
\description{
  Facilitate the creation and modification of audit trails for NIfTI
  class objects.
}
\usage{
enableAuditTrail()
newAuditTrail()
niftiExtensionToAuditTrail(nim, filename = NULL, call = NULL)
niftiAuditTrailToExtension(nim, filename = filename, call = call)
niftiAuditTrailSystemNode(type="system-info", filename = NULL, call = NULL)
niftiAuditTrailSystemNodeEvent(trail, type = NULL, call = NULL,
                               filename = NULL, comment=NULL)
niftiAuditTrailCreated(history = NULL, call = NULL, filename = NULL)
niftiAuditTrailEvent(trail, type = NULL, call = NULL, comment = NULL)
}
\arguments{
  \item{nim}{}
  \item{filename}{}
  \item{call}{}
  \item{type}{}
  \item{trail}{}
  \item{comment}{}
  \item{history}{}
}
\details{
  ...
}
\value{
  ...
}
%\references{}
\author{
  Andrew Thornton <\email{zeripath@users.sourceforge.net}>
}
%\seealso{}
%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
