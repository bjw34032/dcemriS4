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
\alias{getLastCallWithName}
\title{Facilitate the Creation and Modification of Audit Trails}
\description{
  Facilitate the creation and modification of audit trails for NIfTI
  class objects.
}
\usage{
dcemri.info(type)
enableAuditTrail()
newAuditTrail()
niftiExtensionToAuditTrail(nim, filename = NULL, call = NULL)
niftiAuditTrailToExtension(nim, filename = filename, call = call)
niftiAuditTrailSystemNode(type="system-info", filename = NULL, call = NULL)
niftiAuditTrailSystemNodeEvent(trail, type = NULL, call = NULL,
                               filename = NULL, comment=NULL)
niftiAuditTrailCreated(history = NULL, call = NULL, filename = NULL)
niftiAuditTrailEvent(trail, type = NULL, call = NULL, comment = NULL)
getLastCallWithName(functionName)
}
\arguments{
  \item{nim}{is an object of class \code{niftiAuditTrail} or can be converted to such.}
  \item{filename}{The filename associated with the nifti object.}
  \item{call}{A \code{call}, function name in the call-stack or a string.}
  \item{type}{An identifier to add some meaning to the event}
  \item{trail}{The \code{XMLAbstractNode} representing the audit trail or the \code{niftiAuditTrail} object with a trail that will be amended.}
  \item{comment}{Some textual comment}
  \item{history}{An \code{XMLAbstractNode} to take as history to add to the \code{trail}}
  \item{functionName}{The name of a function on the call stack}
}
\details{
  The function \code{dcemri.info} is used to find the \code{ecode} or the XML
\code{namespace} relevant to the audit trail.

  The function \code{enableAuditTrail} is called as the package loads to switch
on the audit trail functionality and declare the class \code{niftiAuditTrail}.
Should you wish to prevent the audit trail functionality loading, you should set
the option \code{NIfTI.audit.trail} to FALSE. Manually calling the function
\code{enableAuditTrail} will override this setting and switch on the
functionality.

  The function \code{newAuditTrail} returns an \code{XMLAbstractNode}
representing the root node of an audit trail. This is mostly intended as an
internal function.

  The function \code{niftiExtensionToAuditTrail} takes an object representing a
NIfTI object, casts it as a \code{niftiAuditTrail} and checks if there is an
extension (a \code{niftiExtensionSection}) with \code{ecode} equal to
\code{dcemri.info("ecode")} i.e. has a extension with data representing a
serialized audit trail. The function will then strip
the object of this extension parsing the serialized \code{edata} into an audit
trail and adding a \code{read} event to the trail.

  The function \code{niftiAuditTrailToExtension} takes a \code{niftiAuditTrail}
and returns a \code{niftiExtensionSection} with \code{edata} containing the
serialized form of the audit trail after adding a \code{saved} event to the
trail.

  The function \code{niftiAuditTrailSystemNodeEvent} adds an element with name
 equal to \code{type} to the \code{trail}. It uses the
\code{niftiAuditTrailSystemNode} function to create the node.
  
  The function \code{niftiAuditTrailSystemNode} is an internal function creating
an \code{XMLAbstractNode} element with name \code{type} and attributes giving
information about the R system and library. The \code{filename} and \code{call}
will also be added as attributes if available.

  The function \code{niftiAuditTrailEvent} adds an element with name
\code{event} to the \code{trail}. The arguments \code{type}, \code{filename},
\code{call} are added as attributes and the \code{comment} is the text value of
the element.

  The function \code{niftiAuditTrailCreated} will create a new audit trail
containing a system node element \code{created} with the child \code{history}
with the contents \code{history}. If the last element of the \code{history}
given is an \code{event} with \code{type="processing"}, then this node will be
removed from the \code{history} and its \code{call} attribute will be used as
the value of the \code{call} attribute on the \code{created} node.

  The function \code{getLastCallWithName} will search the call stack for a call
of the function \code{functionName}, returning last call to that function if
possible. It will default to the call of the function which called the function
which called getLastCallWithName if there was no such call. (And if there was no
such call it will return the call of itself.)
}
\note{
These functions are mostly intended to be used internally in order to document
the changes that occur to NIfTI objects due to functions that are audit-trail
aware. However, as the precise manner in which these functions are used is not
documented anywhere else, I shall proceed to describe which functions are
audit-trail aware and how they interact with the audit trail.

\code{as.nifti} and its S4 alias \code{as(nim, "nifti")} will always produce
\code{niftiAuditTrail} objects if the functionality is on. The function
\code{niftiAuditTrailCreated} will be used and if an exemplar object
is provided e.g. \code{as.nifti(array, niftiExemplar)} then the trail of the
exemplar will be used as the \code{history}. 

\code{readNIfTI} and \code{writeNIfTI} also always produce
\code{niftiAuditTrail} objects if the functionality is on. The functions
\code{niftiExtensionToAuditTrail} and \code{niftiAuditTrailToExtension} are used
internally by these functions to facilitate this behaviour.

}
%\value{ }
%\references{}
\author{
  Andrew Thornton <\email{zeripath@users.sourceforge.net}>
}
%\seealso{}
\examples{
## A good example of the use of these functions is shown by this wrapper
## function which takes a function fun(nim, ...) returning lists of arrays which
## are nifti-ized using as(...)
wrapper <- function(functionToWrap, nameOfCallingFunction, nim, ...) {
  if (!is(nim, "nifti")) 
    nim <- as(nim, "nifti")
  
  if (is(nim, "niftiAuditTrail")) {
    ## This will force as(...) to set the call which created the results to the
    ## calling function's call rather than as(result, nifti) as it would
    ## otherwise do
    nim@trail <- niftiAuditTrailEvent(nim@trail, "processing", nameOfCallingFunction)
  }
  
  result <- functionToWrap(nim, ...)
  as(result, "nifti") <- nim
  return(result)
}

## An example of how wrapper is used follows:
functionToWrap <- function(ignored,x,y) {
  return (array(1, dim=c(x,y)))
}

## The niftiized form
niftiizedForm <- function(nim,...) {
  return(wrapper(functionToWrap, "niftiizedForm", nim, ...))
}

## compare the trails
if (getOption("NIfTI.audit.trail")) {
  print((as(functionToWrap(nifti(), 4, 4), "nifti") <- nifti())@trail)
  print(niftiizedForm(nifti(), 4, 4)@trail)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
