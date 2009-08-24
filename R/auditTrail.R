##
##
## Copyright (c) 2009, Brandon Whitcher and Volker Schmid
## All rights reserved.
## 
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
## 
##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer. 
##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.
##     * The names of the authors may not be used to endorse or promote
##       products derived from this software without specific prior
##       written permission.
## 
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
## 
## $Id: $
##

start.nifti.audit.trail.functionality <- function() {
  if (require("XML")) {
    if (!isClass("niftiAuditTrail")) {
      options("NIfTI.audit.trail"=TRUE)
      setClass("niftiAuditTrail",
	  representation(trail="XMLNode"),
	  prototype(trail=new.audit.trail()),
	  contains="niftiExtension")
    }
  }
}

dcemri.ecode <- 1002
dcemri.namespace <- "http://www.dcemri.org/namespaces/audit-trail/1.0"

new.audit.trail <- function() {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    trail <- xmlNode("audit-trail", attrs=list(xmlns=dcemri.namespace),
                     namespace="")
    return(trail)
  }
}

nifti.extension.to.audit.trail <- function(nim, filename=NULL, call=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    if (!is(nim, "niftiAuditTrail"))
      nim <- as(nim, "niftiAuditTrail")
    ## We enforce that there is a single extension with ecode == dcemri.ecode
    ecodes <- lapply(nim@extensions, function(x) x@ecode)
    oei <- which(ecodes == dcemri.ecode)
    
    if (length(oei) == 0) {
      nim@trail <- nifti.audit.trail.created(filename=filename, call=call)
    } else {
      ## One Trail
      if (length(oei) > 1) {
	warning("Found more than one extension with ecode ==", dcemri.ecode,
                ", Appending to last trail only")
	oei <- oei[length(oei)]
      }
      oe <- nim@extensions[[oei]]@edata
      nim@extensions[[oei]] <- NULL
      nim@trail <- nifti.audit.trail.system.node.event(xmlRoot(xmlTreeParse(oe, asText=TRUE)), type="read", filename=filename, call=call)

    }
    return(nim)
  }
}

nifti.audit.trail.to.extension <- function(nim, filename=filename,
                                              call=call) {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    sec <- new("niftiExtensionSection")
    sec@ecode <- dcemri.ecode
    nim@trail <- nifti.audit.trail.system.node.event(nim@trail, "saved",
                                              filename=filename, call=call)
    ## Serialize the XML to sec@edata
    ## DIRTY DIRTY DIRTY you should wash your eyes out after reading this.
    useFancyQuotes <- getOption("useFancyQuotes")
    options("useFancyQuotes"=FALSE)
    sec@edata <- toString.XMLNode(nim@trail)
    options("useFancyQuotes"=useFancyQuotes)

    ## Fix the esize to be congruent to 0 mod 16
    sec@esize <- nchar(sec@edata, type="bytes") + 8
    sec@esize <- (-sec@esize %% 16) + sec@esize
    return(sec)
  }
}

nifti.audit.trail.system.node <- function(type="system-info", filename=NULL,
                                      call=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    if (is(call, "call"))
      call <- as.character(as.expression(call))
    currentDateTime <- format(Sys.time(), "%a %b %d %X %Y %Z")
    system <- xmlNode(type, attrs=c("filename"=filename, "call"=call),
                      namespace="")
    system <- addAttributes(system,
                            "r-version"=version$version.string,
                            "date"=currentDateTime,
                            "user"=Sys.getenv("LOGNAME"),
                            "dcemri-version"=packageDescription("dcemri")$Version)
    return(system)
  }
}

nifti.audit.trail.created <- function(history=NULL, call=NULL,
                                      filename=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    if (is(history, "niftiAuditTrail")) {
      return(nifti.audit.trail.created(history@trail, call, filename))
    } else {
      require("XML")
      trail <- new.audit.trail()
      created <- nifti.audit.trail.system.node("created", "filename"=filename,
	  "call"=call)
      if (!is.null(history)) {
	historyNode <- xmlNode("history")
	lapply(xmlChildren(history),
	    function(x) { historyNode <<- addChildren(historyNode, x) })
	created <- addChildren(created, historyNode)
      }
      trail <- addChildren(trail, created) 
      return(trail)
    }
  }
}

nifti.audit.trail.event <- function(trail, type=NULL, call=NULL,
                                    comment=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    if (is(trail,"niftiAuditTrail")) {
      return(nifti.audit.trail.event(trail@trail, type, call, comment))
    }
    require("XML")
    if (is(call,"call"))
      call <- as.character(as.expression(call))
    eventNode <- xmlNode("event", attrs=c("type"=type, "call"=call))
    if (!is.null(comment))
      eventNode <- addChildren(eventNode, xmlTextNode(comment))
    trail <- addChildren(trail, eventNode)
    return(trail)
  }
}

nifti.audit.trail.system.node.event <- function(trail, type=NULL, call=NULL,
                                            filename=NULL, comment=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    if (is(trail,"niftiAuditTrail")) {
      return(nifti.audit.trail.system.node.event(trail@trail, type, call, filename, comment))
    }
    require("XML")
    eventNode <- nifti.audit.trail.system.node(type=type, call=call,
	filename=filename)
    if (!is.null(comment))
      eventNode <- addChildren(eventNode, xmlTextNode(comment))
    trail <- addChildren(trail, eventNode)
    return(trail)
  }
}
