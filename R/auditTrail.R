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
## $Id$
##

dcemri.info <- function(type) {
  switch(type,
         ecode = 1002,
         namespace = "http://www.dcemri.org/namespaces/audit-trail/1.0")
}

enableAuditTrail <- function() {
  if (require("XML")) {
    if (!isClass("niftiAuditTrail")) {
      options("NIfTI.audit.trail"=TRUE)
      setClass("niftiAuditTrail",
               representation(trail="XMLAbstractNode"),
               prototype(trail=newAuditTrail()),
               contains="niftiExtension")
    }
  }
}

## Look back through call stack for the last call with the name functionName
## otherwise return the function that called the function that called us

getLastCallWithName <- function(functionName) {
  theCalls <- sys.calls()
  correctCalls <- which(sapply(theCalls, function(x) x[[1]] == functionName))
  if (length(correctCalls) == 0) {
    return(theCalls[max(1, length(theCalls)-2)])
  }
  return(theCalls[[max(correctCalls)]])
}

newAuditTrail <- function() {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    trail <- newXMLNode("audit-trail",
                     namespaceDefinitions=dcemri.info("namespace"))
    return(trail)
  }
}

niftiExtensionToAuditTrail <- function(nim, workingDirectory=NULL,
                                       filename=NULL, call=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    if (!is(nim, "niftiAuditTrail"))
      nim <- as(nim, "niftiAuditTrail")
    ## We enforce that there is a single extension with ecode == dcemri.ecode
    ecodes <- lapply(nim@extensions, function(x) x@ecode)
    oei <- which(ecodes == dcemri.info("ecode"))
    
    if (length(oei) == 0) {
      audit.trail(nim) <-
        niftiAuditTrailCreated(workingDirectory=workingDirectory,
                               filename=filename, call=call)
    } else {
      ## One Trail
      if (length(oei) > 1) {
	warning("Found more than one extension with ecode ==",
                dcemri.info("ecode"), ", Appending to last trail only")
	oei <- oei[length(oei)]
      }
      oe <- nim@extensions[[oei]]@edata
      nim@extensions[[oei]] <- NULL
      audit.trail(nim) <-
        niftiAuditTrailSystemNodeEvent(xmlRoot(xmlParse(oe, asText=TRUE)),
                                       type="read",
                                       workingDirectory=workingDirectory,
                                       filename=filename, call=call)
    }
  }
  return(nim)
}

niftiAuditTrailToExtension <- function(nim, workingDirectory=NULL,
                                       filename=NULL, call=NULL) {
  if (getOption("NIfTI.audit.trail") && is(nim, "niftiAuditTrail")) {
    require("XML")
    sec <- new("niftiExtensionSection")
    sec@ecode <- dcemri.info("ecode")
    audit.trail(nim) <-
      niftiAuditTrailSystemNodeEvent(audit.trail(nim), "saved",
                                     workingDirectory=workingDirectory,
                                     filename=filename, call=call)
    ## Serialize the XML to sec@edata
    sec@edata <- saveXML(audit.trail(nim))

    ## Fix the esize to be congruent to 0 mod 16
    sec@esize <- nchar(sec@edata, type="bytes") + 8
    sec@esize <- (-sec@esize %% 16) + sec@esize
    return(sec)
  }
}

niftiAuditTrailSystemNode <- function(type="system-info",
                                      workingDirectory=NULL, filename=NULL,
                                      call=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    require("XML")
    if (is(call, "character") && is(try(get(call, mode="function"),
                                        silent=TRUE), "function")) 
      call <- as.character(as.expression(getLastCallWithName(call)))
    if (is(call, "call"))
      call <- as.character(as.expression(call))
    currentDateTime <- format(Sys.time(), "%a %b %d %X %Y %Z")
    system <- newXMLNode(type, attrs=c("workingDirectory"=workingDirectory,
                                 "filename"=filename, "call"=call),
                         namespace="")
    system <- addAttributes(system,
                            "r-version"=version$version.string,
                            "date"=currentDateTime,
                            "user"=Sys.getenv("LOGNAME"),
                            "dcemri-version"=packageDescription("dcemriS4")$Version)
    return(system)
  }
}

niftiAuditTrailSystemNodeEvent <- function(trail, type=NULL, call=NULL,
                                           workingDirectory=NULL,
                                           filename=NULL, comment=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    if (is(trail,"niftiAuditTrail")) {
      return(niftiAuditTrailSystemNodeEvent(audit.trail(trail), type, call,
                                            workingDirectory, filename,
                                            comment))
    }
    require("XML")
    eventNode <- niftiAuditTrailSystemNode(type=type, call=call,
                                           workingDirectory=workingDirectory,
                                           filename=filename)
    if (!is.null(comment))
      eventNode <- addChildren(eventNode, newXMLTextNode(comment))
    trail <- addChildren(trail, eventNode)
    return(trail)
  }
}

niftiAuditTrailCreated <- function(history=NULL, call=NULL,
                                   workingDirectory=NULL, filename=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    if (is(history, "niftiAuditTrail")) {
      return(niftiAuditTrailCreated(audit.trail(history), call,
                                    workingDirectory, filename))
    } else {
      require("XML")
      trail <- newAuditTrail()
      if (is.null(history) || length(xmlChildren(history)) == 0) {
	created <-
          niftiAuditTrailSystemNode("created",
                                    "workingDirectory"=workingDirectory,
                                    "filename"=filename, "call"=call)
      } else {
	historyChildren <- xmlChildren(history)

	lastEvent <- historyChildren[[length(historyChildren)]]
	
	if (xmlName(lastEvent) == "event" && xmlAttrs(lastEvent)[["type"]] == "processing") {
	  # We are in some processing history; the given call is not the correct call
	  call <- xmlAttrs(lastEvent)[["call"]]
	  historyChildren <- historyChildren[1:(length(historyChildren) - 1)]
	} 

	created <-
          niftiAuditTrailSystemNode("created",
                                    "workingDirectory"=workingDirectory,
                                    "filename"=filename, "call"=call)
	historyNode <- newXMLNode("history")
	## OK, serialize and reParse the history
	historyChildren <-
          lapply(historyChildren,
                 function(x) xmlRoot(xmlParse(saveXML(x), asText=TRUE)))
	historyNode <- addChildren(historyNode, historyChildren)
	created <- addChildren(created, historyNode)
      }
      trail <- addChildren(trail, created) 
      return(trail)
    }
  }
}

niftiAuditTrailEvent <- function(trail, type=NULL, call=NULL, comment=NULL) {
  if (getOption("NIfTI.audit.trail")) {
    if (is(trail,"niftiAuditTrail")) {
      return(niftiAuditTrailEvent(audit.trail(trail), type, call, comment))
    }
    require("XML")
    if (is(call, "character") && is(try(get(call, mode="function"),
                                        silent=TRUE),"function")) 
      call <- as.character(as.expression(getLastCallWithName(call)))
    if (is(call, "call"))
      call <- as.character(as.expression(call))
    currentDateTime <- format(Sys.time(), "%a %b %d %X %Y %Z")
    eventNode <- newXMLNode("event", attrs=c("type"=type, "call"=call,
                                       "date"=currentDateTime))
    if (!is.null(comment))
      eventNode <- addChildren(eventNode, newXMLTextNode(comment))
    trail <- addChildren(trail, eventNode)
    return(trail)
  }
}

