if (getOption("NIfTI.audit.trail")) {
  nifti.extension.to.audit.trail <<- function(nim, filename=NULL, call=NULL) {
    if (!is(nim, "niftiAuditTrail")) {
      nim <- as(nim, "niftiAuditTrail")
    }
    # We enforce that there is a single extension with ecode == audit.trail.extension.ecode
    ecodes <- lapply(nim@extensions, function(x) { x@ecode })
    oei <- which(ecodes == audit.trail.extension.ecode)

    if (length(oei) == 0) {
      nim@trail <- audit.trail.created(filename=filename,call=call)
    } else {
      # One Trail
      if (length(oei) > 1) {
	warning("Found more than one extension with ecode == ", 
	    audit.trail.extension.ecode, " Appending to last trail only")
	oei <- oei[length(oei)]
      }
      oe <- nim@extensions[[oei]]@edata
      nim@extensions[[oei]] <- NULL

      nim@trail <- audit.trail.systemNode.event(xmlRoot(xmlTreeParse(oe, asText=TRUE)),type="read", filename=filename, call=call)

    }
    return(nim)
  }

  nifti.audit.trail.to.extension <<- function(nim, filename=filename, call=call) {
    sec <- new("niftiExtensionSection")
    sec@ecode <- audit.trail.extension.ecode
    nim@trail <- audit.trail.systemNode.event(nim@trail, "saved", filename=filename, call=call)

    # Serialize the XML to sec@edata
    ## DIRTY DIRTY DIRTY you should wash your eyes out after reading this.
    useFancyQuotes <- getOption("useFancyQuotes")
    options("useFancyQuotes"=FALSE)
    sec@edata <- toString.XMLNode(trail)
    options("useFancyQuotes"=useFancyQuotes)

    # Fix the esize to be congruent to 0 mod 16
    sec@esize <- nchar(sec@edata, type="bytes") + 8
    sec@esize <- (-sec@esize %% 16) + sec@esize
    return(sec)
  }

  audit.trail.system.node <- function(type="system-info", filename=NULL, call=NULL) {
    if (is(call,"call")) {
      call <- as.character(as.expression(call))
    }
    currentDateTime <- format(Sys.time(), "%a %b %d %X %Y %Z")
    system <- xmlNode(type, attrs=c("filename"=filename, "call"=call), namespace="")
    system <- addAttributes(system, "r-version"=version$version.string, "date"=currentDateTime, "user"=Sys.getenv("LOGNAME"), "dcemri-version"=packageDescription("dcemri")$Version)
    return(system)
  }

  audit.trail.created <<- function(history=NULL, call=NULL, filename=NULL) {
    trail <- new.audit.trail()
    created <- audit.trail.system.node("created", "filename"=filename, "call"=call)
    if (!is.null(history)) {
      historyNode <- xmlNode("history")
      lapply(xmlChildren(history), function(x) {
	    historyNode <<- addChildren(historyNode, x)
	  })

      created <- addChildren(created, historyNode)
    }

    trail <- addChildren(trail, created) 

    return(trail)
  }

  audit.trail.event <<- function(trail, type=NULL, call=NULL, comment=NULL) {
    if (is(call,"call")) {
      call <- as.character(as.expression(call))
    }
    eventNode <- xmlNode("event",attrs=c("type"=type, "call"=call))
    if (!is.null(comment)) {
      eventNode <- addChildren(eventNode, xmlTextNode(comment))
    }
    trail <- addChildren(trail, eventNode)
    return(trail)
  }

  audit.trail.systemNode.event <<- function(trail, type=NULL, call=NULL, filename=NULL, comment=NULL) {
    eventNode <- audit.trail.system.node(type=type,call=call,filename=filename)
    if (!is.null(comment)) {
      eventNode <- addChildren(eventNode, xmlTextNode(comment))
    }
    trail <- addChildren(trail, eventNode)
    return(trail)
  }
}
