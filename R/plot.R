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

#############################################################################
## image() for class="nifti"
## promptMethods(image, "image-methods.Rd")
#############################################################################

image.nifti <- function(x, z=1, w=1, col=gray(0:64/64),
                        plot.type=c("multiple","single"), zlim=NULL,
                        xlab="", ylab="", axes=FALSE, ...) {
  ## set dimensions
  X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0)
    stop("size of NIfTI volume is zero, nothing to plot")
  ## check for z-limits; use internal by default
  if (is.null(zlim)) {
    zlim <- c(x@"cal_min", x@"cal_max")
    if (max(zlim) == 0)
      zlim <- c(x@"glmin", x@"glmax")
  }
  ## single or multiple images?
  if (plot.type[1] == "multiple")
    index <- 1:Z
  else
    index <- z
  lz <- length(index)
  ## plotting
  if (is.na(W)) { # three-dimensional array
    if (z < 1 || z > Z)
      stop("slice \"z\" out of range")
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
    for (z in index)
      graphics:::image(1:X, 1:Y, x[,,z], col=col, zlim=zlim,
                       axes=axes, xlab=xlab, ylab=ylab, ...)
  } else { # four-dimensional array
    if (w < 1 || w > W)
      stop("volume \"w\" out of range")
    if (z < 1 || z > Z)
      stop("slice \"z\" out of range")
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
    for (z in index)
      graphics:::image(1:X, 1:Y, x[,,z,w], col=col, zlim=zlim,
                       axes=axes, xlab=xlab, ylab=ylab, ...)
  }
  par(oldpar)
  invisible()
}

setMethod("image", signature(x="nifti"), image.nifti)
setMethod("image", signature(x="anlz"), image.nifti)

#############################################################################
## overlay() for class="nifti"
#############################################################################

overlay.nifti <- function(x, y, z=1, w=1, col.x=gray(0:64/64),
                          col.y=hotmetal(), zlim.x=NULL, zlim.y=NULL,
                          plot.type=c("multiple","single"),
                          xlab="", ylab="", axes=FALSE, ...) {
  ## both volumes must have the same dimension
  if (!all(dim(x)[1:3] == dim(y)[1:3]))
    stop("dimensions of \"x\" and \"y\" must be equal")
  ## set dimensions
  X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0)
    stop("size of NIfTI volume is zero, nothing to plot")
  ## check for z-limits in x; use internal by default
  if (is.null(zlim.x)) {
    zlim.x <- c(x@"cal_min", x@"cal_max")
    if (max(zlim.x) == 0)
      zlim.x <- c(x@"glmin", x@"glmax")
  }
  ## check for z-limits in y; use internal by default
  if (is.null(zlim.y)) {
    zlim.y <- c(y@"cal_min", y@"cal_max")
    if (max(zlim.y) == 0)
      zlim.y <- c(x@"glmin", x@"glmax")
  }
  if (plot.type[1] == "multiple")
    index <- 1:Z
  else
    index <- z
  lz <- length(index)
  if (is.na(W)) { # three-dimensional array
    if (z < 1 || z > Z)
      stop("slice \"z\" out of range")
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z], col=col.x, zlim=zlim.x, axes=axes,
                      xlab=xlab, ylab=ylab, ...)
      graphics::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y, add=TRUE)
    }
  } else { # four-dimensional array
    if (w < 1 || w > W)
      stop("volume \"w\" out of range")
    if (z < 1 || z > Z)
      stop("slice \"z\" out of range")
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
    for (z in index) {
      graphics::image(1:X, 1:Y, x[,,z,w], col=col.x, zlim=zlim.x,
                      axes=axes, xlab=xlab, ylab=ylab, ...)
      graphics::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y, add=TRUE)
    }
  }
  par(oldpar)
  invisible()
}

setGeneric("overlay", function(x, y, ...) standardGeneric("overlay"))
setMethod("overlay", signature(x="nifti", y="nifti"), overlay.nifti)
setMethod("overlay", signature(x="anlz", y="anlz"), overlay.nifti)
setMethod("overlay", signature(x="anlz", y="nifti"), overlay.nifti)
setMethod("overlay", signature(x="nifti", y="anlz"), overlay.nifti)
setMethod("overlay", signature(x="array", y="array"), overlay.nifti)
setMethod("overlay", signature(x="array", y="nifti"), overlay.nifti)
setMethod("overlay", signature(x="nifti", y="array"), overlay.nifti)
setMethod("overlay", signature(x="array", y="anlz"), overlay.nifti)
setMethod("overlay", signature(x="anlz", y="array"), overlay.nifti)

#############################################################################
## orthographic() for class="nifti"
#############################################################################

orthographic.nifti <- function(x, xyz=NULL, crosshairs=TRUE,
                               col.crosshairs="red", w=1, zlim=NULL,
                               col=gray(0:64/64), xlab="", ylab="",
                               axes=FALSE, ...) {
  X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
  ## Center crosshairs if not specified
  if (is.null(xyz))
    xyz <- ceiling(c(X,Y,Z)/2)
  ## check dimensions
  if (X == 0 || Y == 0 || Z == 0)
    stop("size of NIfTI volume is zero, nothing to plot")
  ## check for z-limits in x; use internal by default
  if (is.null(zlim)) {
    zlim <- c(x@"cal_min", x@"cal_max")
    if (max(zlim) == 0)
      zlim <- c(x@"glmin", x@"glmax")
  }
  if (is.na(W)) { # three-dimensional array
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=c(2,2), mar=rep(0,4))
    graphics::image(1:X, 1:Z, x[,xyz[2],], col=col,
                    asp=x@pixdim[4]/x@pixdim[2],
                    xlab=ylab, ylab=xlab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
    graphics::image(1:Y, 1:Z, x[xyz[1],,], col=col,
                    asp=x@pixdim[4]/x@pixdim[3],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
    graphics::image(1:X, 1:Y, x[,,xyz[3]], col=col,
                    asp=x@pixdim[3]/x@pixdim[2],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
  } else { # four-dimensional array
    if (w < 1 || w > W)
      stop("volume \"w\" out of range")
    oldpar <- par(no.readonly=TRUE)
    par(mfrow=c(2,2), mar=rep(0,4))
    graphics::image(1:X, 1:Z, x[,xyz[2],,w], col=col,
                    asp=x@pixdim[4]/x@pixdim[2],
                    xlab=ylab, ylab=xlab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[1], col=col.crosshairs)
    graphics::image(1:Y, 1:Z, x[xyz[1],,,w], col=col,
                    asp=x@pixdim[4]/x@pixdim[3],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[3], v=xyz[2], col=col.crosshairs)
    graphics::image(1:X, 1:Y, x[,,xyz[3],w], col=col,
                    asp=x@pixdim[3]/x@pixdim[2],
                    xlab=xlab, ylab=ylab, axes=axes, ...)
    if (crosshairs)
      abline(h=xyz[2], v=xyz[1], col=col.crosshairs)
  }
  par(oldpar)
  invisible()
}

setGeneric("orthographic", function(x, ...) standardGeneric("orthographic"))
setMethod("orthographic", signature(x="nifti"), orthographic.nifti)
setMethod("orthographic", signature(x="anlz"), orthographic.nifti)
setMethod("orthographic", signature(x="array"),
          function(x) {
            x <- as(x, "nifti")
            orthographic.nifti(x)
          })
