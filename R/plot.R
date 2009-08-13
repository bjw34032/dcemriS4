#############################################################################
## image() for class="nifti"
#############################################################################

setMethod("image", signature(x="nifti"),
          function(x, z=1, w=1, col=gray(0:64/64),
                   plot.type=c("multiple","single"), zlim=NULL,
                   xlab="", ylab="", axes=FALSE, ...) {
            X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
            if (X == 0 || Y == 0 || Z == 0)
              stop("size of NIfTI volume is zero, nothing to plot")
            if (is.null(zlim))
              zlim <- c(x@"cal_min", x@"cal_max")
            if (is.na(W)) {
              ## three-dimensional array
              if (z < 1 || z > Z)
                stop("slice \"z\" out of range")
              oldpar <- par()
              if (plot.type[1] == "multiple")
                index <- 1:Z
              else
                index <- z
              lz <- length(index)
              par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
              for (z in index)
                graphics:::image(1:X, 1:Y, x[,,z], col=col, zlim=zlim,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
            } else {
              ## four-dimensional array
              if (w < 1 || w > W)
                stop("volume \"w\" out of range")
              if (z < 1 || z > Z)
                stop("slice \"z\" out of range")
              if (plot.type[1] == "multiple")
                index <- 1:Z
              else
                index <- z
              lz <- length(index)
              par(mfrow=ceiling(rep(sqrt(lz),2)), mar=rep(0,4))
              for (z in index)
                graphics:::image(1:X, 1:Y, x[,,z,w], col=col, zlim=zlim,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
            }
            
          })

#############################################################################
## overlay() for class="nifti"
#############################################################################

setGeneric("overlay", function(x, y, ...) standardGeneric("overlay"))
setMethod("overlay", signature(x="nifti", y="nifti"),
          function(x, y, z=1, w=1, col.x=gray(0:64/64), col.y=hotmetal(),
                   zlim.x=NULL, zlim.y=NULL, plot.type=c("multiple","single"),
                   xlab="", ylab="", axes=FALSE, ...) {
            ## should do plot.type=c("multiple","single") like plot.ts
            if (!all(dim(x)[1:3] == dim(y)[1:3]))
              stop("dimensions of \"x\" and \"y\" must be equal")
            X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
            if (X == 0 || Y == 0 || Z == 0)
              stop("size of NIfTI volume is zero, nothing to plot")
            if (is.null(zlim.x))
              zlim.x <- c(x@"cal_min", x@"cal_max")
            if (is.null(zlim.y))
              zlim.y <- c(y@"cal_min", y@"cal_max")
            if (is.na(W)) {
              ## three-dimensional array
              oldpar <- par()
              if (plot.type[1] == "single") {
                if (z < 1 || z > Z)
                  stop("slice \"z\" out of range")
                par(mfrow=c(1,1), mar=rep(0,4))
                graphics:::image(1:X, 1:Y, x[,,z], col=col.x, zlim=zlim.x,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
                graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                 add=TRUE)
              } else {
                par(mfrow=c(ceiling(sqrt(Z)), ceiling(sqrt(Z))), mar=rep(0,4))
                for (z in 1:Z) {
                  graphics:::image(1:X, 1:Y, x[,,z], col=col.x, zlim=zlim.x,
                                   axes=axes, xlab=xlab, ylab=ylab, ...)
                  graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                   add=TRUE)
                }
              }
            } else {
              ## four-dimensional array
              oldpar <- par()
              if (w < 1 || w > W)
                stop("volume \"w\" out of range")
              if (z < 1 || z > Z)
                stop("slice \"z\" out of range")
              if (plot.type == "single") {
                par(mfrow=c(1,1), mar=rep(0,4))
                graphics:::image(1:X, 1:Y, x[,,z,w], col=col.x, zlim=zlim.x,
                                 axes=axes, xlab=xlab, ylab=ylab, ...)
                graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                 add=TRUE)
              } else {
                par(mfrow=c(ceiling(sqrt(Z)), ceiling(sqrt(Z))), mar=rep(0,4))
                for (z in 1:Z) {
                  graphics:::image(1:X, 1:Y, x[,,z,w], col=col.x, zlim=zlim.x,
                                   axes=axes, xlab=xlab, ylab=ylab, ...)
                  graphics:::image(1:X, 1:Y, y[,,z], col=col.y, zlim=zlim.y,
                                   add=TRUE)
                }
              }
              graphics:::image(1:X, 1:Y, x[,,z,w], col=colx, zlim=zlimx,
                               axes=axes, xlab=xlab, ylab=ylab, ...)
              graphics:::image(1:X, 1:Y, y[,,z], col=coly, zlim=zlimy,
                               add=TRUE)
            }
            par(oldpar)
            invisible()
          })

#############################################################################
## orthographic() for class="nifti"
#############################################################################

setGeneric("orthographic", function(x, ...) standardGeneric("orthographic"))
setMethod("orthographic", signature(x="nifti"),
          function(x, xyz=NULL, crosshairs=TRUE, col.crosshairs="red", w=1,
                   zlim=NULL, col=gray(0:64/64), xlab="", ylab="",
                   axes=FALSE, ...) {
            X <- nrow(x) ; Y <- ncol(x) ; Z <- dim(x)[3] ; W <- dim(x)[4]
            if (is.null(xyz))
              xyz <- ceiling(c(X,Y,Z)/2)
            if (X == 0 || Y == 0 || Z == 0)
              stop("size of NIfTI volume is zero, nothing to plot")
            if (is.null(zlim))
              zlim <- c(x@"cal_min", x@"cal_max")
            if (is.na(W)) {
              ## three-dimensional array
              oldpar <- par()
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
            } else {
              ## four-dimensional array
              if (w < 1 || w > W)
                stop("volume \"w\" out of range")
              oldpar <- par()
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
          })
