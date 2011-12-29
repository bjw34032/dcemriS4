##
##
## Copyright (c) 2009-2011 Brandon Whitcher and Volker Schmid
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
## $Id:$
##

conv.fft <- function(A, B, C, FFTA=NULL) {
  if (length(dim(A)) == 3) {
    if (length(dim(A)) == length(dim(B)) && length(dim(A)) == length(C)) {
      X <- nrow(A)
      Y <- ncol(A)
      Z <- nsli(A)
      if (is.null(FFTA)) {
        out <- Re(fft(fft(A) * Conj(fft(B)), inverse=TRUE))[X:1,Y:1,Z:1,drop=FALSE] / (X*Y*Z)
      } else {
        out <- Re(fft(FFTA * Conj(fft(B)), inverse=TRUE))[X:1,Y:1,Z:1,drop=FALSE] / (X*Y*Z)
      }
      if (X > 1) {
        out <- out[c((X-C[1]+1):X, 1:(X-C[1])),,,drop=FALSE]
      }
      if (Y > 1) {
        out <- out[,c((Y-C[2]+1):Y, 1:(Y-C[2])),,drop=FALSE]
      }
      if (Z > 1) {
        out <- out[,,c((Z-C[3]+1):Z, 1:(Z-C[3])),drop=FALSE]
      }
      return(out)
    } else {
      stop("Objects are not all the same dimension!")
    }
  } else {
    stop("Only three-dimensional objects are allowed!")
  }
}

find.center <- function(M) {
  if (!is.logical(M)) {
    stop("Object must be logical!")
  }
  if (length(dim(M)) != 3) {
    stop("Object must be three-dimensional!")
  }
  X <- nrow(M)
  Y <- ncol(M)
  Z <- nsli(M)
  xx <- array(1:X, dim(M))
  yy <- array(rep(1:Y, each=X), dim(M))
  zz <- array(rep(1:Z, each=X*Y), dim(M))
  center <- c(mean(xx[M]), mean(yy[M]), mean(zz[M]))
  trunc(center)
}

### FIXME Should this be genericised? 
shift3D <- function(A, s, type, fill=0) {
  if (length(dim(A)) != 3) {
    stop("Object must be three-dimensional!")
  }
  X <- nrow(A)
  Y <- ncol(A)
  Z <- nsli(A)
  if (s != 0) {
    if (type == "LR") {
      ## left-right
      if (s > 0) {
        A <- A[c((X-s+1):X,1:(X-s)),,]
        A[1:s,,] <- fill
      } else {
        A <- A[c((abs(s)+1):X,1:abs(s)),,]
        A[(X-abs(s)+1):X,,] <- fill
      }
    } else {
      if (type == "AP") {
        ## anterior-posterior
        if (s > 0) {
          A <- A[,c((Y-s+1):Y,1:(Y-s)),]
          A[,1:s,] <- fill
        } else {
          A <- A[,c((abs(s)+1):Y,1:abs(s)),]
          A[,(Y-abs(s)+1):Y,] <- fill
        }
      } else {
        if (type == "SI") {
          ## superior-inferior
          if (s > 0) {
            A <- A[,,c((Z-s+1):Z,1:(Z-s))]
            A[,,1:s] <- 0
          } else {
            A <- A[,,c((abs(s)+1):Z,1:abs(s))]
            A[,,(Z-abs(s)+1):Z] <- fill
          }
        } else {
          stop("Type of translation not recognized!")
        }
      }
    }
  }
  return(A)
}

#############################################################################
## setGeneric("ftm")
#############################################################################

setGeneric("ftm", function(input, ...) standardGeneric("ftm"))
setMethod("ftm", signature(input="array"),
          function(input, ...) .dcemriWrapper("ftm", input, ...))

.ftm <- function(input, mask, reference, plot=TRUE, ...) {
  ## Fast template matching via cross-correlation
  W <- ntim(input)
  if (ntim(input) < 1) {
    stop("4D object is assumed!")
  }
  template <- ifelse(mask > 0, reference, 0)
  templateFFT <- fft(template)
  maskFFT <- fft(mask)
  tc <- find.center(ifelse(template > 0, TRUE, FALSE))
  templateSS <- conv.fft(mask, template^2, tc, FFTA=maskFFT)
  ##numerator <- localSS <- array(0, dim(input))
  localCOR <- array(0, dim(input))
  for (w in 1:W) {
    target <- input[,,,w]
    ##numerator[,,,w] <- conv.fft(template, target, tc, FFTA=templateFFT)
    ##localSS[,,,w] <- conv.fft(mask, target^2, tc, FFTA=maskFFT)
    numerator <- conv.fft(template, target, tc, FFTA=templateFFT)
    localSS <- conv.fft(mask, target*target, tc, FFTA=maskFFT)
    localCOR[,,,w] <- numerator / sqrt(templateSS[tc[1],tc[2],tc[3]]) / sqrt(localSS)
  }
  ## localCOR <- numerator / sqrt(templateSS[tc[1],tc[2],tc[3]]) / sqrt(localSS)
  localCOR[!is.finite(localCOR)] <- NA
  wmax.localCOR <- apply(localCOR, 4,
                         function(x) {
                           which(x == max(x,na.rm=TRUE), arr.ind=TRUE)
                         })
  ## which(localCOR[,,,w] == max(localCOR[,,,w], na.rm=TRUE), arr.ind=TRUE)
  offset <- tc - wmax.localCOR
  output <- input@.Data
  for (w in 1:W) {
    output[,,,w] <- shift3D(output[,,,w], offset[1,w], type="LR")
    output[,,,w] <- shift3D(output[,,,w], offset[2,w], type="AP")
    output[,,,w] <- shift3D(output[,,,w], offset[3,w], type="SI")
  }
  storage.mode(output) <- storage.mode(input)
  as(output, "nifti") <- input
  if (plot) {
    matplot(1:ncol(offset), t(offset), type="l", lwd=2, lty=1,
            ylim=range(c(range(offset),-10,10)), xlab="", ylab="voxels",
            main="Motion Correction via Template Matching")
    legend("topleft", c("X","Y","Z"), col=1:3, lty=1, lwd=2, bty="n")
  }
  list(out=output, offset=offset, t.center=tc, localCOR=localCOR, plot=NA)
}

