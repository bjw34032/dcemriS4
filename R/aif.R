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
##
## Time-stamp: <2009-07-14 09:47:30 (bjw34032)>
## $Id: aif.R 332 2010-01-29 16:54:07Z bjw34032 $
##

aif.orton.exp <- function(tt, AB, muB, AG, muG) {
  out <- AB * tt * exp(-muB * tt) + AG * (exp(-muG * tt) - exp(-muB * tt))
  out[tt < 0] <- 0
  return(out)
}

orton.exp.lm <- function(tt, aif,
                         guess=c(log(100), log(10), log(1), log(0.1)),
                         nprint=0) {
  func <- function(x, aparams, aif) {
    AB <- aparams[1]
    muB <- aparams[2]
    AG <- aparams[3]
    muG <- aparams[4]
    return(aif - aif.orton.exp(x, AB, muB, AG, muG))
  }
  out <- minpack.lm::nls.lm(par=guess, fn=func, control=list(nprint=nprint),
                            x=tt, aif=aif)
  list(AB=out$par[1], muB=out$par[2], AG=out$par[3], muG=out$par[4], 
       info=out$info, message=out$message)
}

model.orton.exp <- function(tt, aparams, kparams) {
  ## Extended model using the exponential AIF from Matthew Orton (ICR)
  Cp <- function(tt, ...) {
    AB * tt * exp(-muB * tt) + AG * (exp(-muG * tt) - exp(-muB * tt))
  }
  aparams <- as.numeric(aparams)
  kparams <- as.numeric(kparams)
  AB <- aparams[1]
  muB <- aparams[2]
  AG <- aparams[3]
  muG <- aparams[4]
  vp <- kparams[1]
  ktrans <- kparams[2]
  kep <- kparams[3]
  
  T1 <- AB * kep / (kep - muB)
  T2 <- tt * exp(-muB * tt) -
    (exp(-muB * tt) - exp(-kep * tt)) / (kep - muB)
  T3 <- AG * kep
  T4 <- (exp(-muG * tt) - exp(-kep * tt)) / (kep - muG) -
    (exp(-muB * tt) - exp(-kep * tt)) / (kep - muB)
  
  out <- vp * Cp(tt) + ktrans * (T1 * T2 + T3 * T4)
  out[tt <= 0] <- 0
  return(out)
}

extract.aif <- function(img, x, y, z, thresh=0.9) {
  c.start <- function(ctc) {
    if (sum(is.na(ctc)) > 0) {
      return(0)
    } else {
      if (sd(ctc) == 0) {
        return(0)
      } else {
        out <- cor(ctc, start, use="pairwise.complete.obs")
        return(out)
      }
    }
  }
  
  check <- function(xx, yy, zz, aif.mask, thresh) {
    if (xx != 1 && aif.mask[xx-1,yy,zz] == 0) {
      if (c.test[xx-1,yy,zz] > thresh) {
        aif.mask[xx-1,yy,zz] <- 1
        aif.mask <- check(xx-1, yy, zz, aif.mask, thresh)
      }
    }
    if (xx != X && aif.mask[xx+1,yy,zz] == 0) {
      if (c.test[xx+1,yy,zz] > thresh) {
        aif.mask[xx+1,yy,zz] <- 1
        aif.mask <- check(xx+1, yy, zz, aif.mask, thresh)
      }
    }
    if (yy != 1 && aif.mask[xx,yy-1,zz] == 0) {
      if (c.test[xx,yy-1,zz] > thresh) {
        aif.mask[xx,yy-1,zz] <- 1
        aif.mask <- check(xx, yy-1, zz, aif.mask, thresh)
      }
    }
    if (yy != Y && aif.mask[xx,yy+1,zz] == 0) {
      if (c.test[xx,yy+1,zz] > thresh) {
        aif.mask[xx,yy+1,zz] <- 1
        aif.mask <- check(xx, yy+1, zz, aif.mask, thresh)
      }
    }
    if (zz != 1 && aif.mask[xx,yy,zz-1] == 0) {
      if (c.test[xx,yy,zz-1] > thresh) {
        aif.mask[xx,yy,zz-1] <- 1
        aif.mask <- check(xx, yy, zz-1, aif.mask, thresh)
      }
    }
    if (zz != Z && aif.mask[xx,yy,zz+1] == 0) {
      if (c.test[xx,yy,zz+1] > thresh) {
        aif.mask[xx,yy,zz+1] <- 1
        aif.mask <- check(xx, yy, zz+1, aif.mask, thresh)
      }
    }
    return(aif.mask)
  }

  X <- dim(img)[1]
  Y <- dim(img)[2]
  Z <- dim(img)[3]
  W <- dim(img)[4]
  
  start <- img[x,y,z,]
  c.test <- apply(img, 1:3, c.start)
  
  aif.mask <- array(0, c(X,Y,Z))
  aif.mask[x,y,z] <- 1
  
  aif.mask <- check(x, y, z, aif.mask, thresh)
  n <- sum(aif.mask, na.rm=TRUE)
  test <- array(NA, c(n,W))
  coord <- array(NA, c(n,3))
  l <- 0
  for (i in 1:X) {
    for (j in 1:Y) {
      for (k in 1:Z) {
        if (!is.na(c.test[i,j,k]) && aif.mask[i,j,k] == 1) {
          l <- l + 1
          test[l,] <- img[i,j,k,]
          coord[l,] <- c(i,j,k)
        }
      }
    }
  }
  if (l != n) {
    return(FALSE)
  }
  list("coord"=coord, "conc"=test, "mask"=aif.mask, "cor"=c.test)
}

