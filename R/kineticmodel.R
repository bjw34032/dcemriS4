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
## $Id: kineticmodel.R 332 2010-01-29 16:54:07Z bjw34032 $

kineticModel <- function(time, par, model="extended", aif="fritz.hansen") {

  if (!(is.numeric(par$kep)))
    par$kep <- par$ktrans/par$ve

  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
         },
         fritz.hansen = {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
         },
         stop("ERROR: AIF parameters must be specified!"))
  
  model.weinmann <- function(time, ktrans, kep, ...) {
    ## Convolution of Tofts & Kermode AIF with single-compartment model
    erg <- D * ktrans * ((a1 / (m1 - kep)) *
                         (exp(-(time * kep)) - exp(-(time * m1))) +
                         (a2 / (m2 - kep)) *
                         (exp(-(time * kep)) - exp(-(time * m2))))
    erg[time <= 0] <- 0
    return(erg)
  }
  
  model.extended <- function(time, vp, ktrans, kep, ...) {
    ## Extended Tofts & Kermode model including the concentration of
    ## contrast agent in the blood plasma (vp)
    Cp <- function(tt, D, a1, a2, m1, m2)
      D * (a1 * exp(-m1 * tt) + a2 * exp(-m2 * tt))
    
    erg <- vp * Cp(time, D, a1, a2, m1, m2) +
      D * ktrans * ((a1 / (m1 - kep)) *
                    (exp(-(time * kep)) - exp(-(time * m1))) +
                    (a2 / (m2 - kep)) *
                    (exp(-(time * kep)) - exp(-(time * m2))))
    erg[time <= 0] <- 0
    return(erg)
  }
  
  switch(aif,
         tofts.kermode = {
           D <- 0.1; a1 <- 3.99; a2 <- 4.78; m1 <- 0.144; m2 <- 0.0111
         },
         fritz.hansen = {
           D <- 1; a1 <- 2.4; a2 <- 0.62; m1 <- 3.0; m2 <- 0.016
         },
         stop("ERROR: Model not supported!"))
  
  switch(model,
         weinmann = {
           result <- model.weinmann(time, par$ktrans, par$kep)
         },
         extended = {
           result <- model.extended(time, par$vp, par$ktrans, par$kep)
         })
  
  return(result)
}

