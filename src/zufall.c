/*
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
*/

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

double nulleins() {
  return runif(0.0, 1.0);
}

double RNDGAM(double a, double b) {
  int accept = 0;
  double c, d, u, v, w, x, y, z;
  if (a > 1) {                   /*   Algorithmus S.410 Devroye */
    c = a - 1.0;
    d = 3.0 * a - 3.0/4.0;
    while (accept == 0) {
      u = nulleins();
      v = nulleins();
      w = u * (1 - u);
      y = sqrt(d/w) * (u - 0.5);
      x = c + y;
      if (x >= 0) {
	z = 64.0 * w * w * w * v * v;
	if (z <= (1.0 - (2.0 * y * y / x))) {
	  accept = 1; 
	}
	if (accept == 0) {
	  if (log(z) <= 2.0 * ((c * log(x/c)) - y)) { 
	    accept = 1;
	  }
	}
      }
    }
  } else {                            /* Fall: a<=1; Stuart's theorem */
    x = pow(nulleins(), 1.0/a) * RNDGAM(a + 1.0, 1.0);
  }
  return x/b;
}

double reins() {
  return nulleins() * 2.0 - 1.0;
}

double normal(double m, double s) {
  return rnorm(0.0, 1.0) * sqrt(s) + m;
}

//Erzeugt Normalverteilten Zufallsvektor der Laenge noa
void gausssample(double* temp, int* noa) {
  int i;
  for (i=0; i < *noa; i++) {
    temp[i] = rnorm(0.0, 1.0);
  }
  return;
}
