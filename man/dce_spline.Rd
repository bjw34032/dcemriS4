% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dce_spline.Rd 205 2009-08-26 13:57:14Z bjw34032 $

\name{dcemri.spline}
\alias{dcemri.spline}
\alias{dcemri.spline,array-method}
\alias{dcemri.spline.single}
\title{Bayesian P-Splines for Dynamic Contrasat-Enhanced MRI Data}
\description{
  Quantitative analysis of DCE-MRI typically involves the convolution of
  an arterial input function (AIF) with a nonlinear pharmacokinetic
  model of the contrast agent concentration.  This function takes a
  semi-parametric penalized spline smoothing approach, with which the
  AIF is convolved with a set of B-splines to produce a design matrix
  using locally adaptive smoothing parameters based on Bayesian
  penalized spline models (P-splines).
} 
\usage{
\S4method{dcemri.spline}{array}(conc, time, img.mask, time.input=time,
              model="weinmann", aif="tofts.kermode", user=NULL,
              aif.observed=NULL, nriters=500, thin=5,
              burnin=100, ab.hyper=c(1e-5,1e-5),
              ab.tauepsilon=c(1,1/1000), k=4, p=25, rw=2,
              knots=NULL, nlr=FALSE, t0.compute=FALSE,
              samples=FALSE, multicore=FALSE, verbose=FALSE, ...) 
dcemri.spline.single(conc, time, D, time.input, p, rw, knots, k,
                     A, t0.compute=FALSE, nlr=FALSE, nriters=500,
                     thin=5, burnin=100, ab.hyper=c(1e-5,1e-5),
                     ab.tauepsilon=c(1,1/1000), silent=0,
                     multicore=FALSE, model=NULL,
                     model.func=NULL, model.guess=NULL,
                     samples=FALSE, B=NULL)
}
\arguments{
  \item{conc}{An array of Gd concentration}
  \item{time}{Time points of aquisition of Gd concentration}
  \item{img.mask}{Array of voxels to fit}
  \item{time.input}{Time points of observed Arterial Gd concentration,
    defaults to time}
  \item{model}{The type of compartmental model to be used. Acceptable
    models include: \dQuote{AATH} or \dQuote{weinmann} (default).}
  \item{aif}{Arterial input function to use. Values include:
    \dQuote{tofts.kermode}, \dQuote{fritz.hansen} or \dQuote{observed}.
    If \dQuote{observed} you must provide the observed concentrations in
    aif.observed.}
  \item{aif.observed}{Arterial concentrations observed at timepoints
    time.input}
  \item{multicore}{Use multicore library}
  \item{verbose}{}
  \item{nlr}{Return the generated samples}
  \item{user}{}
  \item{ab.hyper}{}
  \item{ab.tauepsilon}{}
  \item{p}{}
  \item{t0.compute}{}
  \item{samples}{}
  \item{k}{}
  \item{knots}{}
  \item{rw}{}
  \item{nriters}{}
  \item{thin}{}
  \item{burnin}{}
  \item{D}{}
  \item{B}{}
  \item{A}{}
  \item{silent}{}
  \item{model.func}{}
  \item{model.guess}{}
  \item{...}{}
}
\details{
  See Schmid \emph{et al.} (2009) for more details.
}
\value{
  To be added.
}
\references{
  Schmid, V., Whitcher, B., Padhani, A.R. and G.-Z. Yang (2009) A
  semi-parametric technique for the quantitative analysis of dynamic
  contrast-enhanced MR images based on Bayesian P-splines, \emph{IEEE
    Transactions on Medical Imaging}, \bold{28} (6), 789-798.
}
\seealso{
  \code{\link{dcemri.bayes}}, \code{\link{dcemri.lm}},
  \code{\link{dcemri.map}}
}
\examples{
data("buckley")
xi <- seq(5, 300, by=5)
img <- array(t(breast$data)[,xi], c(13,1,1,60))
mask <- array(TRUE, dim(img)[1:3])
time <- buckley$time.min[xi]

## Generate AIF params using the orton.exp function from Buckley's AIF
aif <- buckley$input[xi]

fit.spline <- dcemri.spline(img, time, mask, aif="fritz.hansen",
                            nriters=250, nlr=TRUE)
fit.spline.aif <- dcemri.spline(img, time, mask, aif="observed",
                                aif.observed=aif, nriters=250,
                                nlr=TRUE)

plot(breast$ktrans, fit.spline$ktrans, xlim=c(0,1), ylim=c(0,1),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans})))
points(breast$ktrans, fit.spline.aif$ktrans, pch=2)
abline(0, 1, lwd=1.5, col="red")
legend("right", c("fritz.hansen", "observed"), pch=1:2)

cbind(breast$ktrans, fit.spline$ktrans[,,1], fit.spline.aif$ktrans[,,1])
}
\author{Volker Schmid <\email{volkerschmid@users.sourceforge.net}>}
\keyword{models}
