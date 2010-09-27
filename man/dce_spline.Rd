% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dce_spline.Rd 332 2010-01-29 16:54:07Z bjw34032 $

\name{dcemri.spline}
\alias{dcemri.spline}
\alias{dcemri.spline,array-method}
\alias{dcemri.spline.single}
\title{Bayesian P-Splines for Dynamic Contrast-Enhanced MRI Data}
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
              samples=FALSE, multicore=FALSE, verbose=FALSE, 
              response=FALSE, fitted=FALSE, 
              ...) 
dcemri.spline.single(conc, time, D, time.input, p, rw, knots, k,
                     A, t0.compute=FALSE, nlr=FALSE, nriters=500,
                     thin=5, burnin=100, ab.hyper=c(1e-5,1e-5),
                     ab.tauepsilon=c(1,1/1000), silent=0,
                     multicore=FALSE, model=NULL,
                     model.func=NULL, model.guess=NULL,
                     samples=FALSE, B=NULL)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series (last
    dimension must be time).}
  \item{time}{Time in minutes.}
  \item{img.mask}{Mask matrix or array. Voxels with \code{mask=0} will be
    excluded.}
  \item{time.input}{Time in minutes for observed arterial
    input function (default = \sQuote{time}).}
  \item{aif}{is a character string that identifies the parameters of the arterial input function.  
    Acceptable values are: \dQuote{tofts.kermode}, \dQuote{fritz.hansen} or \dQuote{observed}.
    If \dQuote{observed} you must provide the observed concentrations in
    \sQuote{aif.observed}.}
  \item{aif.observed}{is the user-defined vector of arterial
    concentrations observed at \sQuote{time.input} (only for \sQuote{aif}=observed).}
  \item{multicore}{(logical) use the \pkg{multicore} package.}
  \item{verbose}{(logical) allows text-based feedback during execution
    of the function (default = \code{FALSE}).}
  \item{samples}{If \code{TRUE} output includes samples drawn
    from the posterior distribution for all parameters.}
  \item{nlr}{If \code{TRUE}, a response model is fitted to the estimated response function.}
  \item{model}{Only if \sQuote{nlr}=\code{TRUE} Response model fitted to the estimated response
    function.  Acceptable values include: \dQuote{AATH} or \dQuote{weinmann} (default).}
  \item{ab.hyper}{Hyper priors for adaptive smoothness parameter}
  \item{ab.tauepsilon}{Hyper-prior parameters for observation error
    Gamma prior.}
  \item{p}{Number of knots of B-Spline basis.}
  \item{t0.compute}{If \code{TRUE}, onset will be estimated from response function.}
  \item{k}{Order of B-Splines.}
  \item{knots}{Vector of knots. Use this if you need unequally spaced knots.}
  \item{rw}{Order of random walk prior. Acceptable values are 1 and 2.}
  \item{nriters}{Total number of iterations.}
  \item{thin}{Thining factor.}
  \item{burnin}{Number of iterations for burn-in.}
  \item{response}{If \code{TRUE} response functions per voxel are returned.}
  \item{fitted}{If \code{TRUE} fitted time curved per voxel are returned.}
}
\details{
  See Schmid \emph{et al.} (2009) for more details.
}
\value{
  The maximum of the response function \code{Fp} for the  masked region is provided by default. Where   
  appropriate, response functions, fitted functions, and parameter estimates (along with their standard 
  errors) are provided. All multi-dimensional arrays are provided in \code{nifti} format.
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
