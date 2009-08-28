% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dcemri.Rd 205 2009-08-26 13:57:14Z bjw34032 $

\name{dcemri.lm}
\alias{dcemri.lm}
\alias{dcemri.lm,array-method}
\title{Pharmacokinetic Models for Dynamic Contrast-Enhanced MRI Data}
\description{
  Parameter estimation for single compartment models is performed using
  literature-based or user-specified arterial input functions.  The
  Levenburg-Marquardt algorithm does the heavy lifting.
}
\usage{
\S4method{dcemri.lm}{array}(conc, time, img.mask, model="extended", aif=NULL,
          nprint=0, user=NULL, verbose=FALSE, ...)
}
\arguments{
  \item{conc}{is a multidimensional (1D-4D) array of contrast agent
    concentrations.  The last dimension is assumed to be temporal, while
    the previous dimensions are assued to be spatial.}
  \item{time}{is a vector of acquisition times (in minutes) relative to 
    injection of the contrast agent.  Negative values should be used
    prior to the injection.}
  \item{img.mask}{is a (logical) multidimensional array that identifies
    the voxels to be analyzed. Has to have same dimension as \code{conc}
    minus temporal dimension.}
  \item{model}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \itemize{
      \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with single
	compartment model}
      \item{\dQuote{extended}}{Weinmann model extended with additional
	vascular compartment (default)}
      \item{\dQuote{orton.exp}}{Extended model using Orton's exponential
	AIF} 
      \item{\dQuote{orton.cos}}{Extended model using Orton's raised
	cosine AIF}
    }
  }
  \item{aif}{is a character string that identifies the parameters of the
    type of arterial input function (AIF) used with the above model.
    Acceptable values are: \code{tofts.kermode} (default) or
    \code{fritz.hansen} for the \code{weinmann} and \code{extended}
    models; \code{orton.exp} (default) and \code{user} for the
    \code{orton.exp} model; \code{orton.cos} (default) and \code{user}
    for the \code{orton.cos} model.  All AIF models set the parametric
    form and parameter values -- except \code{user} where a set of
    user-defined parameter values are allowed.
  }
  \item{nprint}{is an integer, that enables controlled printing of
    iterates if it is positive.  In this case, estimates of \code{par}
    are printed at the beginning of the first iteration and every
    \code{nprint} iterations thereafter and immediately prior to return.
    If \code{nprint} is not positive, no tracing information on the
    progress of the optimization is produced.}
  \item{user}{is a list with the following paramters required: D, AB,
    muB, AG, muG.}
  \item{verbose}{.}
  \item{...}{Additional parameters to the function.}
}
\details{
  Compartmental models are the solution to the modified general rate
  equation (Kety 1951).  The specific parametric models considered here
  include the basic Kety model
  \deqn{C_t(t)=K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right],}
  where \eqn{\otimes} is the convoluation operator, and the so-called
  extended Kety model
  \deqn{C_t(t)=v_pC_p(t)+K^{trans}\left[C_p(t)\otimes\exp(-k_{ep}t)\right].}
  The arterial input function must be either literature-based (with
  fixed parameters) or the exponential AIF from Orton \emph{et al.}
  (2008) with user-defined parameters.
}
\value{
  Parameter estimates and their standard errors are provided for the
  masked region of the multidimensional array.  They include
  \item{ktrans}{Transfer rate from plasma to the extracellular,
    extravascular space (EES).}
  \item{kep}{Rate parameter for transport from the EES to plasma.}
  \item{ve}{Fractional occupancy by EES (the ratio between ktrans and
    kep).}
  \item{vp}{Fractional occupancy by plasma.}
  \item{ktranserror}{Standard error for ktrans.}
  \item{keperror}{Standard error for kep.}
  \item{vperror}{Standard error for vp.}
  The residual sum-of-squares is also provided, along with the original
  acquisition times (for plotting purposes).
}
\references{
  Ahearn, T.S., Staff, R.T., Redpath, T.W. and Semple, S.I.K. (2005) The
  use of the Levenburg-Marquardt curve-fitting algorithm in
  pharmacokinetic modelling of DCE-MRI data, \emph{Physics in Medicine
    and Biology}, \bold{50}, N85-N92.

Fritz-Hansen, T., Rostrup, E., Larsson, H.B.W, Sondergaard, L., Ring,
  P. and Henriksen, O. (1993) Measurement of the arterial concentration
  of Gd-DTPA using MRI: A step toward quantitative perfusion imaging,
  \emph{Magnetic Resonance in Medicine}, \bold{36}, 225-231.

  Orton, M.R., Collins, D.J., Walker-Samuel, S., d'Arcy, J.A.,
  Hawkes, D.J., Atkinson, D. and Leach, M.O. (2007) Bayesian
  estimation of pharmacokinetic parameters for DCE-MRI with a robust
  treatment of enhancement onset time, \emph{Physics in Medicine and
  Biology} \bold{52}, 2393-2408.
  
  Orton, M.R., d'Arcy, J.A., Walker-Samuel, S., Hawkes, D.J.,
  Atkinson, D., Collins, D.J. and Leach, M.O. (2008) Computationally
  efficient vascular input function models for quantitative kinetic
  modelling using DCE-MRI, \emph{Physics in Medicine and Biology}
  \bold{53}, 1225-1239.

  Tofts, P.S., Brix, G, Buckley, D.L., Evelhoch, J.L., Henderson, E.,
  Knopp, M.V., Larsson, H.B.W., Lee, T.-Y., Mayr, N.A., Parker, G.J.M.,
  Port, R.E., Taylor, J. and Weiskoff, R. (1999) Estimating kinetic
  parameters from dynamic contrast-enhanced \eqn{T_1}-weighted MRI of a
  diffusable tracer: Standardized quantities and symbols, \emph{Journal
    of Magnetic Resonance}, \bold{10}, 223-232.

  Tofts, P.S. and Kermode, A.G. (1984) Measurement of the blood-brain
  barrier permeability and leakage space using dynamic MR
  imaging. 1. Fundamental concepts, \emph{Magnetic Resonance in
    Medicine}, \bold{17}, 357-367.

  Weinmann, H.J., Laniado, M. and Mutzel, W. (1984) Pharmacokinetics of
  Gd-DTPA/dimeglumine after intraveneous injection into healthy
  volunteers, \emph{Physiological Chemistry and Physics and Medical
    NMR}, \bold{16}, 167-172. 
}
\seealso{
  \code{\link{dcemri.bayes}}, \code{\link{dcemri.map}},
  \code{\link{dcemri.spline}}, \code{\link[minpack.lm]{nls.lm}}
}
\examples{
data("buckley")

## Set up breast data for dcemri
xi <- seq(5, 300, by=5)
img <- array(t(breast$data)[,xi], c(13,1,1,60))
time <- buckley$time.min[xi]
aif <- buckley$input[xi]
mask <- array(TRUE, dim(img)[1:3])

## Generate AIF params using the orton.exp function from Buckley's AIF
aifparams <- orton.exp.lm(time, aif)
fit <- dcemri.lm(img, time, mask, model="orton.exp",
                 aif="user", user=aifparams)

## Scatterplot comparing true and estimated Ktrans values
plot(breast$ktrans, fit$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans})))
abline(0, 1, lwd=1.5, col="red")
cbind(breast$ktrans, fit$ktrans[,,1])
}
\author{Brandon Whitcher, Volker Schmid}
\keyword{models}
