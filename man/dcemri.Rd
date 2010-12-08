% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dcemri.Rd 336 2010-02-19 17:26:47Z volkerschmid $

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
          control=nls.lm.control(), user=NULL, guess=NULL, multicore=FALSE,
          verbose=FALSE, ...)
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
    \describe{
      \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with single
	compartment model}
      \item{\dQuote{extended}}{Weinmann model extended with additional
	vascular compartment (default)}
      \item{\dQuote{orton.exp}}{Extended model using Orton's exponential
	AIF} 
      \item{\dQuote{orton.cos}}{Extended model using Orton's raised
	cosine AIF}
      \item{\dQuote{kety.orton.exp}}{Kety model using Orton's
	exponential AIF}
      \item{\dQuote{kety.orton.cos}}{Kety model using Orton's raised
	cosine AIF}
    }
  }
  \item{aif}{is a character string that identifies the parameters of the
    type of arterial input function (AIF) used with the above model.
    Acceptable values are:
    \itemize{
      \item{\code{tofts.kermode}}{(default) for the \code{weinmann} and
	\code{extended}	models}
      \item{\code{fritz.hansen}}{for the \code{weinmann} and
	\code{extended} models}
      \item{\dQuote{empirical}}{for the \code{weinmann} and
	\code{extended} models}
      \item{\code{orton.exp}}{(default) for the \code{orton.exp} and
	\code{kety.orton.exp} model} 
      \item{\code{orton.cos}}{(default) for the	\code{orton.cos} and
	\code{kety.orton.cos} model.}
      \item{\code{user}}{for the \code{orton.exp} and \code{orton.cos}
	model.}
    }
    All AIF models set the parametric form and parameter values --
    except \code{user}, where a set of user-defined parameter values are
    allowed, and \code{empirical}, where a vector of values that fully
    characterize the empirical AIF.
  }
  \item{control}{is a list of parameters used by \code{nls.lm.control}
    that are set by default, but may be customized by the user.}
  \item{user}{is a list with the following parameters required: D, AB,
    muB, AG, muG.}
  \item{guess}{is a vector of starting values for kinetic parameter
    estimation.  The vector must have length = 3 (with names \code{th0},
    \code{th1} and \code{th3}) when the extended Kety model is used with
    the vascular parameter and length = 2 (with names \code{th1} and
    \code{th3}) otherwise.}
  \item{multicore}{is a logical variable (default = \code{FALSE}) that
    allows parallel processing via \pkg{multicore}.}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
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
  masked region of the multidimensional array.  All multi-dimensional
  arrays are provided in \code{nifti} format.

  They include:
  \item{ktrans}{Transfer rate from plasma to the extracellular,
    extravascular space (EES).}
  \item{kep}{Rate parameter for transport from the EES to plasma.}
  \item{ve}{Fractional occupancy by EES (the ratio between
  \eqn{K^{trans}}{Ktrans} and \eqn{k_{ep}}{kep}).}
  \item{vp}{Fractional occupancy by plasma.}
  \item{ktranserror}{Standard error for \eqn{K^{trans}}{Ktrans}.}
  \item{keperror}{Standard error for \eqn{k_{ep}}{kep}.}
  \item{vperror}{Standard error for \eqn{v_p}{vp}.}
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

## Empirical arterial input function
img <- array(t(breast$data), c(13,1,1,301))
time <- buckley$time.min
mask <- array(TRUE, dim(img)[1:3])

## Estimate kinetic parameters directly from Buckley's empirical AIF
fit1 <- dcemri.lm(img, time, mask, model="weinmann", aif="empirical",
                  user=buckley$input)
fit2 <- dcemri.lm(img, time, mask, model="extended", aif="empirical",
                  user=buckley$input)

## Set up breast data for dcemri
xi <- seq(5, 300, by=3)
img <- array(t(breast$data)[,xi], c(13,1,1,length(xi)))
time <- buckley$time.min[xi]
input <- buckley$input[xi]

## Generate AIF params using the orton.exp function from Buckley's AIF
(aifparams <- orton.exp.lm(time, input))
fit3 <- dcemri.lm(img, time, mask, model="orton.exp", aif="user",
                  user=aifparams)

## Scatterplot comparing true and estimated Ktrans values
plot(breast$ktrans, fit1$ktrans, xlim=c(0,0.75), ylim=c(0,0.75),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans})))
points(breast$ktrans, fit2$ktrans, pch=2)
points(breast$ktrans, fit3$ktrans, pch=3)
abline(0, 1, lwd=1.5, col=2)
legend("bottomright", c("weinmann/empirical", "extended/empirical",
                        "orton.exp/user"), pch=1:3)
cbind(breast$ktrans, fit1$ktrans[,,1], fit2$ktrans[,,1], fit3$ktrans[,,1])

## Scatterplot comparing true and estimated Ktrans values
plot(breast$vp, fit1$vp, type="n", xlim=c(0,0.15), ylim=c(0,0.15),
     xlab=expression(paste("True ", v[p])),
     ylab=expression(paste("Estimated ", v[p])))
points(breast$vp, fit2$vp, pch=2)
points(breast$vp, fit3$vp, pch=3)
abline(0, 1, lwd=1.5, col=2)
legend("bottomright", c("extended/empirical","orton.exp/user"), pch=2:3)
cbind(breast$vp, fit2$vp[,,1], fit3$vp[,,1])
}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net},
  Volker Schmid \email{volkerschmid@users.sourceforge.net}}
\keyword{models}
