% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dce_bayes.Rd 332 2010-01-29 16:54:07Z bjw34032 $

\name{dcemri.bayes}
\alias{dcemri.bayes}
\alias{dcemri.bayes,array-method}
\alias{dcemri.bayes.single}
\title{Bayesian Methods for Pharmacokinetic Modeling of Dynamic
  Contrast-Enhanced MRI Data}
\description{
  Bayesian analysis of contrast agent concentration time curves from
  DCE-MRI.
}
\usage{
\S4method{dcemri.bayes}{array}(conc, time, img.mask, model="extended",
             aif=NULL, user=NULL, nriters=3500, thin=3, burnin=500,
             tune=267, ab.ktrans=c(0,1), ab.kep=ab.ktrans,
             ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000),
             samples=FALSE, multicore=FALSE, verbose=FALSE, dic=FALSE,
             ...)
dcemri.bayes.single(conc, time, nriters=3500, thin=3,
                    burnin=500, tune=267, ab.gamma=c(0,1),
                    ab.theta=c(0,1), ab.vp=c(1,19),
                    ab.tauepsilon=c(1,1/1000), aif.model=0,
                    aif.parameter=c(2.4,0.62,3,0.016), vp=1)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series (last
    dimension must be time).}
  \item{time}{Time in minutes.}
  \item{img.mask}{Mask matrix or array. Voxels with mask=0 will be
    excluded.}
  \item{model}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \itemize{
      \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with single 
	compartment model}
      \item{\dQuote{extended}}{Weinmann model extended with additional 
	vascular compartment (default)}
      \item{\dQuote{orton.exp}}{Extended model using Orton's exponential 
	AIF}
      \item{\dQuote{kety.orton.exp}}{Kety model using Orton's
	exponential AIF}
    }
  }
  \item{aif}{is a character string that identifies the parameters of the
    type of arterial input function (AIF) used with the above model.
    Acceptable values are: \code{tofts.kermode} (default) or
    \code{fritz.hansen} for the \code{weinmann} and \code{extended}
    models; \code{orton.exp} (default) or \code{user} for the
     \code{orton.exp} and \code{kety.orton.exp} model.}
  \item{user}{Vector of AIF parameters.  For Tofts and Kermode:
    \eqn{a_1}{a1}, \eqn{m_1}{m1}, \eqn{a_2}{a2}, \eqn{m_2}{m2}; for
    Orton \emph{et al.}: \eqn{A_b}{Ab}, \eqn{\mu_b}{mub}, \eqn{A_g}{Ag},
    \eqn{\mu_g}{mug}.} 
  \item{nriters}{Total number of iterations.}
  \item{thin}{Thining factor.}
  \item{burnin}{Number of iterations for burn-in.}
  \item{tune}{Number for iterations for tuning.  The algorithm will be 
    tuned to an acceptance rate between 0.3 and 0.6.}
  \item{ab.ktrans}{Mean and variance parameter for Gaussian prior on
    \eqn{\log(K^{trans})}{log(Ktrans)}.}
  \item{ab.kep}{Mean and variance parameter for Gaussian prior on
    \eqn{\log(k_{ep})}{log(kep)}.}
  \item{ab.vp}{Hyper-prior parameters for the Beta prior on
    \eqn{v_p}{vp}.}
  \item{ab.gamma}{\ldots}
  \item{ab.theta}{\ldots}
  \item{ab.tauepsilon}{Hyper-prior parameters for observation error
    Gamma prior.}
  \item{samples}{If \code{TRUE} output includes samples drawn
    from the posterior distribution for all parameters.}
  \item{multicore}{If \code{TRUE} algorithm is parallelized using
    \pkg{multicore}.}
  \item{verbose}{Logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
  \item{dic}{If \code{TRUE}, the deviance information criterion (DIC)
    and effective number of parameters (pD) will be computed.  If
    \dQuote{samples = TRUE}, then samples of the DIC and pD will be
    given.}
  \item{aif.model}{\ldots}
  \item{aif.parameter}{\ldots}
  \item{vp}{\ldots}
  \item{...}{Additional parameters to the function.}
}
\details{
  See Schmid \emph{et al.} (2006) for details.
}
\value{
  Parameter estimates and their standard errors are provided for the
  masked region of the multidimensional array.  All multi-dimensional
  arrays are output in \code{nifti} format.

  They include:
  \item{ktrans}{Transfer rate from plasma to the extracellular,
    extravascular space (EES).}
  \item{ktranserror}{Error on \eqn{K^{trans}}{ktrans}.}
  \item{kep}{Rate parameter for transport from the EES to plasma.}
  \item{keperror}{Error on \eqn{k_{ep}}{kep}.}
  \item{ve}{Fractional occupancy by EES (the ratio between ktrans and
    kep).}
  \item{vperror}{Error on \eqn{v_e}{ve}.}
  \item{vp}{Fractional occupancy by plasma.}
  \item{sigma2}{The residual sum-of-squares from the model fit.}
  \item{time}{Acquisition times (for plotting purposes).}
  \item{DIC}{Deviance information criterion.}
  \item{DIC.map}{Contribution to DIC per voxel.}
  \item{pD}{Effective number of parameters.}
  \item{pD.map}{Constribution to pD per voxel.}
  Note, not all parameters are available under all models choices.
}
\references{
  Schmid, V., Whitcher, B., Padhani, A.R., Taylor, N.J. and Yang, G.-Z. 
  (2006) Bayesian methods for pharmacokinetic models in dynamic
  contrast-enhanced magnetic resonance imaging, \emph{IEEE Transactions
    on Medical Imaging}, \bold{25} (12), 1627-1636.
}
\seealso{
  \code{\link{dcemri.lm}}, \code{\link{dcemri.map}},
  \code{\link{dcemri.spline}}
}
\examples{
data("buckley")
xi <- seq(5, 300, by=5)
img <- array(t(breast$data)[,xi], c(13,1,1,60))
mask <- array(TRUE, dim(img)[1:3])
time <- buckley$time.min[xi]

## Bayesian estimation with Fritz-Hansen default AIF
fit.bayes <- dcemri.bayes(img, time, mask, aif="fritz.hansen")

## Bayesian estimation with "orton.exp" function fit to Buckley's AIF
aif <- buckley$input[xi]
aifparams <- orton.exp.lm(time, aif)
aifparams$D <- 1
fit.bayes.aif <- dcemri.bayes(img, time, mask, model="orton.exp",
                              aif="user", user=aifparams)

plot(breast$ktrans, fit.bayes$ktrans, xlim=c(0,1), ylim=c(0,1),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans}, " (Bayesian)")))
points(breast$ktrans, fit.bayes.aif$ktrans, pch=2)
abline(0, 1, lwd=2, col=2)
legend("right", c("extended/fritz.hansen","orton.exp/user"), pch=1:2)
cbind(breast$ktrans, fit.bayes$ktrans[,,1], fit.bayes.aif$ktrans[,,1])

fit.lm <- dcemri.lm(img, time, mask, aif="fritz.hansen")
fit.lm.aif <- dcemri.lm(img, time, mask, model="orton.exp", aif="user",
                        user=aifparams)

plot(breast$ktrans, fit.bayes$ktrans, xlim=c(0,1), ylim=c(0,1),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans})))
points(breast$ktrans, fit.bayes.aif$ktrans, pch=2)
points(breast$ktrans, fit.lm$ktrans, pch=3)
points(breast$ktrans, fit.lm.aif$ktrans, pch=4)
abline(0, 1, lwd=2, col=2)
legend("bottomright", c("Bayesian Estimation (fritz-hansen)",
                        "Bayesian Estimation (orton.exp)",
                        "Levenburg-Marquardt (weinmann/fritz.hansen)",
                        "Levenburg-Marquardt (orton.exp/user)"), pch=1:4)
}
\author{Volker Schmid \email{volkerschmid@users.sourceforge.net}}
\keyword{models}
