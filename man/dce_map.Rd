% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: dce_map.Rd 191 2009-08-25 15:12:31Z bjw34032 $

\name{dcemri.map}
\alias{dcemri.map}
\alias{dcemri.map,array-method}
\alias{dcemri.map.single}
\title{Pharmacokinetic Modeling of Dynamic Contrast-Enhanced MRI Data}
\description{
  Maximum-a-posteriori (MAP) estimation for single compartment models is
  performed using literature-based or user-specified arterial input
  functions.
}
\usage{
\S4method{dcemri.map}{array}(conc, time, img.mask, model="extended", aif=NULL,
             user=NULL, tau.ktrans=1, tau.kep=tau.ktrans,
             ab.vp=c(1,19), ab.tauepsilon=c(1,1/1000),
             samples=FALSE, multicore=FALSE, verbose=FALSE, ...)
dcemri.map.single(conc, time, posterior, parameter, transform, start, hyper, aif)
}
\arguments{
  \item{conc}{Matrix or array of concentration time series (last
    dimension must be time).}
  \item{time}{Time in minutes.}
  \item{img.mask}{Mask matrix or array. Voxels with \code{mask=0} will be
    excluded.}
  \item{model}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \itemize{
      \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with single
	compartment model}
      \item{\dQuote{extended}}{Weinmann model extended with additional
	vascular compartment (default)}
    }
  }
  \item{aif}{is a character string that identifies the parameters of the
    type of arterial input function (AIF) used with the above model.
    Acceptable values are: \code{tofts.kermode} (default) or
    \code{fritz.hansen} for the \code{weinmann} and \code{extended}
    models; \code{orton.exp} (default) or \code{user} for the
    \code{orton.exp} model.}
  \item{user}{Vector of AIF parameters.  For Tofts and Kermode:
    \eqn{a_1}, \eqn{m_1}, \eqn{a_2}, \eqn{m_2}; for Orton \emph{et al.}:
    \eqn{A_b}, \eqn{\mu_b}, \eqn{A_g}, \eqn{\mu_g}.}
  \item{tau.ktrans}{Variance parameter for prior on
    \eqn{\log(K^{trans})}.}
  \item{tau.kep}{Variance parameter for prior on \eqn{\log(k_{ep})}.}
  \item{ab.vp}{Hyper-prior parameters for the Beta prior on
    \eqn{v_p}{vp}.}
  \item{ab.tauepsilon}{Hyper-prior parameters for observation error
    Gamma prior.} 
  \item{samples}{If \code{TRUE} output includes samples drawn
    from the posterior distribution for all parameters.}
  \item{multicore}{If \code{TRUE} algorithm is parallelized using
    \pkg{multicore}.}
  \item{verbose}{}
  \item{...}{}
  \item{posterior}{}
  \item{parameter}{}
  \item{transform}{}
  \item{start}{}
  \item{hyper}{}
  aif
}
\details{
  Implements maximum-a-posteriori (MAP) estimation for the Bayesian
  model in Schmid \emph{et al.} (2006).
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
  \item{sigma2}{The residual sum-of-squares from the model fit.}
  \item{time}{Acquisition times (for plotting purposes).}
  Note, not all parameters are available under all models choices.
}
\references{
  Schmid, V., Whitcher, B., Padhani, A.R., Taylor, N.J. and Yang, G.-Z. 
  (2006) Bayesian methods for pharmacokinetic models in dynamic
  contrast-enhanced magnetic resonance imaging, \emph{IEEE Transactions
    on Medical Imaging}, \bold{25} (12), 1627-1636.
}
\seealso{
  \code{\link{dcemri.lm}}, \code{\link{dcemri.bayes}}
}
\examples{
data("buckley")
xi <- seq(5, 300, by=5)
img <- array(t(breast$data)[,xi], c(13,1,1,60))
mask <- array(TRUE, dim(img)[1:3])
time <- buckley$time.min[xi]

## MAP estimation with Fritz-Hansen default AIF
fit.map <- dcemri.map(img, time, mask, aif="fritz.hansen",
                      nriters=5000)

plot(breast$ktrans, fit.map$ktrans, xlim=c(0,1), ylim=c(0,1),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans}, " (MAP)")))
abline(0, 1, lwd=1.5, col=2)

\dontrun{
fit.lm <- dcemri.lm(img, time, mask, aif="fritz.hansen")

plot(breast$ktrans, fit.map$ktrans, xlim=c(0,1), ylim=c(0,1),
     xlab=expression(paste("True ", K^{trans})),
     ylab=expression(paste("Estimated ", K^{trans})))
points(breast$ktrans, fit.lm$ktrans, pch=3)
abline(0, 1, lwd=1.5, col="red")
legend("bottomright", c("MAP Estimation (fritz-hansen)",
                        "Levenburg-Marquardt (fritz.hansen)", pch=c(1,3))
}
}
\author{Volker Schmid <\email{volkerschmid@users.sourceforge.net}>}
\keyword{models}
