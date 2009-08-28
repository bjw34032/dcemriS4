% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
% 
% $Id: buckley.Rd 191 2009-08-25 15:12:31Z bjw34032 $

\name{Buckley}
\alias{buckley}
\alias{breast}
\alias{meningioma}
\title{Simulated Data from Buckley (2002)}
\description{
  In Buckley (2002) tissue residue curves for a Meningioma and a Breast
  Cancer were simulated using the MMID4 model. Note, the model is
  described in detail by Bassingthwaighte, J.B. \emph{et al.} (1984) and
  Kroll, K \emph{et al.} (1996).  This model accounts for flow
  dispersion and heterogeneity, and includes capillaries modeled as
  axially distributed blood-tissue exchange units.  A plasma
  concentration-time curve, AKA arterial input function, was simulated
  as an input to the model using measurements made by Fritz-Hansen
  \emph{et al.} (1996).
}
\usage{
data("buckley")
}
\format{
  Two lists are created (breast and meningioma) that contain the
  simulated time curves and all associated kinetic parameter values.
}
\source{
  See below.
}
\references{
  Buckley, D.L. (2002) Uncertainty in the Analysis of Tracer Kinetics
  Using Dynamic Contrast-Enhanced \eqn{T_1}{T1}-weighted MRI,
  \emph{Magnetic Resonance in Medicine} \bold{47}, 601-606.

  Bassingthwaighte, J.B. and Goresky, C.A. (1984) Modelling in the analysis of
  solute and water exchange in the microvasculature. In: Renkin, E.M., Michel,
  C.C. and Geiger, S.R., editors. Handbook of physiology. Section 2. The
  cardiovascular system. Bethesda: American Physiological Society. p549-626.

  Kroll, K., Wilke, N., Jerosch-Herold, M., Wang, Y., Zhang Y., Basche,
  R.J. and Bassingthwaighte, J.B. (1996) Modelling regional myocardial
  flows from residue functions of an intravascular indicator. \emph{Am J
    Physiol} \bold{271}, H1643-H1655.

  Fritz-Hansen, T., Rostrup, E., Larsson, H.B., Sondergaard, L., Ring, P. and
  Hendriksen, O. (1996) Measurement of the arterial concentration Gd-DTPA using
  MRI; a step toward quantitative perfusion imaging. \emph{Magn Reson Med}
  \bold{36}, 347-357.
}
\keyword{datasets}
