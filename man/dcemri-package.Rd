\name{dcemri-package}
\Rdversion{1.1}
\alias{dcemri-package}
\alias{dcemri}
\docType{package}
\title{dcemri: A Package for Medical Image Analysis}
\description{
  A collection of routines and documentation that allows one to
  perform a quantitative analysis of dynamic contrast-enhanced or
  diffusion-weighted MRI data.  Medical imaging data should be organized
  using either the Analyze or NIfTI data formats.
}
\details{
  \tabular{ll}{
    Package: \tab dcemri\cr
    Version: \tab 0.1.4\cr
    Date: \tab 2009-08-25\cr
    Depends: \tab R (>= 2.6.0), graphics, grDevices\cr
    Suggests: \tab bitops, minpack.lm, multicore, splines\cr
    License: \tab BSD\cr
    URL: \tab http://www.dcemri.org/, http://dcemri.sourceforge.net/\cr
    Packaged: \tab 2009-08-25 15:14:01 UTC; bjw34032\cr
  }

  Index:
  \preformatted{
    R10.lm                Estimate Intrinsic Tissue Relaxivity
    adc.lm                Estimate the Apparent Diffusion
                          Coefficient (ADC)
    aif.orton.exp         Arterial Input Functions
    buckley               Simulated Data from Buckley (2002)
    conv.fft              Convolution of 3D Arrays using the
                          Fourier Transform
    dam                   Double-Angle Method for B1+ Mapping
    dcemri.bayes          Bayesian Methods for Pharmacokinetic
                          Modeling of Dynamic Contrast-Enhanced
                          MRI Data
    dcemri.lm             Pharmacokinetic Models for Dynamic
			  Contrast-Enhanced MRI Data
    dcemri.map            Pharmacokinetic Modeling of Dynamic
			  Contrast-Enhanced MRI Data
    dcemri.spline         Bayesian P-Splines for Dynamic Contrast-
                          Enhanced MRI Data
    extract.aif           Seed Growing for a 4D Array
    find.center           Find the Center of a Binary Mask
    ftm                   Fast Template Matching via
                          Cross-Correlation
    hotmetal              Hot Metal Color Table
    kineticModel          Pharmacokinetic Models
    nsli                  Dimension Accessor Functions
    read.analyze.hdr      Read Analyze Format Headers and Binary
    Image Files
    read.hdr              Read Headers and Binary Image Files for
                          Medical Imaging Data
    read.nifti.hdr        Read NIfTI Format Headers and Image
			  Files
    shift3D               Shift a 3D Array in One Dimension
    tim.colors            Tim's Useful Color Table
    write.analyze.img     Write Analyze/NIfTI Format Binary
                          Image Files
  }

  Further information is available in the following vignettes:
  \tabular{ll}{
    \code{dcemri} \tab dcemri(source, pdf)\cr
  }
}
\author{
  Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>\cr
  Volker Schmid <\email{volkerschmid@users.sourceforge.net}>\cr
  Andrew Thornton <\email{zeripath@users.sourceforge.net}> 

  Maintainer: Brandon Whitcher <bjw34032@users.sourceforge.net>
}
\references{
  Schmid, V., Whitcher, B., Padhani, A.R., Taylor, N.J. and Yang, G.-Z. 
  (2006) Bayesian methods for pharmacokinetic models in dynamic
  contrast-enhanced magnetic resonance imaging, \emph{IEEE Transactions
    on Medical Imaging}, \bold{25} (12), 1627-1636.

  Schmid, V., Whitcher, B., Padhani, A.R. and G.-Z. Yang (2009) A
  semi-parametric technique for the quantitative analysis of dynamic
  contrast-enhanced MR images based on Bayesian P-splines, \emph{IEEE
    Transactions on Medical Imaging}, \bold{28} (6), 789-798.
}
\keyword{package}
%\seealso{}
\examples{
\dontrun{
demo(avg152T1)
demo(avg152T1_LR)
demo(avg152T1_RL)
demo(buckley)
demo(filtered_func_data)
demo(zstat1)
}
}
