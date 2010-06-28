\name{dcemriS4-package}
\Rdversion{1.1}
\alias{dcemriS4-package}
\alias{dcemriS4}
\docType{package}
\title{dcemri: A Package for Medical Image Analysis (S4 implementation)}
\description{
  A collection of routines and documentation that allows one to
  perform a quantitative analysis of dynamic contrast-enhanced or
  diffusion-weighted MRI data.  Medical imaging data should be organized
  using either the Analyze or NIfTI data formats.
}
\details{
  \tabular{ll}{
    Package: \tab dcemriS4\cr
    Version: \tab 0.20.1\cr
    Date: \tab 2009-09-02\cr
    Depends: \tab R (>= 2.6.0), grDevices, graphics, methods, utils\cr
    Suggests: \tab bitops, minpack.lm, multicore, splines, XML\cr
    License: \tab BSD\cr
    URL: \tab http://www.dcemri.org/, http://dcemri.sourceforge.net/\cr
  }
  
  Index:
  \preformatted{
    R10.lm                Estimate Intrinsic Tissue Relaxivity
    [<--methods           ~~ Methods for Function [<- in Package 'base' ~~
    adc.lm                Estimate the Apparent Diffusion Coefficient
                          (ADC)
    aif.orton.exp         Arterial Input Functions
    anlz                  Constructor for Analyze
    anlz-class            Class "anlz"
    as.nifti              as.nifti
    audit.trail-methods   Extract or Replace NIfTI Audit Trail
    aux.file-methods      Extract or Replace NIfTI/Analyze Auxiliary
    File
    buckley               Simulated Data from Buckley (2002)
    coerce,array,nifti-method Force an Object to Belong to the NIfTI
                              Class
    conv.fft              Convolution of 3D Arrays using the Fourier
                          Transform
    convert.datatype      Convert NIfTI Codes
    convert.datatype.anlz Convert NIfTI Codes
    dam                   Double-Angle Method for B1+ Mapping
    dcemri.bayes          Bayesian Methods for Pharmacokinetic
                          Modeling of Dynamic Contrast-Enhanced MRI
                          Data
    dcemri.info           Facilitate the Creation and Modification of
                          Audit Trails
    dcemri.lm             Pharmacokinetic Models for Dynamic
                          Contrast-Enhanced MRI Data
    dcemri.map            Pharmacokinetic Modeling of Dynamic
                          Contrast-Enhanced MRI Data
    dcemri.spline         Bayesian P-Splines for Dynamic
                          Contrasat-Enhanced MRI Data
    descrip-methods       Extract or Replace NIfTI/Analyze Description
    extract.aif           Seed Growing for a 4D Array
    find.center           Find the Center of a Binary Mask
    ftm                   Fast Template Matching via Cross-Correlation
    hotmetal              Hot Metal Color Table
    image-methods         Methods for Function 'image' in Package
                          'dcemri'
    integerTranslation    integerTranslation
    kineticModel          Pharmacokinetic Models
    nifti                 Constructor for NIfTI
    nifti-class           Class "nifti"
    niftiAuditTrail-class Class "niftiAuditTrail"
    niftiExtension-class  Class "niftiExtension"
    niftiExtensionSection-class Class "niftiExtensionSection" ~~~
    nsli                  Dimension Accessor Functions
    orthographic-methods  ~~ Methods for Function orthographic ~~
    overlay-methods       ~~ Methods for Function overlay ~~
    performPermutation    performPermutation
    quaternion2rotation   Convert Quaternion into a Rotation Matrix
    readANALYZE           readANALYZE
    readNIfTI             readNIfTI
    reorient              reorient
    shift3D               Shift a 3D Array in One Dimension
    tim.colors            Tim's Useful Color Table
    translateCoordinate   translateCoordinate
    writeANALYZE          writeANALYZE
    writeNIfTI            writeNIfTI
    xyzt2space            Bitwise Conversion Subroutines
  }

  Further information is available in the following vignettes:
  \tabular{ll}{
    \code{dcemriS4} \tab dcemriS4(source, pdf)\cr
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
demo(avg152T1LR)
demo(avg152T1RL)
demo(buckley)
demo("filtered_func_data")
demo(zstat1)
}
}
