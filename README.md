# dcemriS4

##Open-source software has arrived!

**dcemriS4** is a package, written in the <a href="http://www.r-project.org/">R programming environment</a>, for the quantitative analysis of dynamic contrast-enhanced MRI (DCE-MRI) and diffusion-weighted imaging (DWI) for oncology applications. It has been released under the <a href="http://www.opensource.org/licenses/bsd-license.php">BSD license</a>.

The scientific backbone of this software is based on research in the area of parameter estimation and statistical inference over several years. Three scientific publications outline the major techniques available:

* Schmid, Whitcher, et al. (2006), <a href="http://dx.doi.org/10.1109/TMI.2006.884210">Bayesian Methods for Pharmacokinetic Models in Dynamic Contrast-Enhanced Magnetic Resonance Imaging</a>, *IEEE Transactions in Medical Imaging*, **25** (12), 1627-1636.
* Schmid, Whitcher, et al. (2009), <a href="http://dx.doi.org/10.1002/mrm.21807">A Bayesian Hierarchical Model for the Analysis of a Longitudinal Dynamic Contrast-Enhanced MRI Oncology Study</a>, *Magnetic Resonance in Medicine*, **61** (1) 163-174.
* Schmid, Whitcher, et al. (2009), <a href="http://dx.doi.org/10.1109/TMI.2008.2007326">A Semi-parametric Technique for the Quantitative Analysis of Dynamic Contrast-enhanced MR Images Based on Bayesian P-splines</a>, *IEEE Transactions in Medical Imaging*, **28** (6) 789-798. 

An application of both the Bayesian and Frequentist methods to fit compartmental models to DCE-MRI data from a clinical trial may be found at

* Whitcher, Schmid, et al. (2011), <a href="http://dx.doi.org/10.1007/s10334-010-0238-3">A Bayesian Hierarchical Model for DCE-MRI to Evaluate Treatment Response in a Phase II Study in Advanced Squamous Cell Carcinoma of the Head and Neck</a>, *Magnetic Resonance Materials in Physics, Biology and Medicine*, **24** (2), 85-96.

The **dcemriS4** package provides a comprehensive set of R functions (subroutines) that perform all the necessary tasks to produce quantitative estimates of tumor perfusion/permeability using "standard" kinetic models. Data must be in one of two standard formats,	<a href="http://www.mayo.edu/bir/PDF/ANALYZE75.pdf">Analyze</a> or <a href="http://nifti.nimh.nih.gov/">NIfTI</a>, and results may be written out to either format.  If your data are in <a href="http://medical.nema.org">DICOM</a> format, then we suggest using <a href="https://github.com/bjw34032/oro.dicom">oro.dicom</a> in R or feel free to use one of the many OSS implementations that convert to Analyze or NIfTI.

Travis-CI build status for bjw34032/dcemriS4

[![Travis-CI Build Status](https://travis-ci.org/bjw34032/dcemriS4.svg?branch=master)](https://travis-ci.org/bjw34032/dcemriS4)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/dcemriS4)](http://cran.rstudio.com/web/packages/dcemriS4/index.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/dcemriS4)](http://cran.rstudio.com/web/packages/dcemriS4/index.html)
