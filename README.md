# A Package for the Quantitative Analysis of DCE-MRI

![Travis-CI Build Status](https://travis-ci.org/bjw34032/dcemriS4.svg?branch=master)
![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/dcemriS4)
![CRAN Downloads Badge](http://cranlogs.r-pkg.org/badges/grand-total/dcemriS4)
![CRAN Downloads Badge](http://cranlogs.r-pkg.org/badges/dcemriS4)

**dcemriS4** is a package, written in the [R](https://www.r-project.org) programming environment, for the quantitative analysis of dynamic contrast-enhanced MRI (DCE-MRI) and diffusion-weighted imaging (DWI) for oncology applications. It has been released under the [BSD license](https://www.opensource.org/licenses/bsd-license.php).

## Getting started

Install the R package using the following commands on the R console

```
install.packages("devtools")
devtools::install_github("bjw34032/dcemriS4")
library(dcemriS4)
```

or using [NeuroConductor](https://neuroconductor.org)

```
source("https://neuroconductor.org/neurocLite.R")

# Default Install
neuro_install('dcemriS4')

# from GitHub
neuro_install('dcemriS4', release = "stable", release_repo = "github")
neuro_install('dcemriS4', release = "current", release_repo = "github")
```

At this point in time the package is not available on CRAN.  

## References

The scientific backbone of this software is based on research in the area of parameter estimation and statistical inference.  Three scientific publications outline the major techniques available:

* Schmid, Whitcher, et al. (2006), [Bayesian Methods for Pharmacokinetic Models in Dynamic Contrast-Enhanced Magnetic Resonance Imaging](https://dx.doi.org/10.1109/TMI.2006.884210), *IEEE Transactions in Medical Imaging*, **25** (12), 1627-1636.
* Schmid, Whitcher, et al. (2009), [A Bayesian Hierarchical Model for the Analysis of a Longitudinal Dynamic Contrast-Enhanced MRI Oncology Study](https://dx.doi.org/10.1002/mrm.21807), *Magnetic Resonance in Medicine*, **61** (1) 163-174.
* Schmid, Whitcher, et al. (2009), [A Semi-parametric Technique for the Quantitative Analysis of Dynamic Contrast-enhanced MR Images Based on Bayesian P-splines](https://dx.doi.org/10.1109/TMI.2008.2007326), *IEEE Transactions in Medical Imaging*, **28** (6) 789-798. 

An application of both the Bayesian and frequentist methods to fit compartmental models to DCE-MRI data from a clinical trial may be found at

* Whitcher, Schmid, et al. (2011), <a href="http://dx.doi.org/10.1007/s10334-010-0238-3">A Bayesian Hierarchical Model for DCE-MRI to Evaluate Treatment Response in a Phase II Study in Advanced Squamous Cell Carcinoma of the Head and Neck</a>, *Magnetic Resonance Materials in Physics, Biology and Medicine*, **24** (2), 85-96.

The **dcemriS4** package provides a comprehensive set of R functions that perform all the necessary tasks to produce quantitative estimates of tumor perfusion/permeability using "standard" kinetic models. Data must be in one of two standard formats, Analyze or [NIfTI](http://nifti.nimh.nih.gov), and results may be written out to either format.  If your data are in [DICOM](https://medical.nema.org) format, we suggest using [**oro.dicom**](https://github.com/bjw34032/oro.dicom) in R or feel free to use one of the many OSS implementations that convert to Analyze or NIfTI.
