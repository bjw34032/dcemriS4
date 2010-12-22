\name{compartmentalModel}
\alias{compartmentalModel}
\title{
  Compartmental Models for Kinetic Parameter Estimation
}
\description{
  A selection of parametric models are provided that combine a
  compartmental model for tissue and a functional form of the arterial
  input function.
}
\usage{
compartmentalModel(type)
}
\arguments{
  \item{type}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \describe{
      \item{\dQuote{weinmann}}{Weinmann AIF convolved with a
	single compartment (Kety) model}
      \item{\dQuote{extended}}{Kety model extended with additional
	vascular compartment (default)}
      \item{\dQuote{orton.exp}}{Extended model using Orton's exponential
	arterial input function} 
      \item{\dQuote{orton.cos}}{Extended model using Orton's raised
	cosine arterial input function}
      \item{\dQuote{kety.orton.exp}}{Kety model using Orton's
	exponential arterial input function}
      \item{\dQuote{kety.orton.cos}}{Kety model using Orton's raised
	cosine arterial input function}
      \item{\dQuote{weinmann.empirical}}{User-specified empirical AIF
	convolved with a single compartment model}
      \item{\dQuote{extended.empirical}}{Extended model using a
	user-specified empirical arterial input function}
    }
  }
}
\details{
  Parametric models from the DCE-MRI literature are provided to the
  user for kinetic parameter estimation.  All models, with the exception
  of those marked \sQuote{empirical} incorporate a parametric model for
  the arterial input function (AIF).
}
\value{
  A function.
}
%\references{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\seealso{
  \code{\link{aifParameters}}, \code{\link{dcemri.bayes}},
  \code{\link{dcemri.lm}}, \code{\link{dcemri.map}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
