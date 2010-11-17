\name{compartmentalModel}
\alias{compartmentalModel}
\title{
  Compartmental Models for Kinetic Parameter Estimation
}
\description{
  ...
}
\usage{
compartmentalModel(type)
}
\arguments{
  \item{type}{is a character string that identifies the type of
    compartmental model to be used.  Acceptable models include:
    \describe{
      \item{\dQuote{weinmann}}{Tofts & Kermode AIF convolved with a
	single compartment model}
      \item{\dQuote{extended}}{Kety model extended with additional
	vascular compartment (default)}
      \item{\dQuote{orton.exp}}{Extended model using Orton's exponential
	AIF} 
      \item{\dQuote{orton.cos}}{Extended model using Orton's raised
	cosine AIF}
      \item{\dQuote{kety.orton.exp}}{Kety model using Orton's
	exponential AIF}
      \item{\dQuote{kety.orton.cos}}{Kety model using Orton's raised
	cosine AIF}
      \item{weinmann.empirical}{User-specified empirical AIF convolved
	with a single compartment model}
      \item{extended.empirical}{Extended model using a user-specified
	empirical AIF}
    }
  }
}
\details{
...
}
\value{
...
}
% \references{}
\author{Brandon Whitcher \email{bjw34032@users.sourceforge.net}}
\seealso{
  \code{\link{aifParameters}}, \code{\link{dcemri.lm}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
