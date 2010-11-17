\name{aifParameters}
\alias{aifParameters}
\title{
  Parameters for Arterial Input Functions
}
\description{
  Specification of parameters for arterial input functions (AIFs) 
}
\usage{
aifParameters(type, user = NULL)
}
\arguments{
  \item{type}{is one of the following character strings associated with
    an AIF:
    \itemize{
      \item{\code{tofts.kermode}}
      \item{\code{fritz.hansen}}
      \item{\code{orton.exp}}
      \item{\code{orton.cos}}
      \item{\code{user}}
      \item{\code{empirical}}
    }
  }
  \item{user}{is a vector of estimated AIF parameters or the empirical
    AIF values.}
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
  \code{\link{compartmentalModel}}, \code{\link{dcemri.lm}}
}
% \examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
