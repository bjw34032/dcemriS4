% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{nifti}
\Rdversion{1.1}
\alias{nifti}
\alias{is.nifti}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructor for NIfTI}
\description{
  ...
}
\usage{
nifti(img = array(0, dim = rep(1, 4)), dim, datatype = 2, ...)
is.nifti(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{}
  \item{dim}{}
  \item{datatype}{}
  \item{\dots}{}
  \item{x}{}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\seealso{
  \code{\linkS4class{nifti}}, \code{\link{anlz}}
}
\examples{
nim <- nifti() # default
nim
nim <- nifti(datatype=4) # 2-byte integers
nim
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
