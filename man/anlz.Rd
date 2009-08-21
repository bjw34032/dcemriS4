% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id: $

\name{anlz}
\Rdversion{1.1}
\alias{anlz}
\alias{is.anlz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Constructor for Analyze}
\description{
  ...
}
\usage{
anlz(img = array(0, dim = rep(1, 4)), dim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{img}{
%%     ~~Describe \code{img} here~~
}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
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
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(img=array(0, dim=rep(1,4)), dim, ...) {
  if (missing(dim)) {
    if (is.array(img))
      dim <- base::dim(img)
    else
      dim <- c(1, length(img))
  }
  ld <- length(dim)
  if (ld < 3)
    stop(sprintf("length(dim) must be at least 3 and is \%d.", ld))
  
  x <- c(length(dim), dim[1], dim[2], dim[3],
         ifelse(is.na(dim[4]), 1, dim[4]), rep(1,3))
  y <- c(0.0, rep(1.0,length(dim)), rep(0.0,3))
  cal.max <- quantile(img, probs=0.95, na.rm=TRUE)
  cal.min <- quantile(img, probs=0.05, na.rm=TRUE)
  obj <- new("anlz", .Data=array(img, dim=dim), "dim_"=x, "pixdim"=y,
             "cal_max"=cal.max, "cal_min"=cal.min, ...)
  validObject(obj)
  return(obj)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
