% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{writeNIfTI}
\alias{writeNIfTI}
\title{writeNIfTI}
\description{
  ...
}
\usage{
writeNIfTI(nim, filename, gzipped = TRUE, verbose = FALSE, warn = -1)
}
\arguments{
  \item{nim}{}
  \item{filename}{is the path and file name to save the NIfTI file
    (.nii) \bold{without} the suffix.}
  \item{gzipped}{is a character string that enables exportation of
    compressed (.gz) files (default = \code{TRUE}).}
  \item{verbose}{}
  \item{warn}{is a number to regulate the display of warnings (default =
    -1).  See \code{\link{options}} for more details.}
}
\details{
  ...
}
\value{
  ...
}
\references{
  NIfTI-1 \url{http://nifti.nimh.nih.gov/}
}
\seealso{
  \code{\link{writeANALYZE}}
}
\examples{
norm <- dnorm(seq(-5, 5, length=32), sd=2)
norm <- (norm-min(norm)) / max(norm-min(norm))
img <- outer(outer(norm, norm), norm)
img <- round(255*img)
img[17:32,,] <- 255 - img[17:32,,]
X <- nrow(img)
Y <- ncol(img)
Z <- nsli(img)

## NIfTI
nim <- nifti(img)
\dontrun{
writeNIfTI(nim, "test-image-int16")
## These files should be viewable in, for example, FSLview
## Make sure you adjust the min/max values for proper visualization
data <- read.img("test-image-int16.nii.gz")
par(mfrow=c(6,6), mar=rep(0,4))
for (z in 1:32)
  image(img[,,z], zlim=range(img), col=grey(0:255/255), xlab="",
        ylab="", axes=FALSE)
for (z in 1:32)
  image(data[,,z,1], zlim=range(img), col=grey(0:255/255), xlab="",
        ylab="", axes=FALSE)
for (z in 1:32)
  image(abs(data[,,z,1] - img[,,z]), zlim=range(img),
        col=grey(0:255/255), xlab="", ylab="", axes=FALSE)
}

## Loop through all possible data types
datatypes <- c("uint8", "int16", "int32", "float", "double")
equal <- vector("list")
for (x in 1:length(datatypes)) {
  fname <- paste("test-image-", datatypes[x], sep="")
  write.nifti.img(fname, nhdr, img, datatypes[x])
  equal[[x]] <- all(drop(read.img(fname)) == img)
}
names(equal) <- datatypes
unlist(equal)

## Analyze (datatypes are the same as NIfTI)
hdr <- make.hdr(X, Y, Z, 1, "int16")
equal <- vector("list")
for (x in 1:length(datatypes)) {
  fname <- paste("test-analyze-image-", datatypes[x], sep="")
  write.analyze.img(fname, hdr, img, datatypes[x])
  equal[[x]] <- all(drop(read.img(fname)) == img)
}
names(equal) <- datatypes
unlist(equal)
}
\author{Brandon Whitcher and Volker Schmid}
%\keyword{}
