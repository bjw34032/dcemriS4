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
img <- round(255 * img)
img[17:32,,] <- 255 - img[17:32,,]
img.nifti <- nifti(img) # create NIfTI object

writeNIfTI(img.nifti, "test-image-uint8", verbose=TRUE)
## These files should be viewable in, for example, FSLview
## Make sure you adjust the min/max values for proper visualization
data <- readNIfTI("test-image-uint8", verbose=TRUE)
image(img.nifti) # image(img.nifti, z=16, plot.type="single")
image(data)      # image(data, z=16, zlim=range(data), plot.type="single")
abs.err <- abs(data - img.nifti)
image(as(abs.err, "nifti"), zlim=range(img.nifti))

## Loop through all possible data types
datatypes <- list(code=c(2, 4, 8, 16, 64),
                  name=c("uint8", "int16", "int32", "float", "double"))
equal <- vector("list")
for (i in 1:length(datatypes$code)) {
  fname <- paste("test-image-", datatypes$name[i], sep="")
  rm(img.nifti)
  img.nifti <- nifti(img, datatype=datatypes$code[i])
  writeNIfTI(img.nifti, fname)
  equal[[i]] <- all(readNIfTI(fname) == img)
}
names(equal) <- datatypes$name
unlist(equal)
}
\author{Brandon Whitcher \email{} and Volker Schmid \email{}}
\keyword{file}
