% Part of the dcemri package for R
% Distributed under the BSD license: see dcemri/COPYING
%
% $Id$

\name{readNIfTI}
\alias{readNIfTI}
\alias{read.nifti.content}
\title{readNIfTI}
\description{
  These functions read in the header information and multidimensional
  array from a binary file in NIfTI-1 format into a
  \code{\linkS4class{nifti}}-class object.
}
\usage{
readNIfTI(fname, verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)
read.nifti.content(fname, onefile=TRUE, gzipped=TRUE, verbose=FALSE,
                   warn=-1, reorient=FALSE, call=NULL)
}
\arguments{
  \item{fname}{is the file name of the NIfTI file(s).}
  \item{verbose}{is a logical variable (default = \code{FALSE}) that
    allows text-based feedback during execution of the function.}
  \item{warn}{is a number to regulate the display of warnings (default =
    -1).  See \code{options} for more details.}
  \item{reorient}{}
  \item{onefile}{(logical) allows the scanning of single-file format
    (.nii) or dual-file format (.hdr and .img) NIfTI files (default =
    \code{TRUE}).} 
  \item{gzipped}{(logical) allows the importation of compressed (.gz)
    files (default = \code{TRUE}).}
  \item{call}{}
%  \item{ignoreQform}{is a logical variable (default = \code{FALSE}) that
%    ignores any \sQuote{qform} infomration in the NIfTI header.}
%  \item{ignoreSform}{is a logical variable (default = \code{FALSE}) that
%    ignores any \sQuote{sform} infomration in the NIfTI header.}
}
\details{
  The \code{readNIfTI} function utilizes internal methods \code{readBin}
  and \code{readChar} to efficiently extract information from the binary
  file(s).

  Current acceptable data types include
  \describe{
    \item{\dQuote{UINT8}}{DT BINARY (1 bit per voxel)}
    \item{\dQuote{INT16}}{DT SIGNED SHORT (16 bits per voxel)}
    \item{\dQuote{INT32}}{DT SINGED INT (32 bits per voxel)}
    \item{\dQuote{FLOAT32}}{DT FLOAT (32 bits per voxel)}
    \item{\dQuote{DOUBLE64}}{DT DOUBLE (64 bits per voxel)}
    \item{\dQuote{UINT16}}{DT UNSIGNED SHORT (16 bits per voxel)}
    }
}
\value{
  ...
}
\references{
  NIfTI-1\cr
  \url{http://nifti.nimh.nih.gov/}
}
\seealso{
  \code{\link{readANALYZE}}
}
\examples{
nim <- readNIfTI(system.file("nifti/avg152T1_LR_nifti.nii.gz",
                 package="dcemriS4")[1])
image(nim)
}
\author{
  Brandon Whitcher <\email{bjw34032@users.sourceforge.net}>,
  Volker Schmid <\email{volkerschmid@users.sourceforge.net}>,
  Andrew Thornton <\email{zeripath@users.sourceforge.net}>
}
%\keyword{}
