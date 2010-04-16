\name{vecmorlet}
\alias{vecmorlet}
\title{
Morlet Wavelets on a Ridge
}
\description{
Generate Morlet wavelets at specified positions on a ridge.
}
\usage{
vecmorlet(sigsize, nbnodes, bridge, aridge, w0=2 * pi)
}
\arguments{
\item{sigsize}{
Signal size.
}
\item{nbnodes}{
Number of wavelets to be generated.
}
\item{bridge}{
b coordinates of the ridge samples (1D array of length nbnodes).
}
\item{aridge}{
a coordinates of the ridge samples (1D array of length nbnodes).
}
\item{w0}{
Center frequency of the wavelet.
}}
\value{
2D (complex) array containing wavelets located at the specific points.
}
%\details{}
%\references{}
\seealso{
\code{\link{vecgabor}}.
}
\keyword{ts}
