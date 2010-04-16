\name{smoothwt}
\alias{smoothwt}
\title{
Smoothing and Time Frequency Representation
}
\description{
smooth the wavelet (or Gabor) transform in the time direction.
}
\usage{
smoothwt(modulus, subrate, flag=FALSE)
}
\arguments{
\item{modulus}{
Time-Frequency representation (real valued).
}
\item{subrate}{
Length of smoothing window.
}
\item{flag}{
If set to TRUE, subsample the representation.
}}
\value{
2D array containing the smoothed transform.
}
%\details{}
\references{
See discussions in the text of \dQuote{Time-Frequency Analysis}.
}
\seealso{
\code{\link{corona}}, \code{\link{coronoid}}, \code{\link{snake}},
\code{\link{snakoid}}.
}
\keyword{ts}
