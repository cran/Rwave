\name{ext}
\alias{ext}
\title{
Extrema of Dyadic Wavelet Transform
}
\description{
Compute the local extrema of the dyadic wavelet transform modulus.
}
\usage{
ext(wt, scale=FALSE, plot=TRUE)
}
\arguments{
\item{wt}{
dyadic wavelet transform.
}
\item{scale}{
flag indicating if the extrema at each 
resolution will be plotted at the same scale.
}
\item{plot}{
if set to TRUE, displays the transform on the graphics device.
}}
\value{
Structure containing:
\item{original}{
  original signal.
}
\item{extrema}{
  extrema representation.
}
\item{Sf}{
  coarse resolution of signal.
}
\item{maxresoln}{
  number of decomposition scales.
}
\item{np}{
  size of signal.
}
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{mw}}, \code{\link{mrecons}}.
}
\keyword{ts}
