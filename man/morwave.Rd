\name{morwave}
\alias{morwave}
\title{
Ridge Morvelets
}
\description{
Generates the Morlet wavelets at the sample points of the ridge.
}
\usage{
morwave(bridge, aridge, nvoice, np, N, w0=2 * pi)
}
\arguments{
\item{bridge}{
time coordinates of the ridge samples.
}
\item{aridge}{
scale coordinates of the ridge samples.
}
\item{nvoice}{
number of different scales per octave.
}
\item{np}{
number of samples in the input signal.
}
\item{N}{
size of reconstructed signal.
}
\item{w0}{
central frequency of the wavelet.
}}
\value{
Returns the Morlet wavelets at the samples of the time-scale plane
given in the input:
complex array of Morlet wavelets located on the ridge samples
}
%\details{}
\references{
  See discussions in the text of \dQuote{Time-Frequency Analysis}.
}
\seealso{
\code{\link{morwave2}}, \code{\link{gwave}}, \code{\link{gwave2}}.
}
\keyword{ts}
