\name{morwave2}
\alias{morwave2}
\title{
Real Ridge Morvelets
}
\description{
Generates the real parts of the Morlet wavelets at the sample 
points of a ridge
}
\usage{
morwave2(bridge, aridge, nvoice, np, N, w0=2 * pi)
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
Returns the real parts of the Morlet wavelets at the samples of 
the time-scale plane given in the input:
array of Morlet wavelets located on the ridge samples
}
\details{
}
\references{
See discussions in the text of ``Time-Frequency Analysis''.
}
\seealso{
\code{\link{morwave}}, \code{\link{gwave}}, \code{\link{gwave2}}.
}
\keyword{ts}
