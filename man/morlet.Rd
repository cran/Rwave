\name{morlet}
\alias{morlet}
\title{
 Morlet Wavelets
}
\description{
Computes a Morlet wavelet at the point of the time-scale
plane given in the input
}
\usage{
morlet(sigsize, location, scale, w0=2 * pi)
}
\arguments{
\item{sigsize}{
length of the output.
}
\item{location}{
time location of the wavelet.
}
\item{scale}{
scale of the wavelet.
}
\item{w0}{
central frequency of the wavelet.
}}
\value{
Returns the values of the complex Morlet wavelet at the point of
the time-scale plane given in the input
}
\details{
The details of this construction (including the definition formulas)
are given in the text.
}
\references{


See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{gabor}}.
}
\keyword{ts}
