\name{vwt}
\alias{vwt}
\title{
 Voice Wavelet Transform
}
\description{
Compute Morlet's wavelet transform at one scale.
}
\usage{
vwt(input, scale, w0=2 * pi)
}
\arguments{
\item{input}{
Input signal (1D array).
}
\item{scale}{
Scale at which the wavelet transform is to be computed.
}
\item{w0}{
Center frequency of the wavelet.
}}
\value{
1D (complex) array containing wavelet transform at one scale.
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{vgt}}, \code{\link{vDOG}}.
}
\keyword{ts}
