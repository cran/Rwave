\name{cwtp}
\alias{cwtp}
\title{
Continuous Wavelet Transform with Phase Derivative
}
\description{
Computes the continuous wavelet transform with (complex-valued)
Morlet wavelet and its phase derivative.
}
\usage{
cwtp(input, noctave, nvoice=1, w0=2 * pi, twoD=TRUE, plot=TRUE)
}
\arguments{
\item{input}{
input signal (possibly complex-valued)
}
\item{noctave}{
number of powers of 2 for the scale variable
}
\item{nvoice}{
number of scales in each octave (i.e., between two consecutive powers of 2).
}
\item{w0}{
central frequency of the wavelet.
}
\item{twoD}{
logical variable set to \code{T} to organize the output as a 2D array
(signal size \eqn{\times}{x} nb scales), otherwise, the output is a 3D
array (signal size \eqn{\times}{x} noctave \eqn{\times}{x} nvoice).
}
\item{plot}{
if set to \code{T}, display the modulus of the continuous wavelet
transform on the graphic device.
}}
\value{
list containing the continuous (complex) wavelet transform
and the phase derivative
\item{wt}{
  array of complex numbers for the values of the continuous wavelet
  transform.
}
\item{f}{
  array of the same dimensions containing the values of the derivative
  of the phase of the continuous wavelet transform.}
}}
\details{
}
\references{
See discussions in the text of "Practical Time-Frequency Analysis".
}
\seealso{
\code{\link{cgt}}, \code{\link{cwt}}, \code{\link{cwtTh}},
\code{\link{DOG}} for wavelet transform, and \code{\link{gabor}} for
continuous Gabor transform.
}
\keyword{ts}
