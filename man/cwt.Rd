\name{cwt}
\alias{cwt}
\title{
Continuous Wavelet Transform
}
\description{
Computes the continuous wavelet transform with for the (complex-valued) 
Morlet wavelet.
}
\usage{
cwt(input, noctave, nvoice=1, w0=2 * pi, twoD=TRUE, plot=TRUE)
}
\arguments{
\item{input}{
input signal (possibly complex-valued)
}
\item{noctave}{
number of powers of 2 for the scale variable
}
\item{nvoice}{
number of scales in each octave
(i.e. between two consecutive powers of 2).
}
\item{w0}{
central frequency of the wavelet.
}
\item{twoD}{
logical variable set to T to organize the output as a 2D array 
(signal\_size x nb\_scales), otherwise, the output is a 3D array 
(signal\_size x noctave x nvoice).
}
\item{plot}{
if set to T, display the modulus of the
continuous wavelet transform on the graphic device.
}}
\value{
continuous (complex) wavelet transform
}
\details{
The output contains the (complex) values of the wavelet transform of the
input signal. The format of the output can be

2D array (signal\_size x nb\_scales)

3D array (signal\_size x noctave x nvoice)

Since Morlet's wavelet is not strictly speaking a wavelet (it is not of
vanishing integral), artifacts may occur for certain signals.
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwtp}}, \code{\link{cwtTh}}, \code{\link{DOG}},
\code{\link{gabor}}.
}
\keyword{ts}
