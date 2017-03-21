\name{cgt}
\alias{cgt}
\title{
Continuous Gabor Transform
}
\description{
Computes the continuous Gabor transform with Gaussian window.
}
\usage{
cgt(input, nvoice, freqstep=(1/nvoice), scale=1, plot=TRUE)
}
\arguments{
\item{input}{
input signal (possibly complex-valued).
}
\item{nvoice}{
number of frequencies for which gabor transform
is to be computed.
}
\item{freqstep}{
Sampling rate for the frequency axis.
}
\item{scale}{
Size parameter for the window.
}
\item{plot}{
logical variable set to TRUE to display the modulus of the
continuous gabor transform on the graphic device.
}}
\value{
continuous (complex) gabor transform (2D array).
}
\details{
The output contains the (complex) values of the gabor transform of the
input signal. The format of the output is
a 2D array (signal\_size x nb\_scales).
}
\references{
See discussion in text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwt}}, \code{\link{cwtp}}, \code{\link{DOG}} for continuous wavelet transforms.
\code{\link{cwtsquiz}} for synchrosqueezed wavelet transform.
}
\keyword{ts}
\section{Warning}{
freqstep must be less than 1/nvoice to avoid aliasing. freqstep=1/nvoice
corresponds to the Nyquist limit.
}
