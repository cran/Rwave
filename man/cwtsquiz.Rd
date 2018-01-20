\name{cwtsquiz}
\alias{cwtsquiz}
\title{
Squeezed Continuous Wavelet Transform
}
\description{
Computes the synchrosqueezed continuous wavelet transform with
the (complex-valued) Morlet wavelet.
}
\usage{
cwtsquiz(input, noctave, nvoice=1, w0=2 * pi, twoD=TRUE, plot=TRUE)
}
\arguments{
\item{input}{
input signal (possibly complex-valued)
}
\item{noctave}{
number of powers of 2 for the scale variable
}
\item{nvoice}{
number of scales in each octave (i.e. between two consecutive powers of 2).
}
\item{w0}{
central frequency of the wavelet.
}
\item{twoD}{
logical variable set to T to organize the output as a 2D array
(signal size \eqn{\times}{x} nb scales), otherwise, the output is a 3D
array (signal size \eqn{\times}{x} noctave \eqn{\times}{x} nvoice).
}
\item{plot}{
logical variable set to T to T to display the modulus of the squeezed
wavelet transform on the graphic device.
}}
\value{
synchrosqueezed continuous (complex) wavelet transform
}
\details{
The output contains the (complex) values of the squeezed wavelet
transform of the input signal.  The format of the output can be 

2D array (signal size \eqn{\times}{x} nb scales),

3D array (signal size \eqn{\times}{x} noctave \eqn{\times}{x} nvoice).
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwt}}, \code{\link{cwtp}}, \code{\link{DOG}},
\code{\link{cgt}}.
}
\keyword{ts}
