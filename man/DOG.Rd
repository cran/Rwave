\name{DOG}
\alias{DOG}
\title{
Continuous Wavelet Transform with
derivative of Gaussian
}
\description{
Computes the continuous wavelet transform with for (complex-valued) 
derivative of Gaussian wavelets.
}
\usage{
DOG(input, noctave, nvoice=1, moments, twoD=TRUE, plot=TRUE)
}
\arguments{
\item{input}{
input signal (possibly complex-valued).
}
\item{noctave}{
number of powers of 2 for the scale variable.
}
\item{moments}{
number of vanishing moments of the wavelet
(order of the derivative).
}
\item{nvoice}{
number of scales in each octave (i.e. between two
consecutive powers of 2)
}
\item{twoD}{
logical variable set to T to organize the
output as a 2D array 
(signal\_size x nb\_scales), otherwise, the output is a 3D array 
(signal\_size x noctave x nvoice)
}
\item{plot}{
if set to T, display the modulus of the
continuous wavelet transform on the graphic device
}}
\value{
continuous (complex) wavelet transform
}
\details{
The output contains the (complex) values of the wavelet transform of the
input signal.
The format of the output can be 


2D array (signal\_size x nb\_scales)


3D array (signal\_size x noctave x nvoice) 
}
\references{
  See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwt}}, \code{\link{cwtp}}, \code{\link{cwtsquiz}},
\code{\link{cgt}}.
}
\keyword{ts}
