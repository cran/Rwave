\name{cwtimage}
\alias{cwtimage}
\title{
Continuous Wavelet Transform Display
}
\description{
Converts the output (modulus or argument) of cwtpolar to a 
2D array and displays on the graphic device.
}
\usage{
cwtimage(input)
}
\arguments{
\item{input}{
3D array containing a continuous wavelet transform
}}
\value{
2D array continuous (complex) wavelet transform
}
\details{
The output contains the (complex) values of the wavelet transform of the
input signal. The format of the output can be 


2D array (signal\_size x nb\_scales)


3D array (signal\_size x noctave x nvoice) 
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwtpolar}}, \code{\link{cwt}}, \code{\link{DOG}}.
}
\keyword{ts}
