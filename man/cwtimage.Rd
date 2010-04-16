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
% format of the input??

2D array (signal\_size x nb\_scales)

3D array (signal\_size x noctave x nvoice) 
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwtpolar}}, \code{\link{cwt}}, \code{\link{DOG}}.
}
\examples{
    x <- 1:512
    chirp <- sin(2*pi * (x + 0.002 * (x-256)^2 ) / 16)
    retChirp <- cwt(chirp, noctave=5, nvoice=12, twoD=FALSE, plot=FALSE)
    retPolar <- cwtpolar(retChirp)
    retImageMod <- cwtimage(retPolar$modulus)
    retImageArg <- cwtimage(retPolar$argument)
}
\keyword{ts}
