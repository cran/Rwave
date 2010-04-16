\name{cwtpolar}
\alias{cwtpolar}
\title{
Conversion to Polar Coordinates
}
\description{
Converts one of the possible outputs of the function \code{\link{cwt}}
to modulus and phase.
}
\usage{
cwtpolar(cwt, threshold=0)
}
\arguments{
\item{cwt}{
3D array containing the values of a continuous wavelet
transform in the format (signal size \eqn{\times}{x} noctave
\eqn{\times}{x} nvoice) as in the output of the function
\code{\link{cwt}} with the logical flag \code{twodimension} set to
FALSE.
}
\item{threshold}{
value of a level for the absolute value of the modulus below which 
the value of the argument of the output is set to \eqn{-\pi}{-pi}.
}}
\value{
Modulus and Argument of the values of the continuous wavelet transform
\item{output1}{
  3D array giving the values (in the same format as the input) of the
  modulus of the input.
}
\item{output2}{
  3D array giving the values of the argument of the input.
}}
\details{
The output contains the (complex) values of the wavelet transform of the
input signal. The format of the output can be 

2D array (signal size \eqn{\times}{x} nb\_scales)

3D array (signal size \eqn{\times}{x} noctave \eqn{\times}{x} nvoice) 
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{cwt}}, \code{\link{DOG}}, \code{\link{cwtimage}}.
}
\examples{
    x <- 1:512
    chirp <- sin(2*pi * (x + 0.002 * (x-256)^2 ) / 16)
    retChirp <- cwt(chirp, noctave=5, nvoice=12, twoD=FALSE, plot=FALSE)
    retPolar <- cwtpolar(retChirp)
}
\keyword{ts}
