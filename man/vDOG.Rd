\name{vDOG}
\alias{vDOG}
\title{
DOG Wavelet Transform on one Voice
}
\description{
Compute DOG wavelet transform at one scale.
}
\usage{
vDOG(input, scale, moments)
}
\arguments{
\item{input}{
Input signal (1D array).
}
\item{scale}{
Scale at which the wavelet transform is to be computed.
}
\item{moments}{
number of vanishing moments.
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
\code{\link{vgt}}, \code{\link{vwt}}.
}
\keyword{ts}
