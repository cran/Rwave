\name{vgt}
\alias{vgt}
\title{
Gabor Transform on one Voice
}
\description{
Compute Gabor transform for fixed frequency.
}
\usage{
vgt(input, frequency, scale, plot=FALSE)
}
\arguments{
\item{input}{
Input signal (1D array).
}
\item{frequency}{
frequency at which the Gabor transform is to be computed.
}
\item{scale}{
frequency at which the Gabor transform is to be computed.
}
\item{plot}{
if set to TRUE, plots the real part of cgt on the graphic device.
}}
\value{
1D (complex) array containing Gabor transform at specified frequency.
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{vwt}}, \code{\link{vDOG}}.
}
\keyword{ts}
