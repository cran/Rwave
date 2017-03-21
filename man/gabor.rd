\name{gabor}
\alias{gabor}
\title{
Generate Gabor function
}
\description{
Generates a Gabor for given location and frequency.
}
\usage{
gabor(sigsize, location, frequency, scale)
}
\arguments{
\item{sigsize}{
length of the Gabor function.
}
\item{location}{
position of the Gabor function.
}
\item{frequency}{
frequency of the Gabor function.
}
\item{scale}{
size parameter for the Gabor function. See details.
}}
\value{
complex 1D array of size sigsize.
}
\details{The size parameter here corresponds to the
standard deviation for a gaussian.  In the Carmona 1998 book, equation 3.23 has a different scale factor. }
\references{
See discussions in the text of \dQuote{Practical Time-Frequency Analysis}.
}
\seealso{
\code{\link{morlet}}.
}
\keyword{ts}
