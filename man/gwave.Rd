\name{gwave}
\alias{gwave}
\title{
Gabor Functions on a Ridge
}
\description{
Generation of Gabor functions located on the ridge.
}
\usage{
gwave(bridge, omegaridge, nvoice, freqstep, scale, np, N)
}
\arguments{
\item{bridge}{
time coordinates of the ridge samples
}
\item{omegaridge}{
frequency coordinates of the ridge samples
}
\item{nvoice}{
number of different scales per octave
}
\item{freqstep}{
sampling rate for the frequency axis
}
\item{scale}{
scale of the window
}
\item{np}{
size of the reconstruction kernel
}
\item{N}{
number of complex constraints
}}
\value{
Array of Gabor functions located on the ridge samples
}
\details{
}
\references{
See discussions in the text of "Time-Frequency Analysis''.
}
\seealso{
\code{\link{gwave2}}, \code{\link{morwave}}, \code{\link{morwave2}}.
}
\keyword{ts}
