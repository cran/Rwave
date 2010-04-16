\name{gwave2}
\alias{gwave2}
\title{
Real Gabor Functions on a Ridge
}
\description{
Generation of the real parts of gabor functions located on a ridge.
(modification of \code{\link{gwave}}.)
}
\usage{
gwave2(bridge, omegaridge, nvoice, freqstep, scale, np, N)
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
Array of real Gabor functions located on the ridge samples
}
%\details{}
\references{
See discussions in the text of \dQuote{Time-Frequency Analysis}.
}
\seealso{
\code{\link{gwave}}, \code{\link{morwave}}, \code{\link{morwave2}}.
}
\keyword{ts}
