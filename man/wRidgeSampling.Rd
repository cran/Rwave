\name{wRidgeSampling}
\alias{wRidgeSampling}
\title{
Sampling wavelet Ridge
}
\description{
Given a ridge \eqn{\phi}{phi} (for the wavelet transform), returns a
(appropriately) subsampled version with a given subsampling rate.
}
\usage{
wRidgeSampling(phi, compr, nvoice)
}
\arguments{
\item{phi}{
ridge (1D array).
}
\item{compr}{
subsampling rate for the ridge.
}
\item{nvoice}{
number of voices per octave.
}}
\value{
Returns a list containing the discrete values of the ridge.
\item{node}{time coordinates of the ridge samples.}
\item{phinode}{scale coordinates of the ridge samples.}
\item{nbnode}{number of nodes of the ridge samples.}
}
\details{
To account for the variable sizes of wavelets, the sampling
rate of a wavelet ridge is not uniform, and is proportional
to the scale.
}
%\references{}
\seealso{
\code{\link{RidgeSampling}}.
}
\keyword{ts}
