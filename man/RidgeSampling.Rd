\name{RidgeSampling}
\alias{RidgeSampling}
\title{
Sampling Gabor Ridge
}
\description{
Given a ridge phi (for the Gabor transform), returns a 
(regularly) subsampled version of length nbnodes.
}
\usage{
RidgeSampling(phi, nbnodes)
}
\arguments{
\item{phi}{
ridge (1D array).
}
\item{nbnodes}{
number of samples.
}}
\value{
Returns a list containing the discrete values of the ridge.
\item{node}{
  time coordinates of the ridge samples.
}
\item{phinode}{
  frequency coordinates of the ridge samples.
}}
\details{
Gabor ridges are sampled uniformly.
}
\references{
See discussions in the text of "Time-Frequency Analysis''.
}
\seealso{
\code{\link{wRidgeSampling}}.
}
\keyword{ts}
