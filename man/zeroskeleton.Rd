\name{zeroskeleton}
\alias{zeroskeleton}
\title{
Reconstruction from Dual Wavelets
}
\description{
Computes the the reconstructed signal from the ridge when the epsilon
parameter is set to zero
}
\usage{
zeroskeleton(cwtinput, Qinv, morvelets, bridge, aridge, N)
}
\arguments{
\item{cwtinput}{
continuous wavelet transform (as the output of \code{\link{cwt}}).
}
\item{Qinv}{
inverse of the reconstruction kernel (2D array).
}
\item{morvelets}{
array of Morlet wavelets located at the ridge samples.
}
\item{bridge}{
time coordinates of the ridge samples.
}
\item{aridge}{
scale coordinates of the ridge samples.
}
\item{N}{
size of reconstructed signal.
}}
\value{
Returns a list of the elements of the reconstruction of a signal from
sample points of a ridge
\item{sol}{reconstruction from a ridge.}
\item{A}{matrix of the inner products.}
\item{lam}{coefficients of dual wavelets in reconstructed signal. They are the
Lagrange multipliers \eqn{\lambda}{lambda}'s of the text.}
\item{dualwave}{array containing the dual wavelets.}
}
\details{
The details of this reconstruction are the same as for the function
skeleton. They can be found in the text
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
  \code{\link{skeleton}}, \code{\link{skeleton2}},
  \code{\link{zeroskeleton2}}.
}
\keyword{ts}
