\name{skeleton2}
\alias{skeleton2}
\title{
Reconstruction from Dual Wavelet
}
\description{
Computes the reconstructed signal from the ridge in the
case of real constraints.
}
\usage{
skeleton2(cwtinput, Qinv, morvelets, bridge, aridge, N)
}
\arguments{
\item{cwtinput}{
continuous wavelet transform (as the output of cwt).
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
\item{lam}{coefficients of dual wavelets in reconstructed signal. They
  are the Lagrange multipliers \eqn{\lambda}{lambda}'s of the text.}
\item{dualwave}{array containing the dual wavelets.}
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{skeleton}}.
}
\keyword{ts}
