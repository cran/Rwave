\name{ridrec}
\alias{ridrec}
\title{
Reconstruction from a Ridge
}
\description{
Reconstructs signal from sample of a ridge,
in the wavelet case.
}
\usage{
ridrec(cwtinput, node, phinode, noct, nvoice, Qinv, epsilon, np,
w0=2 * pi, check=FALSE, real=FALSE)
}
\arguments{
\item{cwtinput}{
wavelet transform, output of \code{\link{cwt}}.
}
\item{node}{
time coordinates of the ridge samples.
}
\item{phinode}{
scale coordinates of the ridge samples.
}
\item{noct}{
number of octaves (powers of 2).
}
\item{nvoice}{
number of different scales per octave.
}
\item{Qinv}{
inverse of the matrix \eqn{Q} of the quadratic form.
}
\item{epsilon}{
coefficient of the \eqn{Q_2} term in reconstruction kernel
}
\item{np}{
number of samples of the reconstructed signal.
}
\item{w0}{
central frequency of Morlet wavelet.
}
\item{check}{
if set to TRUE, computes \code{\link{cwt}} of reconstructed signal.
}
\item{real}{
if set to TRUE, uses only constraints on the real part of the transform.
}}
\value{
Returns a list containing the reconstructed signal and the chained ridges.
\item{sol}{reconstruction from a ridge}
\item{A}{<wavelets,dualwavelets> matrix }
\item{lam}{coefficients of dual wavelets in reconstructed signal.}
\item{dualwave}{array containing the dual wavelets.}
\item{morvelets}{array of morlet wavelets located on the ridge samples.}
\item{solskel}{wavelet transform of sol, restricted to the ridge}
\item{inputskel}{wavelet transform of signal, restricted to the ridge}
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{sridrec}}, \code{\link{regrec}}, \code{\link{regrec2}}.
}
\keyword{ts}
