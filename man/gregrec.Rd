\name{gregrec}
\alias{gregrec}
\title{
Reconstruction from a Ridge
}
\description{
Reconstructs signal from a ``regularly sampled'' ridge, in the Gabor case.
}
\usage{
gregrec(siginput, gtinput, phi, nbnodes, nvoice, freqstep, scale,
epsilon=0, fast=FALSE, plot=FALSE, para=0, hflag=FALSE, real=FALSE,
check=FALSE)
}
\arguments{
\item{siginput}{
input signal.
}
\item{gtinput}{
Gabor transform, output of \code{\link{cgt}}.
}
\item{phi}{
unsampled ridge.
}
\item{nbnodes}{
number of nodes used for the reconstruction.
}
\item{nvoice}{
number of different scales per octave
}
\item{freqstep}{
sampling rate for the frequency axis
}
\item{scale}{
size parameter for the Gabor function.
}
\item{epsilon}{
coefficient of the \eqn{Q_2} term in reconstruction kernel
}
\item{fast}{
if set to T, the kernel is computed using trapezoidal rule.
}
\item{plot}{
if set to TRUE, displays original and reconstructed signals
}
\item{para}{
scale parameter for extrapolating the ridges.
}
\item{hflag}{
if set to TRUE, uses \eqn{Q_1} as first term in the kernel.
}
\item{real}{
if set to TRUE, uses only real constraints on the transform.
}
\item{check}{
if set to TRUE, computes \code{\link{cwt}} of
reconstructed signal.
}}
\value{
Returns a list containing:
\item{sol}{
  reconstruction from a ridge.
}
\item{A}{
  <gaborlets,dualgaborlets> matrix.
}
\item{lam}{
  coefficients of dual wavelets in reconstructed signal.
}
\item{dualwave}{
  array containing the dual wavelets.
}
\item{gaborets}{
  array containing the wavelets on sampled ridge.
}
\item{solskel}{
  Gabor transform of sol, restricted to the ridge.
}
\item{inputskel}{
  Gabor transform of signal, restricted to the ridge.
}
\item{Q2}{
  second part of the reconstruction kernel.
}}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{regrec}}.
}
\keyword{ts}
