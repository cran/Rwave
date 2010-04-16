\name{regrec}
\alias{regrec}
\title{
Reconstruction from a Ridge
}
\description{
Reconstructs signal from a \dQuote{regularly sampled} ridge, in the wavelet 
case.
}
\usage{
regrec(siginput, cwtinput, phi, compr, noct, nvoice, epsilon=0,
w0=2 * pi, fast=FALSE, plot=FALSE, para=0, hflag=FALSE,
check=FALSE, minnbnodes=2, real=FALSE)
}
\arguments{
\item{siginput}{
input signal.
}
\item{cwtinput}{
wavelet transform, output of \code{\link{cwt}}.
}
\item{phi}{
unsampled ridge.
}
\item{compr}{
subsampling rate for the wavelet coefficients (at scale 1)
}
\item{noct}{
number of octaves (powers of 2)
}
\item{nvoice}{
number of different scales per octave
}
\item{epsilon}{
coefficient of the \eqn{Q_2} term in reconstruction kernel
}
\item{w0}{
central frequency of Morlet wavelet
}
\item{fast}{
if set to TRUE, the kernel is computed using trapezoidal rule.
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
\item{check}{
if set to TRUE, computes \code{\link{cwt}} of reconstructed signal.
}
\item{minnbnodes}{
minimum number of nodes for the reconstruction.
}
\item{real}{
if set to TRUE, uses only real constraints on the transform.
}}
\value{
Returns a list containing:
\item{sol}{
  reconstruction from a ridge.
}
\item{A}{
  <wavelets,dualwavelets> matrix.
}
\item{lam}{
  coefficients of dual wavelets in reconstructed signal.
}
\item{dualwave}{
  array containing the dual wavelets.
}
\item{morvelets}{
  array containing the wavelets on sampled ridge.
}
\item{solskel}{
  wavelet transform of sol, restricted to the ridge.
}
\item{inputskel}{
  wavelet transform of signal, restricted to the ridge.
}
\item{Q2}{
  second part of the reconstruction kernel.
}
\item{nbnodes}{
  number of nodes used for the reconstruction.
}}
%\details{}
\references{
See discussions in the text of \dQuote{Practical Time-Frequency Analysis}.
}
\seealso{
\code{\link{regrec2}}, \code{\link{ridrec}}, \code{\link{gregrec}},
\code{\link{gridrec}}.
}
\keyword{ts}
