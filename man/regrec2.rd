\name{regrec2}
\alias{regrec2}
\title{
Reconstruction from a Ridge
}
\description{
Reconstructs signal from a ``regularly sampled'' ridge,
in the wavelet case, from a precomputed kernel.
}
\usage{
regrec2(siginput, cwtinput, phi, nbnodes, noct, nvoice, Q2,
epsilon=0.5, w0=2 * pi, plot=FALSE)
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
\item{nbnodes}{
number of samples on the ridge
}
\item{noct}{
number of octaves (powers of 2)
}
\item{nvoice}{
number of different scales per octave
}
\item{Q2}{
second term of the reconstruction kernel
}
\item{epsilon}{
coefficient of the \eqn{Q_2} term in reconstruction kernel
}
\item{w0}{
central frequency of Morlet wavelet
}
\item{plot}{
if set to TRUE, displays original and reconstructed signals
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
\item{nbnodes}{
  number of nodes used for the reconstruction.
}}
\details{
The computation of the kernel may be time consuming. This function
avoids recomputing it if it was computed already.
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{regrec}}, \code{\link{gregrec}}, \code{\link{ridrec}},
\code{\link{sridrec}}. 
}
\keyword{ts}
