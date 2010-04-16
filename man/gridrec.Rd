\name{gridrec}
\alias{gridrec}
\title{
Reconstruction from a Ridge
}
\description{
Reconstructs signal from sample of a ridge, in the Gabor case.
}
\usage{
gridrec(gtinput, node, phinode, nvoice, freqstep, scale, Qinv,
epsilon, np, real=FALSE, check=FALSE)
}
\arguments{
\item{gtinput}{
Gabor transform, output of \code{\link{cgt}}.
}
\item{node}{
time coordinates of the ridge samples.
}
\item{phinode}{
frequency coordinates of the ridge samples.
}
\item{nvoice}{
number of different frequencies.
}
\item{freqstep}{
sampling rate for the frequency axis.
}
\item{scale}{
scale of the window.
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
\item{real}{
if set to TRUE, uses only constraints on the real part
of the transform.
}
\item{check}{
if set to TRUE, computes \code{\link{cgt}} of
reconstructed signal.
}}
\value{
Returns a list containing the reconstructed signal and the chained ridges.
\item{sol}{
  reconstruction from a ridge.
}
\item{A}{
  <gaborlets,dualgaborlets> matrix.
}
\item{lam}{
  coefficients of dual gaborlets in reconstructed signal.
}
\item{dualwave}{
  array containing the dual gaborlets.
}
\item{gaborets}{
  array of gaborlets located on the ridge samples.
}
\item{solskel}{
  Gabor transform of sol, restricted to the ridge.
}
\item{inputskel}{
  Gabor transform of signal, restricted to the ridge.
}}
%\details{}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{sridrec}}, \code{\link{gregrec}}, \code{\link{regrec}},
\code{\link{regrec2}}.
}
\keyword{ts}
