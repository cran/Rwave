\name{crcrec}
\alias{crcrec}
\title{
Crazy Climbers Reconstruction by Penalization
}
\description{
Reconstructs a real valued signal from the output of \code{\link{crc}}
(wavelet case) by minimizing an appropriate quadratic form.
}
\usage{
crcrec(siginput, inputwt, beemap, noct, nvoice, compr, minnbnodes=2,
w0=2 * pi, bstep=5, ptile=0.01, epsilon=0, fast=FALSE, para=5, real=FALSE,
plot=2)
}
\arguments{
\item{siginput}{
original signal.
}
\item{inputwt}{
wavelet transform.
}
\item{beemap}{
occupation measure, output of \code{\link{crc}}.
}
\item{noct}{
number of octaves.
}
\item{nvoice}{
number of voices per octave.
}
\item{compr}{
compression rate for sampling the ridges.
}
\item{minnbnodes}{
minimal number of points per ridge.
}
\item{w0}{
center frequency of the wavelet.
}
\item{bstep}{
size (in the time direction) of the steps for chaining.
}
\item{ptile}{
relative threshold of occupation measure.
}
\item{epsilon}{
constant in front of the smoothness term in penalty function.
}
\item{fast}{
if set to TRUE, uses trapezoidal rule to evaluate $Q_2$.
}
\item{para}{
scale parameter for extrapolating the ridges.
}
\item{real}{
if set to TRUE, uses only real constraints.
}
\item{plot}{
1: displays signal,components,and
reconstruction one after another. 2: displays
signal, components and reconstruction.
}}
\value{
Returns a structure containing the following elements:
\item{rec}{
  reconstructed signal.
}
\item{ordered}{
  image of the ridges (with different colors).
}
\item{comp}{
  2D array containing the signals reconstructed from ridges.
}
}
\details{
When ptile is high, boundary effects may appeare.
para controls extrapolation of the ridge.
}
\references{
}
\seealso{
\code{\link{crc}}, \code{\link{cfamily}}, \code{\link{scrcrec}}.
}
\keyword{ts}
