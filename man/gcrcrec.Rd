\name{gcrcrec}
\alias{gcrcrec}
\title{
Crazy Climbers Reconstruction by Penalization
}
\description{
Reconstructs a real valued signal from ridges found by 
crazy climbers on a gabor transform.
}
\usage{
gcrcrec(siginput, inputgt, beemap, nvoice, freqstep, scale, compr,
bstep=5, ptile=0.01, epsilon=0, fast=TRUE, para=5, minnbnodes=3,
hflag=FALSE, real=FALSE, plot=2)
}
\arguments{
\item{siginput}{
original signal.
}
\item{inputgt}{
Gabor transform.
}
\item{beemap}{
occupation measure, output of \code{\link{crc}}.
}
\item{nvoice}{
number of frequencies.
}
\item{freqstep}{
sampling step for frequency axis.
}
\item{scale}{
size of windows.
}
\item{compr}{
compression rate to be applied to the ridges.
}
\item{bstep}{
size (in the time direction) of the steps for chaining.
}
\item{ptile}{
threshold of ridge
}
\item{epsilon}{
constant in front of the smoothness term in penalty function.
}
\item{fast}{
if set to TRUE, uses trapezoidal rule to evaluate \eqn{Q_2}.
}
\item{para}{
scale parameter for extrapolating the ridges.
}
\item{minnbnodes}{
minimal number of points per ridge.
}
\item{hflag}{
if set to FALSE, uses the identity as first term
in the kernel. If not, uses \eqn{Q_1} instead.
}
\item{real}{
if set to TRUE, uses only real constraints.
}
\item{plot}{
  \item{1}{
    displays signal,components, and reconstruction one after another.
  }
  \item{2}{displays signal, components and reconstruction.
}}
}
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
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{crc}}, \code{\link{cfamily}}, \code{\link{crcrec}},
\code{\link{scrcrec}}.
}
\keyword{ts}
