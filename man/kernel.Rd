\name{kernel}
\alias{kernel}
\title{
Kernel for Reconstruction from Wavelet Ridges
}
\description{
Computes the cost from the sample of points on the estimated ridge
and the matrix used in the reconstruction of the original signal
}
\usage{
kernel(node, phinode, nvoice, x.inc=1, x.min=node[1],
x.max=node[length(node)], w0=2 * pi, plot=FALSE)
}
\arguments{
\item{node}{
values of the variable b for the nodes of the ridge.
}
\item{phinode}{
values of the scale variable a for the nodes of the ridge.
}
\item{nvoice}{
number of scales within 1 octave.
}
\item{x.inc}{
step unit for the computation of the kernel.
}
\item{x.min}{
minimal value of x for the computation of \eqn{Q_2}.
}
\item{x.max}{
maximal value of x for the computation of \eqn{Q_2}.
}
\item{w0}{
central frequency of the wavelet.
}
\item{plot}{
if set to TRUE, displays the modulus of the matrix of \eqn{Q_2}.
}}
\value{
matrix of the \eqn{Q_2} kernel
}
\details{
The kernel is evaluated using Romberg's method.
}
\references{
See discussions in the text of "Time-Frequency Analysis".
}
\seealso{
\code{\link{gkernel}}, \code{\link{rkernel}}, \code{\link{zerokernel}}.
}
\keyword{ts}
