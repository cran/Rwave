\name{gkernel}
\alias{gkernel}
\title{
Kernel for Reconstruction from Gabor Ridges
}
\description{
Computes the cost from the sample of points on the estimated ridge
and the matrix used in the reconstruction of the original signal.
}
\usage{
gkernel(node, phinode, freqstep, scale, x.inc=1, x.min=node[1],
x.max=node[length(node)], plot=FALSE)
}
\arguments{
\item{node}{
values of the variable b for the nodes of the ridge.
}
\item{phinode}{
values of the scale variable a for the nodes of the ridge.
}
\item{freqstep}{
sampling rate for the frequency axis.
}
\item{scale}{
size of the window.
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
\item{plot}{
if set to TRUE, displays the modulus of the matrix of \eqn{Q_2}.
}}
\value{
matrix of the \eqn{Q_2} kernel
}
%\details{}
\references{
See discussions in the text of ``Time-Frequency Analysis''.
}
\seealso{
\code{\link{fastgkernel}}, \code{\link{kernel}}, \code{\link{rkernel}},
\code{\link{fastkernel}}, \code{\link{zerokernel}}.
}
\keyword{ts}
