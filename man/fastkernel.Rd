\name{fastkernel}
\alias{fastkernel}
\title{
Kernel for Reconstruction from Wavelet Ridges
}
\description{
Computes the cost from the sample of points on the estimated ridge
and the matrix used in the reconstruction of the original signal,
using simple trapezoidal rule for integrals.
}
\usage{
fastkernel(node, phinode, nvoice, x.inc=1, x.min=node[1],
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
step unit for the computation of the kernel
}
\item{x.min}{
minimal value of x for the computation of \eqn{Q_2}.
}
\item{x.max}{
maximal value of x for the computation of \eqn{Q_2}.
}
\item{w0}{
central frequency of the wavelet
}
\item{plot}{
if set to TRUE, displays the modulus of the matrix of \eqn{Q_2}.
}}
\value{
matrix of the \eqn{Q_2} kernel.
}
\details{
Uses trapezoidal rule (instead of Romberg's method) to evaluate
the kernel.
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.  
}
\seealso{
\code{\link{kernel}}, \code{\link{rkernel}}, \code{\link{gkernel}},
\code{\link{zerokernel}}.
}
\keyword{ts}
\examples{
# The function is currently defined as
function(node, phinode, nvoice, x.inc = 1, x.min = node[1], x.max = node[length(node)], w0 = 2 * pi, plot = F)
{
#########################################################################
#     fastkernel:   
#     -----------
#	   Same as kernel, except that the kernel is computed
#	     using Riemann sums instead of Romberg integration.
#
#     Input:
#     ------
#      node: values of the variable b for the nodes of the ridge
#      phinode: values of the scale variable a for the nodes of the ridge
#      nvoice: number of scales within 1 octave
#      x.inc: step unit for the computation of the kernel
#      x.min: minimal value of x for the computation of Q2
#      x.max: maximal value of x for the computation of Q2
#      w0: central frequency of the wavelet
#      plot: if set to TRUE, displays the modulus of the matrix of Q2
#
#     Output:
#     -------
#      ker: matrix of the Q2 kernel
#
#########################################################################
	lng <- as.integer((x.max - x.min)/x.inc) + 1
	nbnode <- length(node)
	b.start <- x.min - 50
	b.end <- x.max + 50
	ker.r <- matrix(0, lng, lng)
	ker.i <- matrix(0, lng, lng)
	dim(ker.r) <- c(lng * lng, 1)
	dim(ker.i) <- c(lng * lng, 1)
	phinode <- 2 * 2^(phinode/nvoice)
	z <- .C(fastkernel,
		ker.r = as.double(ker.r),
		ker.i = as.double(ker.i),
		as.integer(x.min),
		as.integer(x.max),
		as.integer(x.inc),
		as.integer(lng),
		as.double(node),
		as.double(phinode),
		as.integer(nbnode),
		as.double(w0),
		as.double(b.start),
		as.double(b.end))
	ker.r <- z$ker.r
	ker.i <- z$ker.i
	dim(ker.r) <- c(lng, lng)
	dim(ker.i) <- c(lng, lng)
	ker <- matrix(0, lng, lng)
	i <- sqrt(as.complex(-1))
	ker <- ker.r + i * ker.i
	if(plot == T) {
		par(mfrow = c(1, 1))
		image(Mod(ker))
		title("Matrix of the reconstructing kernel (modulus)")
	}
	ker
}}
%\keyword{~keyword}
% Converted by Sd2Rd version 1.21.
