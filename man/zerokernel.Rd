\name{zerokernel}
\alias{zerokernel}
\title{
Reconstruction from Wavelet Ridges
}
\description{
Generate a zero kernel for reconstruction from ridges.
}
\usage{
zerokernel(x.inc=1, x.min, x.max)
}
\arguments{
\item{x.min}{
minimal value of x for the computation of \eqn{Q_2}.
}
\item{x.max}{
maximal value of x for the computation of \eqn{Q_2}.
}
\item{x.inc}{
step unit for the computation of the kernel.
}}
\value{
matrix of the \eqn{Q_2} kernel
}
%\details{}
%\references{}
\seealso{
\code{\link{kernel}}, \code{\link{fastkernel}}, \code{\link{gkernel}},
\code{\link{gkernel}}.
}
\keyword{ts}
