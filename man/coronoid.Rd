\name{coronoid}
\alias{coronoid}
\title{
Ridge Estimation by Modified Corona Method
}
\description{
Estimate a ridge using the
modified corona method (modified cost function).
}
\usage{
coronoid(tfrep, guess, tfspec=numeric(dim(tfrep)[2]), subrate=1,
temprate=3, mu=1, lambda=2 * mu, iteration=1000000, seed=-7,
stagnant=20000, costsub=1, plot=TRUE)
}
\arguments{
\item{tfrep}{
Estimate for the contribution of the noise to modulus.
}
\item{guess}{
Initial guess for the algorithm.
}
\item{tfspec}{
Estimate for the contribution of the noise to modulus.
}
\item{subrate}{
Subsampling rate for ridge estimation.
}
\item{temprate}{
Initial value of temperature parameter.
}
\item{mu}{
Coefficient of the ridge's derivative in cost function.
}
\item{lambda}{
Coefficient of the ridge's second derivative in cost function.
}
\item{iteration}{
Maximal number of moves.
}
\item{seed}{
Initialization of random number generator.
}
\item{stagnant}{
Maximum number of stationary iterations before stopping.
}
\item{costsub}{
Subsampling of cost function in output.
}
\item{plot}{
When set(default), some results will be shown on the display.
}}
\value{
Returns the estimated ridge and the cost function.
\item{ridge}{
  1D array (of same length as the signal) containing the ridge.
}
\item{cost}{
  1D array containing the cost function.
}
}
\details{
To accelerate convergence, it is useful to preprocess modulus before
running annealing method. Such a preprocessing (smoothing and
subsampling of modulus) is implemented in \code{\link{coronoid}}. The
parameter subrate specifies the subsampling rate.
}
\references{
See discussion in text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{corona}}, \code{\link{icm}}, \code{\link{snake}}, \code{\link{snakoid}}.
}
\keyword{ts}
\section{Warning}{
The returned cost may be a large array.
The argument costsub allows subsampling the cost function.
}
