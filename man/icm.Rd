\name{icm}
\alias{icm}
\title{
Ridge Estimation by ICM Method
}
\description{
Estimate a (single) ridge from a time-frequency representation,
using the ICM minimization method.
}
\usage{
icm(modulus, guess, tfspec=numeric(dim(modulus)[2]), subrate=1,
mu=1, lambda=2 * mu, iteration=100)
}
\arguments{
\item{modulus}{
Time-Frequency representation (real valued).
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
\item{mu}{
Coefficient of the ridge's second derivative in cost function.
}
\item{lambda}{
Coefficient of the ridge's derivative in cost function.
}
\item{iteration}{
Maximal number of moves.
}}
\value{
Returns the estimated ridge and the cost function.
\item{ridge}{
  1D array (of same length as the signal) containing the ridge.
}
\item{cost}{
  1D array containing the cost function.
}}
\details{
To accelerate convergence, it is useful to preprocess modulus before
running annealing method. Such a preprocessing (smoothing and
subsampling of modulus) is implemented in \code{\link{icm}}. The
parameter subrate specifies the subsampling rate.
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{corona}}, \code{\link{coronoid}}, and \code{\link{snake}},
\code{\link{snakoid}}.
}
\keyword{ts}
