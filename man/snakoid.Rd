\name{snakoid}
\alias{snakoid}
\title{
Modified Snake Method
}
\description{
Estimate a ridge from a time-frequency representation, using the
modified snake method (modified cost function).
}
\usage{
snakoid(modulus, guessA, guessB, snakesize=length(guessB),
tfspec=numeric(dim(modulus)[2]), subrate=1, temprate=3, muA=1,
muB=muA, lambdaB=2 * muB, lambdaA=2 * muA, iteration=1000000,
seed=-7, costsub=1, stagnant=20000, plot=TRUE)
}
\arguments{
\item{modulus}{
Time-Frequency representation (real valued).
}
\item{guessA}{
Initial guess for the algorithm (frequency variable).
}
\item{guessB}{
Initial guess for the algorithm (time variable).
}
\item{snakesize}{
The length of the first guess of time variable.
}
\item{tfspec}{
Estimate for the contribution of srthe noise to modulus.
}
\item{subrate}{
Subsampling rate for ridge estimation.
}
\item{temprate}{
Initial value of temperature parameter.
}
\item{muA}{
Coefficient of the ridge's derivative in cost function (frequency
component).
}
\item{muB}{
Coefficient of the ridge's derivative in cost function (time
component).
}
\item{lambdaB}{
Coefficient of the ridge's second derivative in cost function
(time component).
}
\item{lambdaA}{
Coefficient of the ridge's second derivative in cost function
(frequency component).
}
\item{iteration}{
Maximal number of moves.
}
\item{seed}{
Initialization of random number generator.
}
\item{costsub}{
Subsampling of cost function in output.
}
\item{stagnant}{
Maximum number of stationary iterations before stopping.
}
\item{plot}{
when set(default), some results will be displayed
}}
\value{
Returns a structure containing:
\item{ridge}{1D array (of same length as the signal) containing the ridge.}
\item{cost}{1D array containing the cost function.}
\item{plot}{when set(default), some results will be displayed.}
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{corona}}, \code{\link{coronoid}}, \code{\link{icm}},
\code{\link{snake}}.
}
\keyword{ts}
