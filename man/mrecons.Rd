\name{mrecons}
\alias{mrecons}
\title{
Reconstruct from Dyadic Wavelet Transform Extrema
}
\description{
Reconstruct from dyadic wavelet transform modulus extrema.
The reconstructed signal preserves locations and values at extrema.
}
\usage{
mrecons(extrema, filtername="Gaussian1", readflag=TRUE)
}
\arguments{
\item{extrema}{
the extrema representation. 
}
\item{filtername}{
filter used for dyadic wavelet transform.
}
\item{readflag}{
if set to T, read reconstruction kernel from precomputed file.
}}
\value{
Structure containing
\item{f}{
  the reconstructed signal.
}
\item{g}{
  reconstructed signal plus mean of original signal.
}
\item{h}{
  reconstructed signal plus coarse scale component of original signal.
}}
\details{
The reconstruction involves only the wavelet coefficients,
without taking care of the coarse scale component. The latter
may be added a posteriori.
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{mw}}, \code{\link{ext}}.
}
\keyword{ts}
