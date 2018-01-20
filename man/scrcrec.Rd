\name{scrcrec}
\alias{scrcrec}
\title{
Simple Reconstruction from Crazy Climbers Ridges
}
\description{
Reconstructs signal from ridges obtained by \code{\link{crc}},
using the restriction of the transform to the ridge.
}
\usage{
scrcrec(siginput, tfinput, beemap, bstep=5, ptile=0.01, plot=2)
}
\arguments{
\item{siginput}{
input signal.
}
\item{tfinput}{
time-frequency representation (output of \code{\link{cwt}}
or \code{\link{cgt}}.
}
\item{beemap}{
output of crazy climber algorithm
}
\item{bstep}{
used for the chaining (see \code{\link{cfamily}}).
}
\item{ptile}{
threshold on the measure beemap (see \code{\link{cfamily}}).
}
\item{plot}{
  1: displays signal,components, and reconstruction one after another.\cr
  2: displays signal, components and reconstruction.\cr
  Else, no plot.
}}
\value{
Returns a list containing the reconstructed signal and the chained ridges.
\item{rec}{reconstructed signal}
\item{ordered}{image of the ridges (with different colors)}
\item{comp}{2D array containing the signals reconstructed from ridges}
}
%\details{}
\references{
See discussions in the text of \dQuote{Practical Time-Frequency Analysis}.
}
\seealso{
\code{\link{crc}},\code{\link{cfamily}} for crazy climbers method,
\code{\link{crcrec}} for reconstruction methods.
}
\keyword{ts}
