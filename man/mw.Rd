\name{mw}
\alias{mw}
\title{
Dyadic Wavelet Transform
}
\description{
Dyadic wavelet transform, with Mallat's wavelet.
The reconstructed signal preserves locations and values at extrema.
}
\usage{
mw(inputdata, maxresoln, filtername="Gaussian1", scale=FALSE, plot=TRUE)
}
\arguments{
\item{inputdata}{
either a text file or an \R object containing data.
}
\item{maxresoln}{
number of decomposition scales.
}
\item{filtername}{
name of filter (either Gaussian1 for Mallat and Zhong's wavelet or
Haar wavelet).
}
\item{scale}{
when set, the wavelet transform at each scale is plotted 
with the same scale.
}
\item{plot}{
indicate if the wavelet transform at each
scale will be plotted.
}}
\value{
Structure containing
\item{original}{
  original signal.
}
\item{Wf}{
  dyadic wavelet transform of signal.
}
\item{Sf}{
  multiresolution of signal.
}
\item{maxresoln}{
  number of decomposition scales.
}
\item{np}{
  size of signal.
}}
\details{
The decomposition goes from resolution 1 to the given maximum resolution.  
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{dwinverse}}, \code{\link{mrecons}}, \code{\link{ext}}.
}
\keyword{ts}
