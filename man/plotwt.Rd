\name{plotwt}
\alias{plotwt}
\title{
Plot Dyadic Wavelet Transform
}
\description{
Plot dyadic wavelet transform.
}
\usage{
plotwt(original, psi, phi, maxresoln, scale=FALSE, yaxtype="s")
}
\arguments{
\item{original}{
input signal.
}
\item{psi}{
dyadic wavelet transform.
}
\item{phi}{
scaling function transform at last resolution.
}
\item{maxresoln}{
number of decomposition scales.
}
\item{scale}{
when set, the wavelet transform at each scale is plotted with the same
scale.
}
\item{yaxtype}{
 axis type (see \R manual).
}}
\value{

}
\details{
}
\references{
See discussions in the text of ``Time-Frequency Analysis''.
}
\seealso{
\code{\link{plotResult}}, \code{\link{epl}}, \code{\link{wpl}}.
}
\keyword{ts}
