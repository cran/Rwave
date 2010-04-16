\name{dwinverse}
\alias{dwinverse}
\title{
Inverse Dyadic Wavelet Transform
}
\description{
Invert the dyadic wavelet transform.
}
\usage{
dwinverse(wt, filtername="Gaussian1")
}
\arguments{
\item{wt}{
dyadic wavelet transform
}
\item{filtername}{
filters used. ("Gaussian1" stands for the filters
corresponds to those of Mallat and Zhong's wavlet. And
"Haar" stands for the filters of Haar basis.	
}}
\value{
Reconstructed signal
}
%\details{}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{mw}}, \code{\link{ext}}, \code{\link{mrecons}}.
}
\keyword{ts}
