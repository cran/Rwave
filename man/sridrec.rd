\name{sridrec}
\alias{sridrec}
\title{
Simple Reconstruction from Ridge
}
\description{
Simple reconstruction of a real valued signal from a ridge,
by restriction of the transform to the ridge.
}
\usage{
sridrec(tfinput, ridge)
}
\arguments{
\item{tfinput}{
time-frequency representation.
}
\item{ridge}{
ridge (1D array).
}}
\value{
(real) reconstructed signal (1D array)
}
%\details{}
\references{
See discussions in the text of \dQuote{Practical Time-Frequency Analysis}.
}
\seealso{
\code{\link{ridrec}}, \code{\link{gridrec}}.
}
\keyword{ts}
