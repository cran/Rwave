\name{tflmax}
\alias{tflmax}
\title{
Time-Frequency Transform Local Maxima
}
\description{
Computes the local maxima (for each fixed value of the time variable)
of the modulus of a time-frequency transform.
}
\usage{
tflmax(input, plot=TRUE)
}
\arguments{
\item{input}{
time-frequency transform (real 2D array).
}
\item{plot}{
 if set to T, displays the local maxima on the graphic device.
}}
\value{
values of the maxima (2D array).
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{tfgmax}}.
}
\keyword{ts}
