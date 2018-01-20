\name{tfgmax}
\alias{tfgmax}
\title{
Time-Frequency Transform Global Maxima
}
\description{
Computes the maxima (for each fixed value of the time variable)
of the modulus of a continuous wavelet transform.
}
\usage{
tfgmax(input, plot=TRUE)
}
\arguments{
\item{input}{
wavelet transform (as the output of the function \code{\link{cwt}})
}
\item{plot}{
if set to TRUE, displays the values of the energy as
a function of the scale.
}}
\value{
\item{output}{values of the maxima (1D array)}
\item{pos}{positions of the maxima (1D array)}
}
%\details{}
\references{
See discussions in the text of \dQuote{Practical Time-Frequency Analysis}.
}
\seealso{
\code{\link{tflmax}}.
}
\keyword{ts}
