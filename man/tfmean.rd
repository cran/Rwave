\name{tfmean}
\alias{tfmean}
\title{
Average frequency by frequency
}
\description{
Compute the mean of time-frequency representation frequency 
by frequency.
}
\usage{
tfmean(input, plot=TRUE)
}
\arguments{
\item{input}{
time-frequency transform (output of \code{\link{cwt}} or \code{\link{cgt}}).
}
\item{plot}{
if set to T, displays the values of the energy as a 
function of the scale (or frequency).
}}
\value{
1D array containing the noise estimate.
}
%\details{}
\references{
See discussions in the text of \dQuote{Practical Time-Frequency Analysis}.
}
\seealso{
\code{\link{tfpct}},\code{\link{tfvar}}.
}
\keyword{ts}
