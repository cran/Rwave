\name{tfpct}
\alias{tfpct}
\title{
Percentile frequency by frequency
}
\description{
Compute a percentile of time-frequency representation frequency 
by frequency.
}
\usage{
tfpct(input, percent=0.8, plot=TRUE)
}
\arguments{
\item{input}{
time-frequency transform (output of \code{\link{cwt}} or \code{\link{cgt}}).
}
\item{percent}{
percentile to be retained.
}
\item{plot}{
if set to T, displays the values of the energy as a 
function of the scale (or frequency).
}}
\value{
1D array containing the noise estimate.
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{tfmean}},\code{\link{tfvar}}.
}
\keyword{ts}
