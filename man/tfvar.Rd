\name{tfvar}
\alias{tfvar}
\title{
Variance frequency by frequency
}
\description{
Compute the variance of time-frequency representation frequency 
by frequency.
}
\usage{
tfvar(input, plot=TRUE)
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
\code{\link{tfmean}},\code{\link{tfpct}}.
}
\keyword{ts}
