\name{cleanph}
\alias{cleanph}
\title{
Threshold Phase based on Modulus
}
\description{
Sets to zero the phase of time-frequency transform when
modulus is below a certain value.
}
\usage{
cleanph(tfrep, thresh=0.01, plot=TRUE)
}
\arguments{
\item{tfrep}{
continuous time-frequency transform (2D array)
}
\item{thresh}{
(relative) threshold.
}
\item{plot}{
if set to TRUE, displays the maxima of cwt on the graphic
device.
}}
\value{
thresholded phase (2D array)
}
\details{
}
\references{
See discussion in text of ``Practical Time-Frequency Analysis''.
}
\seealso{
}
\keyword{ts}
