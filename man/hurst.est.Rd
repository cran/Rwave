\name{hurst.est}
\alias{hurst.est}
\title{
Estimate Hurst Exponent
}
\description{
Estimates Hurst exponent from a wavelet transform.
}
\usage{
hurst.est(wspec, range, nvoice, plot=TRUE)
}
\arguments{
\item{wspec}{
wavelet spectrum (output of \code{\link{tfmean}})
}
\item{range}{
range of scales from which estimate the exponent.
}
\item{nvoice}{
number of scales per octave of the wavelet transform.
}
\item{plot}{
if set to \code{TRUE}, displays regression line on current plot.
}}
\value{
complex 1D array of size sigsize.
}
\details{
}
\references{
See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
\code{\link{tfmean}}, \code{\link{wspec.pl}}.
}
\keyword{ts}
\examples{
# White Noise Hurst Exponent: The plots on the top row of Figure 6.8
# were produced by the folling S-commands. These make use of the two
# functions Hurst.est (estimation of Hurst exponent from CWT) and
# wspec.pl (display wavelet spectrum).

# Compare the periodogram and the wavelet spectral estimate.
wnoise <- rnorm(8192)
plot.ts(wnoise)
spwnoise <- fft(wnoise)
spwnoise <- Mod(spwnoise)
spwnoise <- spwnoise*spwnoise
plot(spwnoise[1:4096], log="xy", type="l")
lswnoise <- lsfit(log10(1:4096), log10(spwnoise[1:4096]))
abline(lswnoise$coef)
cwtwnoise <- DOG(wnoise, 10, 5, 1, plot=FALSE)
mcwtwnoise <- Mod(cwtwnoise)
mcwtwnoise <- mcwtwnoise*mcwtwnoise
wspwnoise <- tfmean(mcwtwnoise, plot=FALSE)
wspec.pl(wspwnoise, 5)
hurst.est(wspwnoise, 1:50, 5)
}
