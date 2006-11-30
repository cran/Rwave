\name{WV}
\alias{WV}
\title{
  Wigner-Ville function
}
\description{
  Compute the Wigner-Ville transform, without any smoothing.
}
\usage{
WV(input, nvoice, freqstep = (1/nvoice), plot = TRUE)
}
\arguments{
  \item{input}{input signal (possibly complex-valued)}
  \item{nvoice}{number of frequency bands}
  \item{freqstep}{sampling rate for the frequency axis}
  \item{plot}{if set to TRUE, displays the modulus of CWT on the graphic
    device.}
}
\value{
  (complex) Wigner-Ville transform.
}
\details{
}
\references{
  See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
}
\keyword{ts}
