\name{cwtTh}
\alias{cwtTh}
\title{Cauchy's wavelet transform}
\description{
  Compute the continuous wavelet transform with (complex-valued)
  Cauchy's wavelet.
}
\usage{cwtTh(input, noctave, nvoice=1, moments, twoD=TRUE, plot=TRUE)
}
\arguments{
  \item{input}{
    input signal (possibly complex-valued).
  }
  \item{noctave}{
    number of powers of 2 for the scale variable.
  }
  \item{nvoice}{
    number of scales in each octave
    (i.e. between two consecutive powers of 2).
  }
  \item{moments}{
    number of vanishing moments.
  }
  \item{twoD}{
    logical variable set to \code{T} to organize the output as a 2D array
    (signal size x nb scales), otherwise, the output is a 3D array (signal
    size x noctave x nvoice).
  }
  \item{plot}{
    if set to \code{T}, display the modulus of the continuous wavelet
    transform on the graphic device.
}}
\value{
  \item{tmp}{continuous (complex) wavelet transform.}
}
\details{
  The output contains the (complex) values of the wavelet transform of
  the input signal.  The format of the output can be

  2D array (signal size \eqn{\times}{x} nb scales)

  3D array (signal size \eqn{\times}{x} noctave \eqn{\times}{x} nvoice)
}
\references{
  See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\seealso{
  \code{\link{cwt}}, \code{\link{cwtp}}, \code{\link{DOG}},
  \code{\link{gabor}}.
}
\keyword{ts}
