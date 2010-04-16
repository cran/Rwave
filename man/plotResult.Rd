\name{plotResult}
\alias{plotResult}
\title{
Plot Dyadic Wavelet Transform Extrema
}
\description{
Plot extrema of dyadic wavelet transform.
}
\usage{
plotResult(result, original, maxresoln, scale=FALSE, yaxtype="s")
}
\arguments{
  \item{result}{
    result.
  }
  \item{original}{
    input signal.
  }
  \item{maxresoln}{
    number of decomposition scales.
  }
  \item{scale}{
    when set, the extrema at each scale is plotted withe the same scale.
  }
  \item{yaxtype}{
    y axis type (see \R manual).
  }
}
%\value{}
%\details{}
\references{
  See discussions in the text of ``Time-Frequency Analysis''.
}
\seealso{
  \code{\link{plotwt}}, \code{\link{epl}}, \code{\link{wpl}}.
}
\keyword{ts}
