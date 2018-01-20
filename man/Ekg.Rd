\name{Ekg}
\alias{Ekg}
\title{Heart Rate Data}
\description{Successive beat-to-beat intervals for a normal patient.
}
\usage{data(Ekg)
}
\format{A vector containing 16,042 observations.
}
\source{
  See discussions in the text of ``Practical Time-Frequency Analysis''.
}
\references{
  Carmona, R. A., W. L. Hwang and B Torresani (1998, eBook ISBN:978008053942)
  \emph{Practical Time-Frequency Analysis: Gabor and Wavelet Transforms
    with an Implementation in S}, Academic Press, San Diego.
}
\examples{
data(Ekg)
plot.ts(Ekg)
}
\keyword{datasets}
