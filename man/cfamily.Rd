\name{cfamily}
\alias{cfamily}
\title{
Ridge Chaining Procedure
}
\description{
Chains the ridge estimates produced by the function \code{\link{crc}}.
}
\usage{
cfamily(ccridge, bstep=1, nbchain=100, ptile=0.05)
}
\arguments{
\item{ccridge}{
unchained ridge set as the output of the function \code{\link{crc}}
}
\item{bstep}{
maximal length for a gap in a ridge.
}
\item{nbchain}{
maximal number of chains produced by the function.
}
\item{ptile}{
relative threshold for the ridges.
}}
\value{
  Returns the results of the chaining algorithm
  \item{ordered map}{
    image containing the ridges (displayed with different colors)
  }
  \item{chain}{
    2D array containing the chained ridges, according to the chain data
    structure\cr
    chain[,1]: first point of the ridge\cr
    chain[,2]: length of the chain\cr
    chain[,3:(chain[,2]+2)]: values of the ridge\cr
  }
  \item{nbchain}{
    number of chains produced by the algorithm
}}
\details{
  \code{\link{crc}} returns a measure in time-frequency (or time-scale)
  space. \code{\link{cfamily}} turns it into a series of one-dimensional
  objects (ridges). The measure is first thresholded, with a relative
  threshold value set to the input parameter ptile. During the chaining
  procedure, gaps within a given ridge are allowed and filled in. The
  maximal length of such gaps is the input parameter bstep.
}
\references{
  See discussion in text of ``Practical Time-Frequency Analysis''.
}
\seealso{
  \code{\link{crc}} for the ridge estimation, and \code{\link{crcrec}},
  \code{\link{gcrcrec}} and \code{\link{scrcrec}} for corresponding
  reconstruction functions.
}
\keyword{ts}
