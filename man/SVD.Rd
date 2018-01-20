\name{SVD}
\alias{SVD}
\title{
Singular Value Decomposition
}
\description{
Computes singular value decomposition of a matrix.
}
\usage{
SVD(a)
}
\arguments{
\item{a}{
input matrix.
}}
\value{
a structure containing the 3 matrices of the singular value
decomposition of the input.
}
\details{
\R interface for Numerical Recipes singular value decomposition routine.
}
\references{
See discussions in the text of \dQuote{Time-Frequency Analysis}.
}
%\seealso{}
\examples{
   hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
     X <- hilbert(6)
     z  = SVD(X)
     z
}
\keyword{ts}
