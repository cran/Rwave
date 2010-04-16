\name{snakeview}
\alias{snakeview}
\title{
Restriction to a Snake
}
\description{
Restrict time-frequency transform to a snake.
}
\usage{
snakeview(modulus, snake)
}
\arguments{
\item{modulus}{
Time-Frequency representation (real valued).
}
\item{snake}{
Time and frequency components of a snake.
}}
\value{
2D array containing the restriction of the transform modulus
to the snake.
}
\details{
Recall that a snake is a (two components) \R structure.
}
\references{
See discussions in the text of \dQuote{Time-Frequency Analysis}.
}
%\seealso{}
\keyword{ts}
