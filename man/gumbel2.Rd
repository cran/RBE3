\name{gumbel2}
\alias{dgumbel2}
\alias{pgumbel2}
\alias{qgumbel2}
\title{The Gumbel2 distribution
}
\description{
Density, distribution function and quantile function for the Gumbel2 distribution.
}
\usage{
dgumbel2(x, log=FALSE)
pgumbel2(q)
qgumbel2(p)
}
\arguments{
  \item{x, q}{
Vector of quantiles.
}
\item{p}{
Vector of probabilities.
}
\item{log}{
logical; if TRUE, probabilities p are given as log(p).
}
}
\details{
The cumulative distribution function for the Gumbel2 distribution is given by \eqn{F(x)=1-exp(-exp(x)).}
}
\value{
dgumbel2 gives the density, pgumbel2 gives the distribution function and qgumbel2 gives the quantile function.

The length of the result is determined by the maximum of the lengths of the numerical arguments.
}

\author{
Diego Gallardo and Marcelo Bourguignon.
}
\examples{
dgumbel2(c(4,10))
pgumbel2(c(4,10))
qgumbel2(c(0.1, 0.5))
}


