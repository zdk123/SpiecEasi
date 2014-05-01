\name{rzipois}
\alias{rzipois}
\title{Draw samples from a zero-inflated poisson distribution}
\usage{
rzipois(n, lambda, pstr0 = 0)
}
\arguments{
  \item{n}{the number of samples to draw}

  \item{lambda}{The poisson rate parameter}

  \item{pstr0}{probability of drawing a zero}
}
\value{
Poisson counts of length \eqn{n}
}
\description{
Draw samples from a zero-inflated poisson distribution
}

