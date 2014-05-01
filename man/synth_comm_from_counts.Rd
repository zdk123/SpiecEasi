\name{synth_comm_from_counts}
\alias{synth_comm_from_counts}
\title{from count data (ex HMP) fit parameters to OTU margins
and simulate a new community with those properties}
\usage{
synth_comm_from_counts(comm, mar = 2, distr, Sigma = cov(comm), params,
  n = ncol(comm), retParams = FALSE, ...)
}
\arguments{
  \item{comm}{community: matrix of counts}

  \item{distr}{distribution to fit (see fitdistr)}

  \item{Sigma}{covariance structure (defaults to empirical
  cov of comm)}

  \item{params}{optionally supply already fitted
  parameters}

  \item{n}{number of samples (defaults to comm samples)}

  \item{retParams}{if TRUE, return the fitted parameters}

  \item{...}{additional parameters to parameter fitting}
}
\value{
community
}
\description{
from count data (ex HMP) fit parameters to OTU margins and
simulate a new community with those properties
}

