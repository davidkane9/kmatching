\name{kmatch}
\alias{kmatch}
\title{Create matching vectors based on data}
\usage{
  kmatch(data, match.var, weight.var, sumlimit = NULL, ...)
}
\arguments{
  \item{data}{A data frame containing data we want to
  match}

  \item{match.var}{A list of names of columns of 'data' we
  wish to match}

  \item{weight.var}{The name of the column of orignal set
  of weights.}

  \item{sumlimit}{An optional parameter specifying what the
  weights in the vector should sum to. This cannot be used
  if there are discrete matching variables because it ends
  up being redundant or contradictory.}

  \item{...}{Parameters to be passed to the sampling
  methods}
}
\description{
  An abstraction on top of the random sampling methods.
  This function allows a user to specify which variables
  they would like to match exposures to an original vector
  of weights.
}
\examples{
data = data.frame(size = rnorm(50), weight = rep(.02, 50))
weights = kmatch(data = data, match.var = "size", weight.var = "weight", n = 100)
}
\author{
  Mike Flynn \email{mjf2@williams.edu}
}

