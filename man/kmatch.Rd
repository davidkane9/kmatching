\name{kmatch}
\alias{kmatch}
\title{Create matching vectors based on data Creates a matrix of weights which}
\usage{
  kmatch(x, weight.var, match.var, n = 1, replace = FALSE,
    ...)
}
\arguments{
  \item{x}{data frame containing needed input data}

  \item{weight.var}{character name of the column of the
  input weights}

  \item{match.var}{character vector of names of columns of
  'data' we wish to match on}

  \item{n}{numeric number of weight vectors desired.
  Default is 1}

  \item{replace}{logical indicating whether or not
  bservations weighted in the original weight.var are
  allowed positive weight in the output. Default is FALSE}

  \item{...}{parameters to be passed to the sampling
  methods}
}
\value{
  Returns a matrix of 'n' sets of weights that match the
  given set of weights in terms of weighted averages to the
  'match.var' factors. The columns are the sets of weights.
}
\description{
  Create matching vectors based on data Creates a matrix of
  weights which
}
\examples{
set.seed(40)
x <- data.frame(size = rnorm(50), weight = rep(.02, 50))
weights <- kmatch(x, weight.var = "weight", match.var = "size", n = 100, replace = TRUE)
}
\author{
  Mike Flynn \email{mflynn210@gmail.com}
}

