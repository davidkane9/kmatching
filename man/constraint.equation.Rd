\name{constraint.equation}
\alias{constraint.equation}
\title{Creates a constraint equation based on the input data frame.}
\usage{
  constraint.equation(x, weight.var, match.var, replace)
}
\arguments{
  \item{x}{data frame containing needed input data}

  \item{weight.var}{character name of the column of the
  input weights}

  \item{match.var}{character vector of names of columns of
  'data' we wish to match on}

  \item{replace}{logical indicating whether or not
  bservations weighted in the original weight.var are
  allowed positive weight in the output.}
}
\value{
  A list with two components: A and b, representing the
  components of the constraint equation \eqn{Ax = b}
}
\description{
  Creates a constraint equation based on the input data
  frame.
}
\author{
  David Kane \email{<dave.kane@gmail.com>}
}

