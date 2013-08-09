\name{kmatch}
\alias{kmatch}
\title{Create matching vectors based on data}
\usage{
  kmatch(x, weight.var, match.var, n = 1, chains = 1,
    replace = FALSE, verbose = FALSE, skiplength = 5, ...)
}
\arguments{
  \item{x}{data frame containing needed input data}

  \item{weight.var}{character name of the column of the
  input weights}

  \item{match.var}{character vector of names of columns of
  'data' we wish to match on}

  \item{n}{numeric number of weight vectors desired.
  Default is 1}

  \item{chains}{number of different chains, starting from
  different starting points. Default is 1}

  \item{replace}{if FALSE then any rows that had non-zero
  weight in 'weight.var' will be set to 0 in the output}

  \item{verbose}{TRUE to give verbose output, including
  gelman-rubin analysis on the random walk}

  \item{...}{parameters to be passed to the sampling
  methods}
}
\value{
  Returns a list of "chains" matrices of 'n' sets of
  weights that match the given set of weights in terms of
  weighted averages to the 'match.var' factors. The columns
  of each matrix are the sets of weights.
}
\description{
  Creates a matrix of weights which match a given set of
  weights in exposures to designated factors
}
\examples{
set.seed(40)
x <- data.frame(size = rnorm(50), weight = rep(.02, 50))
weights <- kmatch(x, weight.var = "weight", match.var = "size", n = 100, replace = TRUE)

x2 <- data.frame(var1 = c(-2,-1, 0, 1, 2), var2 = c(1, -1, 2, -1, 0), w = rep(.2, 5))
matchvars = c("var1", "var2")
weights2 = kmatch(x=x2, weight.var = "w", match.var = matchvars, n= 1000, skiplength = 100, replace= TRUE)[[1]]
## for interactive graphics:
## library(rgl)
## plot3d(x = weights2[1,], y = weights2[2,], z = weights2[3,])
}
\author{
  Mike Flynn \email{mflynn210@gmail.com}
}

