\name{overlap}
\alias{overlap}
\title{Checks the similarity of a group of vectors}
\usage{
  overlap(mat, compare = FALSE)
}
\arguments{
  \item{mat}{The matrix of vectors to compare}

  \item{compare}{Whether to compare to vectors distributed
  uniformly on a simplex scaled by the sum of the rightmost
  vector}
}
\description{
  This function takes a matrix as an argument where the
  columns are vectors we want to compare and calculates the
  correlation and Euclidean distance between these vectors
}
\examples{
data = data.frame(size = rnorm(50), weight = rep(.02, 50))
weights = kmatch(data = data, match.var = "size", weight.var = "weight", n = 100)
overlap(weights)
}
\author{
  Mike Flynn \email{<mflynn210@gmail.com>}
}

