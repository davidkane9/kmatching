\name{mirror}
\alias{mirror}
\title{Generates weights to a new portfolio using mirror algorithm
Fulfills equality constraints while maintaining randomness by
using a Monte Carlo Random Walk reflecting at the boundaries. Based
on xsample() function in limSolve package.}
\usage{
  mirror(Emat, x0, n, verbose = FALSE, numjump = 20)
}
\arguments{
  \item{Emat}{This is the matrix of the equality constraint
  coefficients}

  \item{x0}{An original solution to the constraints}

  \item{n}{Number of random solutions to output}

  \item{verbose}{Give verbose output describing the
  progress of the function}

  \item{numjump}{The number of jumps to scatter around the
  direction given by the difference from zero}
}
\description{
  Generates weights to a new portfolio using mirror
  algorithm Fulfills equality constraints while maintaining
  randomness by using a Monte Carlo Random Walk reflecting
  at the boundaries. Based on xsample() function in
  limSolve package.
}
\examples{
Emat = matrix(1, ncol = 3, nrow = 1)
x0 = c(.3, .3, .4)
mirror(Emat, x0, 1)
}
\author{
  Mike Flynn \email{<mflynn210@gmail.com>}
}
\references{
  Van Den Meershe, Karel, Karline Soetaert, and Dick Van
  Oevelen. "Xsample(): An R Function for Sampling Linear
  Inverse Problems." Journal of Statistical Software 30
  (2009): 1-15. Print.
  \url{http://cran.cermin.lipi.go.id/web/packages/limSolve/vignettes/xsample.pdf}
}

