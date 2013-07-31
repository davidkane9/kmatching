\name{hitandrun}
\alias{hitandrun}
\title{Samples from \eqn{Ax=b} and \eqn{x>0}.}
\usage{
  hitandrun(A, b, n, discard = 0, skiplength = 5,
    chains = 1, verbose = FALSE)
}
\arguments{
  \item{A}{Matrix of constraint coefficients, rows should
  correspond to each constraint. A must not have collinear
  rows}

  \item{b}{A vector corresponding to the right hand side of
  the constraints}

  \item{n}{The number of output vectors desired}

  \item{discard}{A burninlength, how many vectors should be
  discarded before recording}

  \item{skiplength}{Only 1 out of every 'skiplength'
  vectors will be recorded}

  \item{chains}{number of different chains, starting from
  different starting points}

  \item{verbose}{Give verbose output of how the function is
  progressing.}
}
\value{
  In the even chains = 1, gives back a matrix with 'n'
  columns corrresponding to n solutions of Ax = b, if
  chains > 1, gives back a list of 'chains' matrices with
  'n' columns.
}
\description{
  Randomly samples from a convex polytope given by linear
  equalities in the parameters using a hit-and-run
  algorithm. Given constraints: \eqn{Ax = b} and \eqn{x \ge
  0} the algorithm finds a point on the interior of the
  constraints. From there it picks a direction in the
  k-plane defined by \eqn{Ax = b} and then calculates the
  maximum and minimum distances (\code{tmin} and
  \code{tmax}) it can move in that direction. It picks a
  random distance to travel between \code{tmin} and
  \code{tmax} and this is used as the next point. This
  algorithm is useful because each sample is made in
  constant time.
}
\details{
  To find \code{tmin} and \code{tmax}, we must find the
  first component that becomes zero when going in the
  negative or positive direction (t). We have that \eqn{x_i
  + u_it = 0} or \eqn{t = =-\frac{x_i}{u_i}}. We must find
  the minimum and maximum for positive and negative t,
  respectively, in i. Once those bounds are found, a value
  of t is picked uniformly on the interval between
  \code{tmin} and \code{tmax}.
}
\examples{
A <- matrix(1, ncol = 5, nrow = 1)
b <- 1
w <- hitandrun(A, b, n = 100)

A <- matrix(1, ncol = 100, nrow = 1)
b <- 50
A <- rbind(A, rnorm(100))
b <- c(b,0)
w <- hitandrun(A, b, n = 100)

##2 chains
chains.2 <- hitandrun(A, b, n = 10, chains = 2)
}
\author{
  Mike Flynn \email{mflynn210@gmail.com}
}

