\name{hitandrun}
\alias{hitandrun}
\title{Uniformly samples from {A*x=A*x0}, {x>0} or {Ax=b} & {x>0}}
\usage{
  hitandrun(A, b = NULL, x0 = NULL, n, discard = 0,
    skiplength = 5, verbose = FALSE)
}
\arguments{
  \item{A}{Matrix of constraint coefficients, rows should
  correspond to each constraint. A must not have collinear
  rows.}

  \item{b}{A vector of exposures that correspond to the
  right hand side of the constraints. Should not be used at
  the same time as x0}

  \item{x0}{An original solution we want to match, should
  not be used at the same time as b.}

  \item{n}{The number of output vectors desired}

  \item{discard}{A burninlength, how many vectors should be
  discarded before recording}

  \item{skiplength}{Only 1 out of every 'skiplength'
  vectors will be recorded, in order to optimally spread
  out the output}

  \item{verbose}{Give verbose output of how the function is
  progressing.}
}
\description{
  Randomly samples uniformly from a convex polytope given
  by linear equalities in the parameters. Uses a
  hit-and-run algorithm. Given constraints: \eqn{Ax = b,
  x\ge 0} or \eqn{Ax = Ax_0, x>0} the algorithm finds a
  point on the interior of the constraints. From there it
  picks a direction in the k-plane defined by $Ax=b$, it
  then calculates the maximum and minimum distances it can
  move the point in that direction, called \verb+tmin+ and
  \verb+tmax+. It pick a random distance to travel between
  \verb+tmin+ and \verb+tmax+ and this is used as the next
  point. This algorithm is useful because each sample is
  made in constant time.
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
}
\author{
  Mike Flynn \email{mflynn210@gmail.com}
}

