#' Samples from \eqn{Ax=b} and \eqn{w>0}.
#' 
#' Uniformly samples from a convex polytope given by linear equalities in the 
#' parameters using a hit-and-run algorithm. Given constraints: \eqn{Ax = b} and
#' \eqn{x \ge 0} the algorithm finds a point on the interior of the constraints.
#' From there it picks a direction in the k-plane defined by \eqn{Ax = b} and
#' then calculates the maximum and minimum distances (\code{tmin} and
#' \code{tmax}) it can move in that direction. It picks a random
#' distance to travel between \code{tmin} and \code{tmax} and this is used as
#' the next point. This algorithm is useful because each sample is made in
#' constant time.
#' 
#' To find \code{tmin} and \code{tmax}, we must find the first component that becomes zero
#' when going in the negative or positive direction (t). We have that \eqn{x_i + u_it  = 0} or 
#' \eqn{t = =-\frac{x_i}{u_i}}. We must find the minimum and maximum for positive and negative t, respectively, in i.
#' Once those bounds are found, a value of t is picked uniformly on the interval between \code{tmin} and \code{tmax}.
#' 
#' @param A Matrix of constraint coefficients, rows should correspond to each 
#'   constraint. A must not have collinear rows
#' @param b A vector corresponding to the right hand side of the constraints
#' @param n The number of output vectors desired
#' @param discard A burninlength, how many vectors should be discarded before 
#'   recording
#' @param skiplength Only 1 out of every 'skiplength' vectors will be recorded
#' @param chains number of different chains, starting from different starting points
#' @param verbose Give verbose output of how the function is progressing.
#' @param achr Whether to use "accelerated convergence hit and run" algorithm proposed in
#' paper: <Insert Kaufman, Smith citation>. "discard" will be used to collect presamples to
#' estimate the expected value of the span.
#' 
#' @return Gives back a list of matrices with 'n' columns corrresponding
#' to n uniformly sampled solutions of Ax = b. The number of lists = "chains" variable. 
#' 'n' columns.  
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#'   
#' @examples
#' A <- matrix(1, ncol = 5, nrow = 1)
#' b <- 1
#' w <- hitandrun(A, b, n = 100)
#' 
#' A <- matrix(1, ncol = 100, nrow = 1)
#' b <- 50
#' A <- rbind(A, rnorm(100))
#' b <- c(b,0)
#' w <- hitandrun(A, b, n = 100) 
#' 
#' ##2 chains
#' chains.2 <- hitandrun(A, b, n = 10, chains = 2)

hitandrun <- function(A, b, n, discard = 1000, skiplength = 5, chains = 1, verbose = FALSE, achr = TRUE) {
    
    if(n <= 0 || n %% 1 != 0) {
      stop("n must be a positive integer")
    }
    
    if(any(is.na(A))) 
      stop(paste("'A' cannot have NA's in it, has ", sum(is.na(A)), sep = ""))
  
    
    if(any(is.na(b)))
      stop(paste("'b' cannot have NA's in it, has ", sum(is.na(b)), sep = ""))
    
    ## overdetermined, more constraints than degrees of freedom
    dimen = dim(A)
    if(dimen[1] > dimen[2]) {  
      stop("Problem is overdetermined, more constraints than degrees of freedom")
    }
    
    #unique solution
    if(dimen[1] == dimen[2]) {
      if(det(A) != 0) {
        warning("Solution to Ax=b is unique")
        sol = solve(A) %*% b
        if(all(sol > 0)) {
          return(list(sol))
        } else {
          stop("There is no positive solution")
        }
      }
    }
    
    ## resolve weird quirk in Null() function
    if(ncol(A) ==1) {
      Z <- Null(A)
    } else {
      Z <- Null(t(A))
    }
    
    
    chainlist <- 
    llply(1:chains, function(chainnum) {
      str <- "Finding an intial solution..."
      if(verbose) cat(str)
      X <- matrix(0, nrow = ncol(A), ncol = (n + discard))
      index <- 1
      ## make initial solution = Ap %*% b where 
      ## Ap is the psuedoinverse of A
      ## We want to get to solution x of Ax = b
      ## if A was invertible we could just use
      ## x = Ai*b, however, A is not square and therefore
      ## not invertible. Luckily, there exists an object
      ## called the psuedoinverse Ap such that Ap * A * b = b.
      ## Ap can be constructed from the SVD of A, where
      ## SVD(A) = U D V', Ap = V (1/D)' U'. 1/D simply means take
      ## 1/d for each non-zero entry d, and zero all the others.
      SVD <- svd(A)
      d <- SVD[['d']]
      V <- SVD[['v']]
      U <- SVD[['u']]
      ## get rid of division errors in d because they will
      ## mess up 1/d (we keep doing 1/0, so R gives us Inf as well
      ## we must zero these)
      d[d < 1e-10] <- 0
      di <- 1/d
      di[di == Inf] <- 0
      if(length(di) <= 1) Dt <- matrix(di, ncol = 1, nrow =1) else Dt <- t(diag(di))
      Ap <- V %*% Dt %*% t(U)
      ## l is initial solution
      ## http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Obtaining_all_solutions_of_a_linear_system
      l <- Ap %*% b + (diag(nrow(Ap)) - Ap %*% A) %*% rnorm(nrow(Ap))
      
      ## if l isn't in feasible space use mirror algorithm to find viable solution
      ## in the interior
      if(!(all(l > 0))) {
        if(verbose) for(i in 1:nchar(str)) cat("\b")
        str <- "Using mirror algorithm to find inner solution...\n"
        if(verbose) cat(str)
        y <- mirror(A, l, 1, verbose)
      } else {
        y <- l
      }
      
      if(verbose) for(i in 1:nchar(str)) cat("\b")
      return(hnr_loop(y, Z, n, skiplength, discard, achr))
    }, .parallel = FALSE, .paropts = list(.packages = "kmatching", .export = "hnr_loop"))
    return(chainlist)
}
