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

hitandrun <- function(A, b, n, discard = 0, skiplength = 5, chains = 1, verbose = FALSE, cpp= FALSE) {
    
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

      if(verbose) cat(paste("Random Walk: Chain ", chainnum, "\nDone with: "))
      str <- "0"
      if(verbose) cat(str)
      if(cpp) return(hnr_loop(y, Z, n, skiplength, discard))
      for(i in 1:(n*skiplength+discard)) {
        tmin<-0;tmax<-0;
        ## runs counts how many times tried to pick a direction, if
        ## too high fail.
        runs = 0
        while(tmin ==0 && tmax ==0) {
          ## r is a random unit vector in with basis in Z
          r <- rnorm(ncol(Z))
          r <- r/sqrt(sum(r^2))
          
          ## u is a unit vector in the appropriate k-plane pointing in a
          ## random direction Z %*% r is the same as in mirror
          u <- Z%*%r
          c <- y/u
          ## determine intersections of x + t*u with walls
          ## the limits on how far you can go backward and forward
          ## i.e. the maximum and minimum ratio y_i/u_i for negative and positive u.
          tmin <- max(-c[u>0]); tmax <- min(-c[u<0]);
          ## unboundedness
          if(tmin == -Inf || tmax == Inf){
            stop("problem is unbounded")
          }
          ## if stuck on boundary point
          if(tmin==0 && tmax ==0) {
            runs = runs + 1
            if(runs >= 1000) stop("hitandrun found can't find feasible direction, cannot generate points")
          }
        }
        
        ## chose a point on the line segment
        y <- y + (tmin + (tmax - tmin)*runif(1))*u;
        
        ## choose a point every 'skiplength' samples
        if(i %% skiplength == 0) {
          X[,index] <- y
          index <- index + 1
        }
        if(verbose) for(j in 1:nchar(str)) cat("\b")
        str <- paste(i)
        if(verbose) cat(str)
      }
      if(verbose) {
        for(i in 1:nchar(str)) cat("\b")
        for(i in 1:nchar(paste("Random Walk: Chain ", chainnum, "\nDone with: "))) cat("\b")
      }
      return(X[,(discard+1):ncol(X)])
    }, .parallel = TRUE, .paropts = list(.packages = "kmatching"))
    return(chainlist)
}
