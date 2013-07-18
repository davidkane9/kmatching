#' Uniformly samples from {A*x=A*x0}, {x>0} or {Ax=b} & {x>0}
#' 
#' Randomly samples uniformly from a convex polytope given by linear equalities 
#' in the parameters. Uses a hit-and-run algorithm. Given constraints:
#' $$Ax = b, x\ge 0$$ or $$Ax = Ax_0, x>0$$ the algorithm finds a point on
#' the interior of the constraints. From there it picks a direction in the k-plane
#' defined by $Ax=b$, it then calculates the maximum and minimum distances it can
#' move the point in that direction, called \verb+tmin+ and \verb+tmax+. It pick a random
#' distance to travel between \verb+tmin+ and \verb+tmax+ and this is used as the next point.
#' This algorithm is useful because each sample is made in constant time. 
#' 
#' @param A Matrix of constraint coefficients, rows should correspond to each constraint. A must not 
#' have collinear rows.
#' @param b A vector of exposures that correspond to the right hand side of the constraints. Should
#' not be used at the same time as x0
#' @param x0 An original solution we want to match, should not be used at the same time as b.
#' @param n The number of output vectors desired
#' @param discard A burninlength, how many vectors should be discarded before recording
#' @param skiplength Only 1 out of every 'skiplength' vectors will be recorded, in order to
#' optimally spread out the output
#' @param verbose Give verbose output of how the function is progressing.
#' 
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

hitandrun <- function(A, b = NULL, x0 = NULL, n, discard = 0, skiplength = 5, verbose = FALSE) {
    
    if(n <= 0 || n %% 1 != 0) {
      stop("n must be a positive integer")
    }
    ## should use x0 or b, not both
    stopifnot(!is.null(x0) || !is.null(b))
    
    if(is.null(x0)) {
        str = "Finding an intial solution..."
        if(verbose) cat(str)

        ## make initial solution = Ap %*% b where 
        ## Ap is the psuedoinverse of A
        SVD = svd(A)
        d = SVD[['d']]
        V = SVD[['v']]
        U = SVD[['u']]
        ## get rid of division errors in d because they will
        ## mess up 1/d
        d[d < 1e-10] = 0
        di = 1/d
        di[di == Inf] = 0
        if(length(di) <= 1) Dt = matrix(di, ncol = 1, nrow =1) else Dt = t(diag(di))
        Ap = V %*% Dt %*% t(U)
        ## l is initial solution
        l = Ap %*% b
        
        ## if l isn't in feasible space mirror it
        if(!(all(l > 0))) {
            if(verbose) for(i in 1:nchar(str)) cat("\b")
            str = "Using mirror algorithm to find inner solution...\n"
            if(verbose) cat(str)
            y = mirror(A, l, 1, verbose)
        } else {
          y = l
        }

        if(verbose) for(i in 1:nchar(str)) cat("\b")
    } else {
        y = x0;
        if(!is.null(b)) stopifnot(A %*% x0 == b)
    }
    ## resolve weird quirk in Null() function
    if(ncol(A) ==1) {
        Z = Null(A)
    } else {
        Z = Null(t(A))
    }
    X = matrix(0, nrow = length(y), ncol = n + discard)
    if(verbose) cat("Random Walk\nDone with: ")
    str = "0"
    if(verbose) cat(str)
    index = 1
    for(i in 1:(n*skiplength+discard)) {
        tmin=0;tmax=0;
        while(tmin ==0 && tmax ==0) {
            ## r is a random unit vector in with basis in Z
            r = rnorm(ncol(Z))
            r = r/sqrt(sum(r^2))

            ## u is a unit vector in the appropriate k-plane pointing in a
            ## random direction
            u = Z%*%r
            c = y/u
            ## determine intersections of x + t*u with walls
            ## the limits on how far you can go backward and forward
            tmin = max(-c[u>0]); tmax = min(-c[u<0]);
            if(tmin == -Inf || tmax == Inf){
              stop("problem is unbounded")
            }
            if(tmin==0 && tmax ==0) {
                stop("hitandrun found can't find feasible direction, cannot generate points")
            }
        }

        ## chose a point on the line segment
        y = y + (tmin + (tmax - tmin)*runif(1))*u;
        
        ## choose a point every 'skiplength' samples
        if(i %% skiplength == 0) {
            X[,index] = y
            index = index + 1
        }
        if(verbose) for(j in 1:nchar(str)) cat("\b")
        str = paste(i)
        if(verbose) cat(str)
    }
    if(verbose) cat("\n")
    return(X[,(discard+1):ncol(X)])
}
