#' Uniformly samples from {A*x=A*x0} U {x>0}
hitandrun <- function(A, b = NULL, x0 = NULL, n, discard = 0, skiplength = 5, verbose = FALSE) {
  if(n <= 0 || n %% 1 != 0){
    stop("n must be a positive integer")
  }
  
    stopifnot(!is.null(x0) || !is.null(b))
    if(is.null(x0)) {
        str = "Finding an intial solution..."
        if(verbose) cat(str)
        ##l = lsei(E = A, F = b, G = diag(ncol(A)), H = rep(0, ncol(A)))$X
        l = t(A) %*% solve(A %*% t(A)) %*% b
        if(!(all(l > 0))) {
            if(verbose) for(i in 1:nchar(str)) cat("\b")
            str = "Using mirror algorithm to find inner solution...\n"
            if(verbose) cat(str)
        }
        y = mirror(A, l, 1, verbose)
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
            ## r is a random unit vector
            r = rnorm(ncol(Z))
            r = r/sqrt(sum(r^2))

            ## d is a unit vector in the appropriate k-plane pointing in a
            ## random direction
            u = Z%*%r
            c = y/u
            ## determine intersections of x + t*u with
            ## not sure if the following code is necessary
            if(tmin == -Inf || tmax == Inf){
              stop("problem is unbounded")
            }
            tmin = max(-c[u>0]); tmax = min(-c[u<0]);
            if(tmin==0 && tmax ==0) {
                stop("found bad direction")
            }
        }

        ##writeLines(paste("tmin: ", tmin, "\ntmax: ", tmax, "\n", sep = ""))
        ## chose a point on the line segment
        y = y + (tmin + (tmax - tmin)*runif(1))*u;
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
