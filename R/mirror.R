#' Generates weights to a new portfolio using mirror algorithm

#' Fulfills equality constraints while maintaining randomness by
#' using a Monte Carlo Random Walk reflecting at the boundaries. Based
#' on xsample() function in limSolve package.
#' 
#' @param Emat This is the matrix of the equality constraint coefficients
#' @param x0 An original solution to the constraints
#' @param n Number of random solutions to output
#' @param verbose Give verbose output describing the progress of the function
#' @param numjump The number of jumps to scatter around the direction given by the difference from zero
#' 
#' @author Mike Flynn \email{<mflynn210@@gmail.com>}
#' @export
#' 
#' @references Van Den Meershe, Karel, Karline Soetaert, and Dick Van Oevelen. "Xsample(): An R Function for Sampling Linear Inverse Problems." Journal of Statistical Software 30 (2009): 1-15. Print. \url{http://cran.cermin.lipi.go.id/web/packages/limSolve/vignettes/xsample.pdf}
#' 
#' @examples
#' Emat = matrix(1, ncol = 3, nrow = 1)
#' x0 = c(.3, .3, .4)
#' mirror(Emat, x0, 1)
mirror <- function(Emat, x0, n, verbose = FALSE, numjump= 20) {
    ## columns of Z are orthoganal, unit basis of null space of Emat
    Z = Null(t(Emat))
    ## initialize return matrix
    ret = matrix(0, nrow = length(x0), ncol = n + 1)
    nc = ncol(Z)
    mn = mean(x0)
    ## jump from initial point, distance normally distributed
    ret[, 1] = x0 + Z %*% rnorm(nc, 0, abs(mn))/sqrt(nc)
    k = 0
    bestjump = 0
    for (i in 2:(n + 1)) {
        ret[, i] = ret[, i - 1] + Z %*% rnorm(nc, 0, abs(mn))/sqrt(nc)
        olddist = Inf
        ## if any of the components is negative, mirror component back
        while(any(ret[, i] < 0)) {
            ## project vector into infeasible space
            reflection = rep(0, ncol(Emat))
            overdist = rep(0, ncol(Emat))
            overdist[which(ret[, i] < 0)] = ret[, i][which(ret[, i] < 0)]
            dist = sqrt(sum(overdist^2))
            if(olddist <= dist) {
              stop("mirror failing to converge")
            }
            if(verbose) str = paste("Distance from walls: ", dist, "\nBest jump: ", bestjump, sep = "" )
            if(verbose) cat(str)
            ## project the infeasible vector back into feasible space
            for (j in 1:ncol(Z)) {
                proj =  Z[, j] * (overdist %*% Z[, j])/Z[,j] %*% Z[, j]
                reflection = reflection  - proj
                overdist = overdist - proj
            }
            ## randomly generate jump lengths, pick best, converges faster
            jumps = matrix(abs(rnorm(numjump, sd = 2)), nrow = 1)
            dists = apply(jumps, 2, function(x) {
              point = ret[,i] + x*reflection
              sqrt(sum(point[which(point<0)]^2))
            })
            bestjump = jumps[1, which.min(dists)]
            ret[,i] = ret[,i] + bestjump*reflection
            if(verbose) for(j in 1:nchar(str))  cat("\b")
            olddist = dist
        }
    }
    ret = ret[, 2:(n + 1)]
    return(ret)
}
