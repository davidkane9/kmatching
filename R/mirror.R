#' Generates weights to a new portfolio using mirror algorithm

#' Fulfills equality constraints while maintaining randomness by
#' using a Monte Carlo Random Walk reflecting at the boundaries.

#' @param Emat This is the matrix of the equality constraint coefficients
#' @param x0 An original solution to the constraints
#' @param n Number of random solutions to output
#' 
#' @author Mike Flynn \email{<mjf2@@willaims.edu>}
#' @export
#' @examples
#' Emat = matrix(1, ncol = 3, nrow = 1)
#' x0 = c(.3, .3, .4)
#' mirror(Emat, x0, 1)
mirror <- function(Emat, x0, n, verbose = FALSE, numjump= 20) {
    Z = Null(t(Emat))
    ret = matrix(0, nrow = length(x0), ncol = n + 1)
    ## Would it be better to use apply here?
    nc = ncol(Z)
    mn = mean(x0)
    ret[, 1] = x0 + Z %*% rnorm(nc, 0, abs(mn))/sqrt(nc)
    k = 0
    bestjump = 0
    for (i in 2:(n + 1)) {
        ret[, i] = ret[, i - 1] + Z %*% rnorm(nc, 0, abs(mn))/sqrt(nc)
        while(any(ret[, i] < 0)) {
            reflection = rep(0, ncol(Emat))
            overdist = rep(0, ncol(Emat))
            overdist[which(ret[, i] < 0)] = ret[, i][which(ret[, i] < 0)]
            dist = sqrt(sum(overdist^2))
            if(verbose) str = paste("Distance from walls: ", dist, "\nBest jump: ", bestjump, sep = "" )
            if(verbose) cat(str)
            for (j in 1:ncol(Z)) {
                proj =  Z[, j] * (overdist %*% Z[, j])/Z[,j] %*% Z[, j]
                reflection = reflection  - proj
                overdist = overdist - proj
            }
            jumps = matrix(abs(rnorm(numjump, sd = 2)), nrow = 1)
            dists = apply(jumps, 2, function(x) {
              point = ret[,i] + x*reflection
              sqrt(sum(point[which(point<0)]^2))
            })
            bestjump = jumps[1, which.min(dists)]
            ret[,i] = ret[,i] + bestjump*reflection
            if(verbose) for(j in 1:nchar(str))  cat("\b")
        }
    }
    ret = ret[, 2:(n + 1)]
    return(ret)
}
