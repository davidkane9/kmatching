#' Generates weights using mirror algorithm.

#' Fulfills equality constraints while maintaining randomness by
#' using a random Walk reflecting at the boundaries. Based
#' on xsample() function in limSolve package. Given a set of constraints:
#' \eqn{Ex = Ex_0, x \ge 0} mirror starts at \eqn{x_0} and repeatedly jumps from
#' the point in a random direction in the k-plane that defines \eqn{Ax=b}. It then
#' checks against \eqn{x\ge 0}. If it has violated this constraint, it projects onto 
#' the violating components and projects the resulting vector back into the plane.
#' This final vector is subtracted from the violating jump, with the length scaled by
#' a random number that is calculated to maximally reduce the distance from the walls
#' (helps it converge faster). This process is repeated until there are no components 
#' violating the constraints. In practice this process generates points in time that is
#' exponential in n, the number of components of x.
#' 
#' @param Amat This is the matrix of the equality constraint coefficients
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
#' Amat <- matrix(1, ncol = 3, nrow = 1)
#' x0 <- c(.3, .3, .4)
#' mirror(Amat, x0, 1)

mirror <- function(Amat, x0, n, verbose = FALSE, numjump= 20) {
  
    ## columns of Z are orthogonal, unit basis of null space of Amat
    ## a.k.a. vectors in the plane defined by Ax=b
    Z = Null(t(Amat))
    
    ## initialize return matrix
    ret = matrix(0, nrow = length(x0), ncol = n + 1)
    
    ## number of cols in Z and mean of x0 used to normalize jumps else
    ## the convergence time grows much faster for higher n
    nc = ncol(Z)
    mn = mean(x0)
    ## jump from initial point, distance normally distributed
    ret[, 1] = x0 + Z %*% rnorm(nc, 0, abs(mn))/sqrt(nc)
    
    ## bestjump will eventually be used to store the optimal length to scale the
    ## reflection
    bestjump = 0
    for (i in 2:(n + 1)) {
        ## jump
        ret[, i] = ret[, i - 1] + Z %*% rnorm(nc, 0, abs(mn))/sqrt(nc)
        
        ## we will compare olddist to dist, if olddist < dist, then we have
        ## moved away from feasible space with a jump, and are not converging
        olddist = Inf
        ## if any of the components is negative, mirror component back
        while(any(ret[, i] < 0)) {
            ## intialize the reflection
            reflection = rep(0, ncol(Amat))
            
            ## overdist is the vector in infeasible space
            overdist = rep(0, ncol(Amat))
            overdist[which(ret[, i] < 0)] = ret[, i][which(ret[, i] < 0)]
            ## measure distance of negative components from x ==0 
            dist = sqrt(sum(overdist^2))
            
            ## throw error if mirror not converging
            if(olddist <= dist) {
              stop("mirror failing to converge, possibly no solution")
            }
            if(verbose) str = paste("Distance from walls: ", dist, "\nBest jump: ", bestjump, sep = "" )
            if(verbose) cat(str)
            ## project the infeasible vector back into feasible space
            for (j in 1:ncol(Z)) {
                ## projection = u * (u*v)/(v*V)
                proj =  Z[, j] * (overdist %*% Z[, j])/(Z[,j] %*% Z[, j])
                ## add projection to reflection
                reflection = reflection  - proj
                ## remove component from "overdist"
                overdist = overdist - proj
            }
            ## randomly generate jump lengths, pick best, converges faster
            jumps = matrix(abs(rnorm(numjump, sd = 2)), nrow = 1)
            
            ## find distances when reflection is scaled by jumps
            dists = apply(jumps, 2, function(x) {
              point = ret[,i] + x*reflection
              sqrt(sum(point[which(point<0)]^2))
            })
            
            ## pick closest distance to zero and use that jump
            bestjump = jumps[1, which.min(dists)]
            ret[,i] = ret[,i] + bestjump*reflection
            if(verbose) for(j in 1:nchar(str))  cat("\b")
            olddist = dist
        }
    }
    ret = ret[, 2:(n + 1)]
    return(ret)
}
