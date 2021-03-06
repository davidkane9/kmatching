#' Plots 3-d mirror output in readable format
#' 
#' Plot the steps of a random walk on the simplex in 3 varibles, using arrows to direct
#' 
#' @param m A matrix with 3 rows, corresponding to the output of a call to mirror
#' with includeInfeasible = TRUE.
#' @return A ggplot object gives the arrow plot
#' 
#' @export
#' @examples
#' A <- matrix(1, ncol = 3)
#' x0 <- c(.2, -.2, 1)
#' 
#' m <- mirror(A, x0, n = 1, includeInfeasible = TRUE)
#' plotmirror(m)
plotmirror <- function(m) {
  ## map to a traingle with vertices at (0,0), (1,0) and (1/2, sqrt(3)/2).
  ## an equilateral triangle that matches the 3-d triangle on (1,0,0), (0,1,0)
  ## and (0,0,1)
  points <- data.frame(x = m[2,] + 1/2*m[3,], y = sqrt(3)/2*m[3,])
  points$xend <- c(points$x[2:nrow(points)], NA)
  points$yend <- c(points$y[2:nrow(points)], NA)
  points <- points[-nrow(points),]
  
  p <- ggplot(data = points) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), arrow = arrow(length = unit(.2, "cm"))) +
    geom_segment(aes(x = 0, y = 0, xend = .5, yend = sqrt(3)/2)) +
    geom_segment(aes(x = .5, y = sqrt(3)/2, xend = 1, yend = 0)) + 
    geom_segment(aes(x = 1, y = 0, xend = 0, yend = 0))
  
  p
}