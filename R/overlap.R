#' Checks the similarity of a group of vectors
#' 
#' This function takes a matrix as an argument 
#' where the columns are vectors we want to compare
#' and calculates the correlation and Euclidean
#' distance between these vectors
#' 
#' @param mat The matrix of vectors to compare
#' @param compare Whether to compare to vectors distributed uniformly on a simplex scaled
#' by the sum of the rightmost vector
#' 
#' @author Mike Flynn \email{<mflynn210@@gmail.com>}
#' @export
#' @examples
#' data <- data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights <- kmatch(data = data, match.var = "size", weight.var = "weight", n = 100, replace = TRUE)
#' overlap(weights)

overlap <- function(mat, compare = FALSE) {
  
  correlations = cor(mat)
  cat("Overlap:\n");
  cat("Mean Correlaton:  ", paste(mean(correlations[lower.tri(correlations)])), "\n")
  cat("Sd Correlations:  ", paste(sd(correlations[lower.tri(correlations)])), "\n\n")
  
  dists = apply(mat, 2, function(m) {
    return(apply(mat, 2, function(x) sqrt(sum((x-m)^2))))    
  })
  cat("Mean Distance:    ", paste(mean(dists[lower.tri(dists)])), "\n")
  cat("Sd Distance:      ", paste(sd(dists[lower.tri(dists)])), "\n\n")
  
  cat("Average weight:   ", paste(mean(mat)), "\n\n")
  
  if(compare) {
    cat("Compare to Uniform on Scaled Simplex:\n")
    matU = matrix(rexp(nrow(mat)*ncol(mat)), ncol = ncol(mat), nrow = nrow(mat))
    matU = apply(matU, 2, function(x) x/sum(x)*sum(mat[,1]))
    mat = matU
    correlations = cor(mat)
    
    cat("Mean Correlaton:  ", paste(mean(correlations[lower.tri(correlations)])), "\n")
    cat("Sd Correlations:  ", paste(sd(correlations[lower.tri(correlations)])), "\n\n")
    
    dists = apply(mat, 2, function(m) {
      return(apply(mat, 2, function(x) sqrt(sum((x-m)^2))))    
    })
    cat("Mean Distance:    ", paste(mean(dists[lower.tri(dists)])), "\n")
    cat("Sd Distance:      ", paste(sd(dists[lower.tri(dists)])), "\n\n")
    
    cat("Average weight:   ", paste(mean(mat)), "\n\n")
  }
  
}
