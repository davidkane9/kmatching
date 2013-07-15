overlap <-
function(mat, orig = NULL, compare = FALSE) {
  
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
