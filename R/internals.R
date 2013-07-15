#' Creates a dummy matrix for a vector of categorical variables.
#' @keywords internal
.dummy <- function(vec) {
  names = sort(unique(vec))
  mat = matrix(rep(0, length(vec)*length(names)), nrow = length(names))
  for(i in 1:nrow(mat)) {
    mat[i,][which(vec == names[i])] = 1
  }
  return(mat)
}