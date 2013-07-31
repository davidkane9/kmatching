#' Creates a dummy matrix for a vector of character/factor variables.
#' 
#' @param vec vector of character or factor variable 
#' 
#' @return Matrix of dummy variables
#' 
#' @author David Kane \email{<dave.kane@@gmail.com>}
#' @export
#' 
#' @examples
#' dummy(letters[1:3])

dummy <- function(vec){
  

  ## sort in order to match up match up exposures when printing
  names <- sort(unique(vec))
  mat <- matrix(rep(0, length(vec)*length(names)), nrow = length(names))
  for(i in 1:nrow(mat)) {
    mat[i,][which(vec == names[i])] <- 1
  }
  return(mat)
}

