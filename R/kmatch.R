#' Create matching vectors based on data Creates a matrix of weights which
#' 
#' @param x data frame containing needed input data
#' @param weight.var character name of the column of the input weights
#' @param match.var character vector of names of columns of 'data' we wish to
#'   match on
#' @param n numeric number of weight vectors desired. Default is 1
#' @param replace logical indicating whether or not bservations weighted in the
#'   original weight.var are allowed positive weight in the output. Default is
#'   FALSE
#' @param ... parameters to be passed to the sampling methods
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#'   
#'   
#' @examples
#' set.seed(40)
#' x <- data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights <- kmatch(x, weight.var = "weight", match.var = "size", n = 100, replace = TRUE)

kmatch <- function(x, weight.var, match.var,  n = 1, replace = FALSE, ...) {
    
  stopifnot(n > 0)
  stopifnot(all(c(weight.var, match.var) %in% names(x)))
  
  equation <- constraint.equation(x, weight.var, match.var, replace)
  
  weights <- hitandrun(equation$A, equation$b, n = n, ...)
  
  ret = matrix(0, nrow = nrow(x), ncol = n)
  if(!replace) {
    ret[which(x[[weight.var]] == 0),] = weights
  } else {
    ret = weights
  }  
  return(ret)
}
