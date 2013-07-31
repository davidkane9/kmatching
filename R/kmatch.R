#' Create matching vectors based on data Creates a matrix of weights which
#' 
#' @param x data frame containing needed input data
#' @param weight.var character name of the column of the input weights
#' @param match.var character vector of names of columns of 'data' we wish to
#'   match on
#' @param n numeric number of weight vectors desired. Default is 1
#' @param chains number of different chains, starting from different starting points. Default is 1
#' @param replace logical indicating whether or not bservations weighted in the
#'   original weight.var are allowed positive weight in the output. Default is
#'   FALSE
#' @param ... parameters to be passed to the sampling methods
#' 
#' @return Returns a matrix of 'n' sets of weights that match the given set of
#' weights in terms of weighted averages to the 'match.var' factors. The columns
#' are the sets of weights.
#' @export
#' @author Mike Flynn \email{mflynn210@@gmail.com}
#'   
#'   
#' @examples
#' set.seed(40)
#' x <- data.frame(size = rnorm(50), weight = rep(.02, 50))
#' weights <- kmatch(x, weight.var = "weight", match.var = "size", n = 100, replace = TRUE)

kmatch <- function(x, weight.var, match.var,  n = 1, chains = 1, replace = FALSE, ...) {
  
  ## Input checking.
  
  stopifnot(is.data.frame(x))
  stopifnot(is.character(weight.var))
  stopifnot(is.character(match.var))
  stopifnot(n > 0)
  stopifnot(all(c(weight.var, match.var) %in% names(x)))
  
  ## Might consider moving the bulk of these inline comments to the description
  ## above, with better math formatting. For now, I will just leave them here.
  
  ## First step is to build the constraint equation that samples will need to 
  ## satisfy: Ax = b. Note that this equation does not include the inequality 
  ## constraint that all elements of x must be >= 0. A is a matrix which 
  ## includes the match.var variable exposures as rows. Multiplying these rows 
  ## by the column vector x equals the appropriate element in b. 
  
  ## The number of constraints, and hence the number of rows in A and elements 
  ## in b, equals the number of levels of all factor/character variables in 
  ## match.var plus the number of numeric variables plus one. The last row is 
  ## due to the requirement that the sum of the output weights must equal the 
  ## sum of the input weight.var. So, if match.var includes country (factor 
  ## variable with 10 levels), size and value (both numeric), there would be 13
  ## rows in A and elements in b.
  
  equation <- constraint.equation(x, weight.var, match.var, replace)
  
  ## Generate the weights. Are we dealing with number of chains in a sensible
  ## way? Surely, this needs to be a choice variable for kmatch.
  
  weights <- hitandrun(equation$A, equation$b, n = n, chains = chains, ...)
  
  ## Now that we have the weights, we ceate the return matrix, taking account of
  ## the value of replace.
  
  ret <- matrix(0, nrow = nrow(x), ncol = n)
  
  if(!replace){
    ret[which(x[[weight.var]] == 0),] <- weights
  } 
  else{
    ret <- weights
  }  
  
  return(ret)
}
